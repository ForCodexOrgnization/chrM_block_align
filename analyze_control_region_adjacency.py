#!/usr/bin/env python3
"""QC of circular TRNP/TRNF gene adjacency in the same input species as the block pipeline.

Reports:
  - all-rank summary of annotation completeness, circular adjacency and control-region eligibility
  - a per-species table including the intervening *annotated canonical* genes and exclusion reasons

"Adjacent" means no other RESOLVED/ANNOTATED one of the canonical 37 mitochondrial
genes between the two genes. It does NOT establish physical/biological absence of
unannotated homologs or identify the true control region by sequence homology.

Use --pipeline genbank for the GenBank-first + MITOS2-fallback feature selection
used by run_genbank_mitos_blocks.py, or --pipeline mitos2 for the older
MITOS2 + tRNAscan-SE annotation selection of chrM_block_align_padded.py.

No MAFFT or existing stitched results are needed.
"""
import argparse
import csv
from pathlib import Path
from collections import Counter

import chrM_block_align as base
import chrM_block_align_padded as padded
import run_genbank_mitos_blocks as gb_pipeline
from human_genbank import load_human_genbank

RANKS=("kingdom","phylum","subphylum","superclass","class",
       "superorder","order","suborder","infraorder","parvorder")
TARGET=("TRNP","TRNF")

def parse_args():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--pipeline",choices=["genbank","mitos2"],default="genbank")
    p.add_argument("--fasta-dir",default="data/unrot_fasta")
    p.add_argument("--human-fasta",default="data/Humans_mtdna.fasta")
    p.add_argument("--human-genbank",default="data/sequence.gb")
    p.add_argument("--genbank-annotations",default="data/kingdom_genes_v2.csv")
    p.add_argument("--gene-aliases",default="data/genefiltfinal5_2026_09_24.csv")
    p.add_argument("--mitos2-summary",default="data/mitos2_species_summary.tsv")
    p.add_argument("--trnascan-csv",default="data/consolidated_trnascan-2.csv")
    p.add_argument("--max-control-bp",type=int,default=2500,
                   help="Replicate padded pipeline's maximum intergenic interval length")
    p.add_argument("--ranks",nargs="+",choices=RANKS,default=list(RANKS))
    p.add_argument("--output-dir",default=None)
    return p.parse_args()

def load_gff_summary(path):
    with open(path,newline="") as handle:
        return {base.norm(r["species"]):r["result_gff"]
                for r in csv.DictReader(handle,delimiter="\t")
                if r.get("species") and r.get("result_gff")}

def features_for_rank(rank,records,a,genbank=None,mitos=None,trnas=None):
    if a.pipeline=="genbank":
        # Same feature-merging routine as the GenBank-first stitcher. It picks
        # GenBank annotations first, then MITOS2 for unresolved gene identities,
        # and intentionally does not use tRNAscan-SE.
        fasta_path=Path(a.fasta_dir)/f"{rank}_chrm.faa"
        selected,_,_=gb_pipeline.build_feature_sets(rank,fasta_path,genbank,mitos)
        return selected
    selected={}
    for name,seq in records:
        key=base.norm(name)
        if key=="homosapiens":
            continue
        gff=mitos.get(key)
        if not gff or not Path(gff).is_file():
            selected[key]={}
            continue
        fs,_=base.pick_features(base.read_gff(gff,len(seq))+trnas.get(key,[]),len(seq))
        selected[key]=fs
    return selected

def qc_one(rank,name,seq,fs,human_ref,max_control_bp):
    L=len(seq)
    order=padded.ordered(fs,L)
    present_p="TRNP" in fs
    present_f="TRNF" in fs
    row=dict(rank=rank,species=name,sequence_length=L,n_annotated_canonical_genes=len(order),
             trnp_present=int(present_p),trnf_present=int(present_f),
             both_present=int(present_p and present_f),
             trnp_trnf_forward_adjacent=0,either_direction_adjacent=0,
             n_intervening_annotated_genes="",intervening_annotated_genes="",
             trnp_strand=fs["TRNP"].strand if present_p else "",
             trnf_strand=fs["TRNF"].strand if present_f else "",
             strand_matches_human=0,control_interval_bp="",
             within_max_control_bp=0,control_alignment_eligible=0,
             classification="",note="")
    if not (present_p and present_f):
        if not present_p and not present_f:status="missing_both"
        elif not present_p:status="missing_trnp"
        else:status="missing_trnf"
        row["classification"]=status
        return row
    if len(order)<2:
        row["classification"]="insufficient_resolved_genes"
        return row
    i=order.index("TRNP")
    j=order.index("TRNF")
    # Circular forward direction is that used by the MAFFT control-region
    # extraction: TRNP -> TRNF in ascending (unrotated) genome coordinates.
    intervening=[order[(i+k)%len(order)] for k in range(1,(j-i)%len(order))]
    forward=(j==(i+1)%len(order))
    reverse=(i==(j+1)%len(order))
    row["trnp_trnf_forward_adjacent"]=int(forward)
    row["either_direction_adjacent"]=int(forward or reverse)
    row["n_intervening_annotated_genes"]=len(intervening)
    row["intervening_annotated_genes"]=",".join(intervening)
    strand=(fs["TRNP"].strand==human_ref["TRNP"].strand
            and fs["TRNF"].strand==human_ref["TRNF"].strand)
    row["strand_matches_human"]=int(strand)
    gap=(fs["TRNF"].start-fs["TRNP"].end)%L
    row["control_interval_bp"]=gap
    row["within_max_control_bp"]=int(gap<=max_control_bp)
    eligible=forward and strand and gap<=max_control_bp
    row["control_alignment_eligible"]=int(eligible)
    if not forward:
        row["classification"]="reverse_only" if reverse else "intervening_annotated_genes"
    elif not strand:
        row["classification"]="strand_mismatch"
    elif gap>max_control_bp:
        row["classification"]="interval_too_long"
    else:
        row["classification"]="eligible"
    # Some gene annotations may be missing: observed adjacency may not be
    # true biological adjacency. The two properties are not independent.
    if len(order)<37:
        row["note"]="Observed adjacency is conditional on resolved gene annotations"
    return row

def write_tsv(path,rows,fields):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="") as handle:
        out=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",extrasaction="ignore")
        out.writeheader()
        out.writerows(rows)

def summarize(rank,rows):
    n=len(rows)
    count=lambda key:sum(int(r[key]) for r in rows)
    both=count("both_present")
    forward=count("trnp_trnf_forward_adjacent")
    any_adj=count("either_direction_adjacent")
    eligible=count("control_alignment_eligible")
    cls=Counter(r["classification"] for r in rows)
    return dict(rank=rank,n_species=n,
        n_with_trnp=count("trnp_present"),n_with_trnf=count("trnf_present"),
        n_with_both=both,pct_with_both_of_all=round(100*both/n,2) if n else "",
        n_forward_adjacent_no_intervening_annotated_genes=forward,
        pct_forward_adjacent_of_both=round(100*forward/both,2) if both else "",
        pct_forward_adjacent_of_all=round(100*forward/n,2) if n else "",
        n_adjacent_either_direction=any_adj,
        n_with_intervening_annotated_genes=cls["intervening_annotated_genes"],
        n_reverse_only=cls["reverse_only"],
        n_forward_adjacent_strand_mismatch=cls["strand_mismatch"],
        n_forward_adjacent_interval_too_long=cls["interval_too_long"],
        n_control_alignment_eligible=eligible,
        pct_eligible_of_both=round(100*eligible/both,2) if both else "",
        pct_eligible_of_all=round(100*eligible/n,2) if n else "",
        n_missing_trnp=cls["missing_trnp"],n_missing_trnf=cls["missing_trnf"],
        n_missing_both=cls["missing_both"])

def main():
    a=parse_args()
    if a.max_control_bp<0:raise ValueError("--max-control-bp must not be negative")
    out=Path(a.output_dir or f"control_region_adjacency_qc_{a.pipeline}")
    out.mkdir(parents=True,exist_ok=True)
    human_records=list(base.fasta(a.human_fasta))
    if len(human_records)!=1 or len(human_records[0][1])!=16569:
        raise ValueError("Human FASTA must be exactly one 16569-bp sequence")
    human_features=load_human_genbank(a.human_genbank,human_records[0][1],
                                     base.Feature,base.EXPECTED)
    mitos=load_gff_summary(a.mitos2_summary)
    if a.pipeline=="genbank":
        aliases=gb_pipeline.load_aliases(a.gene_aliases)
        reforder=[(gene,human_features[gene].strand)
                  for gene in padded.ordered(human_features,16569)]
        gb,_,_,_=gb_pipeline.load_genbank_csv(a.genbank_annotations,aliases,reforder)
        trnas=None
    else:
        gb=None
        trnas=base.read_trnas(a.trnascan_csv)
    details=[]
    summaries=[]
    for rank in a.ranks:
        fa=Path(a.fasta_dir)/f"{rank}_chrm.faa"
        records=list(base.fasta(fa))
        seen=set()
        for name,_ in records:
            key=base.norm(name)
            if key in seen:raise ValueError(f"Duplicate normalized species name in {fa}: {name}")
            seen.add(key)
        selected=features_for_rank(rank,records,a,gb,mitos,trnas)
        per_rank=[]
        for name,seq in records:
            key=base.norm(name)
            if key=="homosapiens":
                continue
            # All nonhuman FASTA records are in denominator, including those
            # without GFF or resolved gene annotations.
            row=qc_one(rank,name,seq,selected.get(key,{}),human_features,a.max_control_bp)
            per_rank.append(row)
        if not per_rank:raise ValueError(f"No nonhuman species in {fa}")
        fields=list(per_rank[0])
        write_tsv(out/f"{rank}_species_adjacency.tsv",per_rank,fields)
        details.extend(per_rank)
        summary=summarize(rank,per_rank)
        summaries.append(summary)
        print(f"{rank:12s} {summary['n_forward_adjacent_no_intervening_annotated_genes']:5d}/"
              f"{summary['n_with_both']:5d} forward-adjacent with both genes; "
              f"{summary['n_control_alignment_eligible']:5d}/{summary['n_species']:5d} "
              f"eligible among all species")
    write_tsv(out/"all_ranks_species_adjacency.tsv",details,list(details[0]))
    write_tsv(out/"all_ranks_adjacency_summary.tsv",summaries,list(summaries[0]))
    print("Wrote",out.resolve())

if __name__=="__main__":
    main()
