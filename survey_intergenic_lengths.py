#!/usr/bin/env python3
"""Survey physical intergenic distances before choosing padding for chrM blocks.

Uses existing rank FASTAs, MITOS2 CDS/rRNA annotation and tRNAscan tRNA
annotation, with the canonical human GenBank genes as reference gene order.
This is annotation-based descriptive QC, not an alignment or evidence of
homology. Lengths are measured only for adjacent, unique annotated features.
For a pair (A,B) in human order, counts other-species A->B in its physical
circular gene order; reverse/rearranged/nonadjacent pairs are flagged.

Positive length = unannotated interval; negative length = gene overlap.
Coordinates are circular and use each species' own chrM length.
"""
import argparse
import collections
import csv
from pathlib import Path
import numpy as np
import pandas as pd
from chrM_block_align import (Feature, EXPECTED, fasta, norm, read_gff,
                              read_trnas, pick_features)
from human_genbank import load_human_genbank

def physical_pair_distance(a,b,L):
    """Signed distance from end of a to start of b for circular adjacent genes.

    Inputs features are stored in 0-based forward-genome coordinate order.
    Return negative for overlap and positive for gaps. The wrap-around pair
    (last feature, first feature) has b.start + L as the downstream start.
    """
    right=b.start + (L if b.start<=a.start else 0)
    return right-a.end

def human_pairs(features,L):
    order=sorted(features,key=lambda g:features[g].start%L)
    return [(order[i-1],order[i]) for i in range(len(order))]

def survey(samples,reference,rank):
    """Return raw observations for every human-adjacent gene pair and species."""
    pairs=human_pairs(reference,16569)
    records=[]
    for species,(sequence,features) in samples.items():
        L=len(sequence)
        ordered=sorted(features.values(),key=lambda f:f.start%L)
        idx={f.gene:i for i,f in enumerate(ordered)}
        for left,right in pairs:
            a=features.get(left); b=features.get(right)
            if a is None or b is None:
                status="missing_annotation"
                gap=None
            elif (idx[left]+1)%len(ordered)!=idx[right]:
                status="not_adjacent"
                gap=None
            elif a.strand!=reference[left].strand or b.strand!=reference[right].strand:
                status="strand_mismatch"
                gap=None
            else:
                gap=physical_pair_distance(a,b,L)
                # Very long genes/annotation anomalies can yield negative
                # distances even for circular-adjacent features.
                status=("overlap" if gap<0 else "touching" if gap==0 else "gap")
            records.append({"rank":rank,"species":species,"upstream_gene":left,
                            "downstream_gene":right,"status":status,
                            "intergenic_bp":gap,"genome_bp":L})
    return pd.DataFrame(records)

def summarize(observations):
    rows=[]
    for (rank,left,right),g in observations.groupby(
            ["rank","upstream_gene","downstream_gene"],sort=False):
        positive=g.loc[g.status=="gap","intergenic_bp"]
        signed=g.loc[g.intergenic_bp.notna(),"intergenic_bp"]
        entry={"rank":rank,"upstream_gene":left,"downstream_gene":right,
               "n_species":len(g),"n_adjacent":int(g.status.isin(["gap","touching","overlap"]).sum()),
               "n_missing":int((g.status=="missing_annotation").sum()),
               "n_not_adjacent":int((g.status=="not_adjacent").sum()),
               "n_strand_mismatch":int((g.status=="strand_mismatch").sum()),
               "n_gap":int((g.status=="gap").sum()),
               "n_touching":int((g.status=="touching").sum()),
               "n_overlap":int((g.status=="overlap").sum()),
               "adjacent_fraction":float(g.status.isin(["gap","touching","overlap"]).mean()),
               "median_gap_bp":positive.median(),"p90_gap_bp":positive.quantile(.9),
               "p95_gap_bp":positive.quantile(.95),"p99_gap_bp":positive.quantile(.99),
               "max_gap_bp":positive.max(),
               "median_signed_bp":signed.median()}
        for padding in (50,100):
            # How often would the non-overlapping intergenic interval be
            # fully spanned by either flank (padding) or both (2*padding)?
            entry[f"gap_le_{padding}_bp_fraction"]=float((positive<=padding).mean()) if len(positive) else float("nan")
            entry[f"gap_le_{padding*2}_bp_fraction"]=float((positive<=padding*2).mean()) if len(positive) else float("nan")
        rows.append(entry)
    return pd.DataFrame(rows)

def load_samples(fasta_path,gffs,trnas,reference):
    samples={};skipped=[]
    for name,seq in fasta(fasta_path):
        key=norm(name)
        if key in samples:
            skipped.append({"species":name,"reason":"duplicate_species"});continue
        if key=="homosapiens":
            if len(seq)!=16569:raise ValueError("Human FASTA must be 16569 bp")
            features=reference
        else:
            gf=gffs.get(key)
            if not gf or not Path(gf).exists():
                skipped.append({"species":name,"reason":"missing_gff"});continue
            features,amb=pick_features(read_gff(gf,len(seq))+trnas.get(key,[]),len(seq))
        samples[key]=(seq,features)
    if "homosapiens" not in samples:raise ValueError("Human missing in rank FASTA")
    return samples,pd.DataFrame(skipped,columns=["species","reason"])

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta-dir",default="data/unrot_fasta")
    p.add_argument("--annotation-summary",default="data/mitos2_species_summary.tsv")
    p.add_argument("--trnascan-csv",default="data/consolidated_trnascan-2.csv")
    p.add_argument("--human-fasta",default="data/Humans_mtdna.fasta")
    p.add_argument("--human-genbank",default="data/sequence.gb")
    p.add_argument("--ranks",nargs="+",default=["class"])
    p.add_argument("--out-dir",default="block_results_v6/intergenic_length_survey")
    args=p.parse_args()
    out=Path(args.out_dir);out.mkdir(parents=True,exist_ok=True)
    human_records=list(fasta(args.human_fasta))
    if len(human_records)!=1 or len(human_records[0][1])!=16569:
        raise ValueError("Need a single 16569-bp human reference")
    reference=load_human_genbank(args.human_genbank,human_records[0][1],Feature,EXPECTED)
    with open(args.annotation_summary) as inp:
        gffs={norm(r["species"]):r["result_gff"] for r in csv.DictReader(inp,delimiter="\t") if r.get("result_gff")}
    trnas=read_trnas(args.trnascan_csv)
    observed=[]
    for rank in args.ranks:
        path=Path(args.fasta_dir)/f"{rank}_chrm.faa"
        if not path.exists():raise FileNotFoundError(path)
        samples,skipped=load_samples(path,gffs,trnas,reference)
        df=survey(samples,reference,rank)
        df.to_csv(out/f"{rank}_intergenic_observations.tsv",sep="\t",index=False)
        summary=summarize(df)
        summary.to_csv(out/f"{rank}_intergenic_summary.tsv",sep="\t",index=False)
        skipped.to_csv(out/f"{rank}_skipped_species.tsv",sep="\t",index=False)
        observed.append(df)
        gaps=df[df.status=="gap"].intergenic_bp
        print(f"{rank}: {len(samples)} species; {len(gaps)} positive adjacency gaps; "
              f"median={gaps.median() if len(gaps) else 'NA'}, "
              f"p95={gaps.quantile(.95) if len(gaps) else 'NA'}, "
              f"max={gaps.max() if len(gaps) else 'NA'} bp")
    if observed:
        all_obs=pd.concat(observed,ignore_index=True)
        summarize(all_obs).to_csv(out/"all_ranks_intergenic_summary.tsv",sep="\t",index=False)
    print("Reports:",out)

if __name__=="__main__":main()
