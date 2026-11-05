#!/usr/bin/env python3
"""Read-only survey of candidate mtDNA control-region intervals and MITOS2 OriH.

Uses the EXACT GenBank CDS/rRNA + tRNAscan-SE tRNA feature selection of
run_genbank_mitos_blocks.build_feature_sets(). OriH comes from the separate
MITOS2 GFF, and optional D-loop annotations are read from the exported
GenBank CSV (if they are present there). No MAFFT/stitching is performed.

Reports both all ranks and the bird-like TRNP > ND6 > TRNE > TRNF path.
A species can contribute multiple candidate intervals. All locations are
in the original, unrotated FASTA coordinate system (1-based inclusive).
"""
import argparse
import collections
import csv
import re
from pathlib import Path

import chrM_block_align as base
import chrM_block_align_padded as padded
import run_genbank_mitos_blocks as pipeline
from human_genbank import load_human_genbank

RANKS=("kingdom","phylum","subphylum","superclass","class","superorder",
       "order","suborder","infraorder","parvorder")
BIRD_PATH=("TRNP","ND6","TRNE","TRNF")
ORIGIN_NAMES={"OH","ORIH","ORIGINOFREPLICATIONH","HEAVYSTRANDORIGIN"}
DLOOP_REGEX=re.compile(r"\bd[\s_-]*loop\b|\bcontrol[\s_-]*region\b",re.I)

def read_tsv_map(path):
    result={}
    with open(path,newline="",encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh,delimiter="\t"):
            key=base.norm(row.get("species",""))
            gff=(row.get("result_gff") or "").strip()
            if key and gff:
                if key in result and result[key]!=gff:
                    raise ValueError("Conflicting MITOS2 GFF paths for "+key)
                result[key]=gff
    return result

def attrs_dict(text):
    return {k:v for segment in text.split(";") if "=" in segment
            for k,v in [segment.split("=",1)]}

def is_heavy_origin(ftype,attrs):
    """Only top-level replication-origin annotations; never OH-b exons."""
    if ftype.lower() not in ("origin_of_replication","rep_origin"):
        return False
    labels=[attrs.get(k,"") for k in ("Name","gene_id","ID","product")]
    for label in labels:
        compact=re.sub(r"[^A-Z0-9]","",label.upper())
        if compact in ORIGIN_NAMES or compact in ("GENEOH","TRANSCRIPTOH"):
            return True
    return False

def read_oh_gff(path,L):
    """Convert GFF 1-based inclusive feature coordinates to modulo-zero intervals."""
    if not path or not Path(path).is_file():
        return None
    intervals=[]
    with open(path,encoding="utf-8",errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():continue
            fields=line.rstrip("\n").split("\t")
            if len(fields)<9:continue
            if not is_heavy_origin(fields[2],attrs_dict(fields[8])):continue
            try:
                start=int(fields[3]);end=int(fields[4])
            except ValueError:continue
            if not (1<=start<=L and 1<=end<=L):continue
            # Allow circular wrap only when explicitly represented by start>end.
            intervals.append((start-1,end if end>=(start) else end+L))
    return intervals

def read_dloops_csv(path):
    """Potential D-loop features in the GenBank annotation CSV export.

    Important: absence here does not mean absence in the original GenBank record.
    """
    found=collections.defaultdict(list)
    with open(path,newline="",encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            species=row.get("Organism") or row.get("organism") or row.get("species") or ""
            key=base.norm(species)
            fields=("Gene","G_Symbol","G_locus_tag","alt_name","alt_synonym",
                    "alt_Product","product","feature","type","Feature")
            descriptor=" ".join(str(row.get(k) or "") for k in fields)
            if not DLOOP_REGEX.search(descriptor):continue
            try:
                start=int(row["start"]);end=int(row["end"])
            except (ValueError,TypeError,KeyError):continue
            found[key].append((start,end,descriptor.strip()))
    return found

def circular_segments(start,end,L):
    """Half-open circle interval, endpoints may cross zero once."""
    if end<start:
        end+=L
    length=end-start
    if length<0 or length>L: return []
    left=start%L
    right=left+length
    if right<=L:return [(left,right)]
    return [(left,L),(0,right-L)]

def bp_overlap(a_start,a_end,b_start,b_end,L):
    aa=circular_segments(a_start,a_end,L)
    bb=circular_segments(b_start,b_end,L)
    return sum(max(0,min(x2,y2)-max(x1,y1))
               for x1,x2 in aa for y1,y2 in bb)

def collect_gaps(fs,L):
    """Only between consecutive *selected canonical* genes, in FASTA direction."""
    order=padded.ordered(fs,L)
    gaps=[]
    if len(order)<2:return gaps
    for idx,left in enumerate(order):
        right=order[(idx+1)%len(order)]
        f=fs[left];g=fs[right]
        start=f.end
        # unwrap the successor start relative to the current gene start
        next_start=g.start + (L if idx==len(order)-1 else 0)
        # accommodate origin-wrapping gene annotations
        if next_start<=f.start:
            next_start+=L
        gap=next_start-start
        if gap<0: # overlaps are not noncoding gaps
            continue
        gaps.append(dict(left_gene=left,right_gene=right,
                         start0=start%L,end0=next_start%L,
                         gap_bp=gap,start_unwrapped=start,end_unwrapped=next_start,
                         left_strand=f.strand,right_strand=g.strand))
    return gaps

def path_from_to(order,start,end):
    if start not in order or end not in order:return []
    a=order.index(start);b=order.index(end)
    return [order[(a+i)%len(order)] for i in range((b-a)%len(order)+1)]

def write_tsv(path,rows,fields):
    with open(path,"w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",extrasaction="ignore")
        w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ranks",nargs="+",default=list(RANKS),choices=RANKS)
    p.add_argument("--fasta-dir",default="data/unrot_fasta")
    p.add_argument("--genbank-annotations",default="data/kingdom_genes_v2.csv")
    p.add_argument("--gene-aliases",default="data/genefiltfinal5_2026_09_24.csv")
    p.add_argument("--trnascan-csv",default="data/consolidated_trnascan-2.csv")
    p.add_argument("--mitos2-summary",default="data/mitos2_species_summary.tsv")
    p.add_argument("--human-fasta",default="data/Humans_mtdna.fasta")
    p.add_argument("--human-genbank",default="data/sequence.gb")
    p.add_argument("--min-gap-bp",type=int,default=300,
                   help="Exploratory threshold for large noncoding intervals; not a CR definition")
    p.add_argument("--output-dir",default="control_region_candidate_qc")
    a=p.parse_args()
    if a.min_gap_bp<0:p.error("--min-gap-bp must be nonnegative")
    out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    human=list(base.fasta(a.human_fasta))
    if len(human)!=1:raise ValueError("Expected one human reference FASTA")
    refs=load_human_genbank(a.human_genbank,human[0][1],base.Feature,base.EXPECTED)
    ref_order=[(g,refs[g].strand) for g in padded.ordered(refs,len(human[0][1]))]
    aliases=pipeline.load_aliases(a.gene_aliases)
    gb,_,_,_=pipeline.load_genbank_csv(a.genbank_annotations,aliases,ref_order)
    trnas=base.read_trnas(a.trnascan_csv)
    gff_map=read_tsv_map(a.mitos2_summary)
    dloops=read_dloops_csv(a.genbank_annotations)
    origin_cache={}
    species_rows=[];gap_rows=[];summary_rows=[]
    for rank in a.ranks:
        print("Processing",rank,flush=True)
        fasta_path=Path(a.fasta_dir)/f"{rank}_chrm.faa"
        selected,_,_=pipeline.build_feature_sets(rank,fasta_path,gb,trnas)
        per_rank=[]
        seen=set()
        for name,seq in base.fasta(fasta_path):
            key=base.norm(name)
            if key=="homosapiens":continue
            if key in seen:raise ValueError(f"Duplicate species in {rank}: {key}")
            seen.add(key);L=len(seq);fs=selected.get(key,{})
            order=padded.ordered(fs,L)
            path=path_from_to(order,"TRNP","TRNF")
            bird=path==list(BIRD_PATH)
            direct=path==["TRNP","TRNF"]
            reverse=path_from_to(order,"TRNF","TRNP")==["TRNF","TRNP"]
            if key not in origin_cache:
                origin_cache[key]=read_oh_gff(gff_map.get(key),L)
            oh=origin_cache[key]
            gb_dloop=[]
            for begin,end,desc in dloops.get(key,[]):
                if 1<=begin<=L and 1<=end<=L:
                    gb_dloop.append((begin-1,end if end>=begin else end+L))
            gaps=collect_gaps(fs,L)
            longest=max((g["gap_bp"] for g in gaps),default=0)
            large=[g for g in gaps if g["gap_bp"]>=a.min_gap_bp]
            oh_hits=0;dloop_hits=0
            bird_gap=None
            for ix,g in enumerate(gaps,1):
                overlaps_oh=sum(bp_overlap(g["start_unwrapped"],g["end_unwrapped"],
                                          start,end,L)>0 for start,end in (oh or []))
                overlaps_dloop=sum(bp_overlap(g["start_unwrapped"],g["end_unwrapped"],
                                             start,end,L)>0 for start,end in gb_dloop)
                oh_hits+=int(overlaps_oh>0)
                dloop_hits+=int(overlaps_dloop>0)
                if g["left_gene"]=="TRNE" and g["right_gene"]=="TRNF":
                    bird_gap=(g,overlaps_oh,overlaps_dloop)
                gap_rows.append(dict(rank=rank,species=name,interval_id=ix,
                    left_gene=g["left_gene"],right_gene=g["right_gene"],
                    start_1based=(g["start0"]+1),end_1based=g["end0"],
                    gap_bp=g["gap_bp"],large_gap=int(g["gap_bp"]>=a.min_gap_bp),
                    left_strand=g["left_strand"],right_strand=g["right_strand"],
                    n_OriH_overlaps=overlaps_oh,n_GenBank_Dloop_overlaps=overlaps_dloop,
                    is_bird_TRNE_TRNF_gap=int(bird and g["left_gene"]=="TRNE" and g["right_gene"]=="TRNF")))
            row=dict(rank=rank,species=name,genome_bp=L,n_selected_genes=len(fs),
                trnp_present=int("TRNP" in fs),trnf_present=int("TRNF" in fs),
                human_like_forward=int(direct and fs["TRNP"].strand==refs["TRNP"].strand
                        and fs["TRNF"].strand==refs["TRNF"].strand),
                reverse_adjacent=int(reverse),
                trnp_to_trnf_path=">".join(path),
                bird_like_path=int(bird),
                OriH_gff_status="missing_gff" if oh is None else "read",
                n_OriH_features="" if oh is None else len(oh),
                n_GenBank_Dloop_export_features=len(gb_dloop),
                n_intergenic_gaps=len(gaps),largest_gap_bp=longest,
                n_gaps_ge_threshold=len(large),
                n_gaps_overlapping_OriH=oh_hits,
                n_gaps_overlapping_GenBank_Dloop=dloop_hits,
                bird_TRNE_TRNF_gap_bp=bird_gap[0]["gap_bp"] if bird_gap else "",
                bird_TRNE_TRNF_overlaps_OriH=int(bird_gap[1]>0) if bird_gap else "",
                bird_TRNE_TRNF_overlaps_Dloop=int(bird_gap[2]>0) if bird_gap else "",
                candidate_status=("unique_large_gap_with_OriH" if sum(
                    g["gap_bp"]>=a.min_gap_bp and any(bp_overlap(
                    g["start_unwrapped"],g["end_unwrapped"],x,y,L)>0
                    for x,y in (oh or [])) for g in gaps)==1 else
                    "no_large_gap" if not large else "multiple_or_unvalidated"))
            per_rank.append(row);species_rows.append(row)
        n=len(per_rank);bird=[r for r in per_rank if r["bird_like_path"]]
        summary_rows.append(dict(rank=rank,n_species=n,
            n_human_like_forward=sum(r["human_like_forward"] for r in per_rank),
            n_bird_like_path=len(bird),
            n_bird_with_TRNE_TRNF_gap=sum(r["bird_TRNE_TRNF_gap_bp"]!="" for r in bird),
            n_bird_with_large_TRNE_TRNF_gap=sum(
                isinstance(r["bird_TRNE_TRNF_gap_bp"],int) and
                r["bird_TRNE_TRNF_gap_bp"]>=a.min_gap_bp for r in bird),
            n_bird_with_OriH_gff=sum(r["OriH_gff_status"]=="read" and
                r["n_OriH_features"]>0 for r in bird),
            n_bird_TRNE_TRNF_gap_overlaps_OriH=sum(
                r["bird_TRNE_TRNF_overlaps_OriH"]==1 for r in bird),
            n_bird_TRNE_TRNF_gap_overlaps_Dloop=sum(
                r["bird_TRNE_TRNF_overlaps_Dloop"]==1 for r in bird),
            n_any_large_gap=sum(r["n_gaps_ge_threshold"]>0 for r in per_rank),
            n_any_gap_overlaps_OriH=sum(r["n_gaps_overlapping_OriH"]>0 for r in per_rank),
            n_any_gap_overlaps_Dloop=sum(r["n_gaps_overlapping_GenBank_Dloop"]>0 for r in per_rank)))
    write_tsv(out/"all_ranks_summary.tsv",summary_rows,list(summary_rows[0]))
    write_tsv(out/"all_ranks_species.tsv",species_rows,list(species_rows[0]))
    write_tsv(out/"all_ranks_intergenic_intervals.tsv",gap_rows,
        ["rank","species","interval_id","left_gene","right_gene",
         "start_1based","end_1based","gap_bp","large_gap",
         "left_strand","right_strand","n_OriH_overlaps",
         "n_GenBank_Dloop_overlaps","is_bird_TRNE_TRNF_gap"])
    write_tsv(out/"bird_like_TRNP_ND6_TRNE_TRNF_species.tsv",
        [r for r in species_rows if r["bird_like_path"]],list(species_rows[0]))
    for row in summary_rows:
        print(row["rank"],"bird-like",row["n_bird_like_path"],
              "TRNE-TRNF + OriH",row["n_bird_TRNE_TRNF_gap_overlaps_OriH"],
              "other OriH-supported gap",row["n_any_gap_overlaps_OriH"],flush=True)
    print("Results written to",out.resolve())

if __name__=="__main__":
    main()
