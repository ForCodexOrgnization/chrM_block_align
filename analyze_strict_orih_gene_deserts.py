#!/usr/bin/env python3
"""Strict OriH-contained gene-desert survey, one interval per species per mode.

A: no annotated human-homologous mitochondrial genes (37 canonical genes).
B: no annotated genes of any kind (GenBank, MITOS2, tRNAscan-SE).

For each species and mode, enumerate ALL nonzero gaps in the union of gene
annotations on the circular FASTA; retain only gaps fully containing every
MITOS2 OH origin_of_replication feature; choose the longest surviving gap,
tie-breaking by starting coordinate. No length cutoff or tRNA flanking rule.
A species with missing OH annotations is not a success. The interval may still
contain non-gene regulatory elements, and genes missing from the annotations
cannot be excluded. OriH-b exons are not independently counted.
"""
import argparse
import collections
import csv
from pathlib import Path
import re

import chrM_block_align as base
from analyze_longest_all_gene_desert import (
    RANKS, TYPES, fasta, norm, number, spans, union, overlap,
    gff_paths, trnascan_genes, write,
)
import run_genbank_mitos_blocks as pipeline
import resolve_gene_order_annotations as go

OH_LABELS={"OH","ORIH","GENEOH","TRANSCRIPTOH","HEAVYSTRANDORIGIN"}
# Explicit names for MITOS2 non-gene labels which can have type 'gene'.
SKIP_LABELS={"OH","ORIH","OL","ORIL","GENEOH","GENEOL",
             "TRANSCRIPTOH","TRANSCRIPTOOL","DLOOP","CONTROLREGION"}
CANONICAL=base.EXPECTED

def cleaned(value):
    return re.sub(r"[^A-Z0-9]","",str(value or "").upper())

def read_genbank(path, aliases):
    """Retain all GenBank gene coordinates; tag canonical rows for mode A.

    A row is considered human-homologous if ANY of its alias fields resolves
    to a canonical gene, even when a gene order ambiguity remains unresolved.
    """
    result=collections.defaultdict(list)
    with open(path,newline="",encoding="utf-8-sig") as fh:
        reader=csv.DictReader(fh)
        for row in reader:
            key=norm(row.get("Organism") or row.get("organism") or row.get("species") or "")
            if not key:continue
            a=number(row.get("start"));b=number(row.get("end"))
            if a is None or b is None:continue
            type_=str(row.get("feature") or row.get("Feature") or row.get("type") or "").strip().lower()
            # If an explicit type exists, discard non-gene annotations.
            if type_ and type_ not in TYPES:continue
            labels=[row.get(k,"") for k in ("Gene","G_Symbol","G_locus_tag","alt_name","alt_Product")]
            if any(cleaned(v) in SKIP_LABELS for v in labels):continue
            is_canonical=bool(go.candidates(row,aliases)&CANONICAL)
            result[key].append((a,b,is_canonical))
    return result

def canonical_mitos_trna(attrs):
    # Mitochondrial gene labels usually use trnA / trnL1 / trnS2.
    def decode(value):
        s=str(value or "")
        s=re.sub(r"^(gene_|transcript_|mt[_-])","",s,flags=re.I)
        s=s.split("(")[0]
        s=re.sub(r"_[0-9]+$","",s)
        s=cleaned(s)
        if s in CANONICAL and s in base.TRNA:return s
        if re.fullmatch(r"TRN[ARNDCQEGHIKMFPSTWYV]",s):return s
        return None
    for field in ("name","gene_id","id"):
        gene=decode(attrs.get(field,""))
        if gene:return gene
    return None

def mitos2_features(path,L):
    """Return all gene intervals, human-homologous intervals, and OH loci.

    Only GFF top-level origins are counted. GFF 'gene' annotations identified
    as OH/OL or D-loop are excluded from the occupied gene intervals.
    """
    all_genes=[];canonical_genes=[];origins=[]
    with open(path,encoding="utf-8",errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):continue
            parts=line.rstrip("\n").split("\t")
            if len(parts)<9:continue
            ft=parts[2].lower()
            attrs={k.casefold():v for item in parts[8].split(";") if "=" in item
                   for k,v in [item.split("=",1)]}
            labels=[cleaned(attrs.get(k,"")) for k in ("name","gene_id","id","product")]
            pieces=spans(number(parts[3]),number(parts[4]),L)
            if not pieces:continue
            if ft in ("origin_of_replication","rep_origin"):
                if any(label in OH_LABELS for label in labels):
                    origins.append(pieces)
                continue
            if ft not in TYPES or any(label in SKIP_LABELS for label in labels):
                continue
            all_genes.extend(pieces)
            candidate=None
            for k in ("name","gene_id","id"):
                candidate=base.canonical_mitos(attrs.get(k,""))
                if candidate:break
            if candidate is None:
                candidate=canonical_mitos_trna(attrs)
            if candidate in CANONICAL:
                canonical_genes.extend(pieces)
    return all_genes,canonical_genes,origins

def circular_gaps(occupied,L):
    """All positive circular gaps, represented as unwrapped [start,end)."""
    blocks=union(occupied)
    if not blocks:return []
    result=[]
    for i,(_,right) in enumerate(blocks):
        next_left=blocks[(i+1)%len(blocks)][0]+(L if i==len(blocks)-1 else 0)
        if next_left>right:result.append((right,next_left))
    return result

def origin_union(origins):
    return union([piece for group in origins for piece in group])

def select_interval(occupied,origins,L):
    """Require 100% containment of ALL distinct annotated OH bases."""
    gaps=circular_gaps(occupied,L)
    unique_origin=origin_union(origins)
    origin_bp=sum(b-a for a,b in unique_origin)
    if not origin_bp:return gaps,None,0,0,0
    valid=[]
    for start,end in gaps:
        covered=sum(overlap(start,end,a,b,L) for a,b in unique_origin)
        if covered==origin_bp:
            valid.append((start,end))
    best=min(valid,key=lambda x:(-(x[1]-x[0]),x[0]%L)) if valid else None
    return gaps,best,len(valid),origin_bp,len(unique_origin)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta-dir",default="data/unrot_fasta")
    p.add_argument("--genbank-annotations",default="data/kingdom_genes_v2.csv")
    p.add_argument("--gene-aliases",default="data/genefiltfinal5_2026_09_24.csv")
    p.add_argument("--trnascan-csv",default="data/consolidated_trnascan-2.csv")
    p.add_argument("--mitos2-summary",default="data/mitos2_species_summary.tsv")
    p.add_argument("--ranks",nargs="+",choices=RANKS,default=list(RANKS))
    p.add_argument("--out-dir",default="strict_orih_gene_desert_qc")
    a=p.parse_args();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    aliases=pipeline.load_aliases(a.gene_aliases)
    gb=read_genbank(a.genbank_annotations,aliases)
    tr=trnascan_genes(a.trnascan_csv)
    mitos_files=gff_paths(a.mitos2_summary)
    cache={};rows=[];summary=[]
    for rank in a.ranks:
        by_rank=collections.defaultdict(list);seen=set()
        for name,seq in fasta(Path(a.fasta_dir)/(rank+"_chrm.faa")):
            key=norm(name)
            if key=="homosapiens":continue
            if key in seen:raise ValueError(f"Duplicate species in {rank}: {name}")
            seen.add(key);L=len(seq)
            if (key,L) not in cache:
                gff=mitos_files.get(key)
                cache[key,L]=mitos2_features(gff,L) if gff and Path(gff).is_file() else ([],[],[])
            m_all,m_canon,origins=cache[key,L]
            gb_all=[];gb_canon=[]
            for first,last,is_canonical in gb.get(key,[]):
                pieces=spans(first,last,L)
                gb_all.extend(pieces)
                if is_canonical:gb_canon.extend(pieces)
            tRNA=[]
            for first,last in tr.get(key,[]):
                tRNA.extend(spans(first,last,L))
            has_gff=bool(mitos_files.get(key) and Path(mitos_files[key]).is_file())
            for mode,occupied in (
                ("no_human_homologous_genes",gb_canon+m_canon+tRNA),
                ("no_genes_at_all",gb_all+m_all+tRNA),
            ):
                gaps,winner,n_valid,origin_bp,n_oh_parts=select_interval(occupied,origins,L)
                best=max((e-s for s,e in gaps),default=0)
                row=dict(rank=rank,species=name,mode=mode,genome_bp=L,
                    n_genbank_features=len(gb.get(key,[])),
                    n_mitos2_gene_interval_pieces=len(m_all),
                    n_trnascan_features=len(tr.get(key,[])),
                    mitos2_gff_available=int(has_gff),
                    n_OriH_features=len(origins),n_OriH_union_parts=n_oh_parts,
                    OriH_annotated_bp=origin_bp,n_gene_free_gaps=len(gaps),
                    largest_unfiltered_gap_bp=best,
                    n_fully_containing_gaps=n_valid,
                    selected=int(winner is not None),
                    selected_gap_bp=winner[1]-winner[0] if winner else "",
                    selected_start_1based=(winner[0]%L)+1 if winner else "",
                    selected_end_1based=winner[1]%L if winner else "",
                    selected_wraps_origin=int(winner[1]>L) if winner else "",
                    status=("missing_mitos2_gff" if not has_gff else
                        "no_OriH_prediction" if not origins else
                        "no_annotated_genes" if not occupied else
                        "no_gene_free_gap" if not gaps else
                        "no_full_OriH_containment" if not winner else
                        "selected_100pct_containment"))
                rows.append(row);by_rank[mode].append(row)
        for mode,group in by_rank.items():
            n=len(group);chosen=[r for r in group if r["selected"]]
            summary.append(dict(rank=rank,mode=mode,n_total_species=n,
                n_with_OriH_prediction=sum(r["n_OriH_features"]>0 for r in group),
                n_with_any_gene_free_gap=sum(r["n_gene_free_gaps"]>0 for r in group),
                n_selected_full_containment=len(chosen),
                pct_all_species_selected=round(100*len(chosen)/n,3) if n else "",
                n_selected_ge300=sum(r["selected_gap_bp"]>=300 for r in chosen),
                median_selected_length_bp=(
                    sorted(r["selected_gap_bp"] for r in chosen)[len(chosen)//2]
                    if chosen else ""),
                n_no_full_containment=sum(r["status"]=="no_full_OriH_containment" for r in group)))
            print(rank,mode,len(chosen),"/",n,flush=True)
    write(out/"all_ranks_species_strict_OriH.tsv",rows)
    write(out/"all_ranks_summary_strict_OriH.tsv",summary)
    print("Results:",out.resolve())

if __name__=="__main__":main()
