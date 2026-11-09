#!/usr/bin/env python3
"""Coordinate-level audit of strict OriH candidate failures.

Run on McCleary, where MITOS2 GFF paths in mitos2_species_summary.tsv exist.
Produces one row per failed rank/species/mode plus original OriH GFF features
and per-source gene conflicts. Does not modify alignments or prior outputs.

Important: GFF start/end checks prove annotation geometry, not which
coordinate system/assembly produced the file. Multiple OH features are
required to fit the same gap, matching strict pipeline semantics.
"""
import argparse
import collections
import csv
from pathlib import Path

import chrM_block_align as base
import resolve_gene_order_annotations as go
import run_genbank_mitos_blocks as pipeline
from analyze_longest_all_gene_desert import fasta, norm, number, spans, union, overlap, gff_paths, write
from analyze_strict_orih_gene_deserts import (
    CANONICAL, OH_LABELS, TYPES, SKIP_LABELS, cleaned,
    read_genbank, canonical_mitos_trna, circular_gaps,
)

def fmt(start,end,L):
    return f"{start%L+1}..{end%L}" + (" (wrap)" if end>L else "")

def pieces_overlap(piecesA,piecesB,L):
    return sum(overlap(a,b,c,d,L) for a,b in piecesA for c,d in piecesB)

def annotated_features(gff,L):
    all_gene=[];canonical=[];origin=[];raw_origin=[];gff_gene=[]
    with open(gff,encoding="utf-8",errors="replace") as f:
        for line_number,line in enumerate(f,1):
            if not line.strip() or line.startswith("#"):continue
            p=line.rstrip("\n").split("\t")
            if len(p)<9:continue
            ftype=p[2].lower()
            attrs={k.casefold():v for item in p[8].split(";") if "=" in item
                   for k,v in [item.split("=",1)]}
            labels=[cleaned(attrs.get(k,"")) for k in ("name","gene_id","id","product")]
            a,b=number(p[3]),number(p[4])
            loci=spans(a,b,L)
            oh=(ftype in ("origin_of_replication","rep_origin")
                and any(label in OH_LABELS for label in labels))
            if oh:
                raw_origin.append(dict(
                    gff_line=line_number,seqid=p[0],source=p[1],
                    type=p[2],start=a,end=b,strand=p[6],
                    raw_length_bp=(b-a+1 if a is not None and b is not None and b>=a else
                                   L-a+1+b if a is not None and b is not None else ""),
                    coordinate_valid=int(bool(loci)),raw_attributes=p[8]))
                if loci:origin.extend(loci)
                continue
            if ftype not in TYPES or any(label in SKIP_LABELS for label in labels) or not loci:
                continue
            all_gene.extend(loci)
            gene=None
            for field in ("name","gene_id","id"):
                gene=base.canonical_mitos(attrs.get(field,""))
                if gene:break
            if gene is None:gene=canonical_mitos_trna(attrs)
            canonical_gene=gene in CANONICAL
            if canonical_gene:canonical.extend(loci)
            gff_gene.append((gene or attrs.get("name","") or attrs.get("id","") or "unknown",
                             loci,canonical_gene,p[2],p[3],p[4]))
    return all_gene,canonical,origin,raw_origin,gff_gene

def read_genbank_labeled(path,aliases):
    """Recompute exactly the same inclusion as strict workflow, with labels."""
    result=collections.defaultdict(list)
    with open(path,newline="",encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            k=norm(row.get("Organism") or row.get("organism") or row.get("species") or "")
            a,b=number(row.get("start")),number(row.get("end"))
            if not k or a is None or b is None:continue
            typ=str(row.get("feature") or row.get("Feature") or row.get("type") or "").strip().lower()
            if typ and typ not in TYPES:continue
            labels=[row.get(x,"") for x in ("Gene","G_Symbol","G_locus_tag","alt_name","alt_Product")]
            if any(cleaned(v) in SKIP_LABELS for v in labels):continue
            is_canonical=bool(go.candidates(row,aliases)&CANONICAL)
            gene=next((str(v) for v in labels if v),"unknown")
            result[k].append((a,b,is_canonical,gene))
    return result

def read_trnas_labeled(path):
    result=collections.defaultdict(list)
    with open(path,newline="",encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            k=norm(row.get("Organism",""))
            a,b=number(row.get("Begin")),number(row.get("End"))
            if k and a is not None and b is not None:
                result[k].append((min(a,b),max(a,b),
                                  str(row.get("Type","")+"/"+row.get("AntiCodon",""))))
    return result

def audit(occupied,oh,L,genes):
    occupied_union=union(occupied)
    gaplist=circular_gaps(occupied,L)
    oh_union=union(oh)
    oh_bp=sum(b-a for a,b in oh_union)
    gap_best=max(gaplist,key=lambda t:t[1]-t[0]) if gaplist else None
    # Each gene entry has (label, pieces, origin, canonical)
    conflicts=[]
    for label,loci,source,canonical in genes:
        bp=pieces_overlap(oh_union,loci,L)
        if bp:
            conflicts.append(dict(gene=label,source=source,
                                  canonical=canonical,overlap_bp=bp))
    conflict_by_source=collections.Counter()
    for piece in oh_union:
        for label,loci,source,canonical in genes:
            if pieces_overlap([piece],loci,L)>0:
                conflict_by_source[source]+=1
    # Per OH base: "covered" by merged occupied union (not double counted).
    oh_gene_overlap=pieces_overlap(oh_union,occupied_union,L)
    best_intersection=0;closest=None
    for a,b in gaplist:
        covered=pieces_overlap([(a,b)],oh_union,L)
        if covered>best_intersection:
            best_intersection=covered;closest=(a,b)
    return dict(n_gene_blocks=len(occupied_union),
                n_gene_free_gaps=len(gaplist),
                OriH_union_bp=oh_bp,OriH_gene_overlap_bp=oh_gene_overlap,
                OriH_gene_overlap_fraction=round(oh_gene_overlap/oh_bp,6) if oh_bp else "",
                longest_gap_bp=(gap_best[1]-gap_best[0] if gap_best else 0),
                longest_gap_coordinates=fmt(*gap_best,L) if gap_best else "",
                best_OriH_gap_overlap_bp=best_intersection,
                best_OriH_gap_coverage_fraction=round(best_intersection/oh_bp,6) if oh_bp else "",
                best_OriH_gap_coordinates=fmt(*closest,L) if closest else "",
                n_annotated_gene_features_overlapping_OriH=len(conflicts),
                n_GenBank_conflict=sum(c["source"]=="GenBank" for c in conflicts),
                n_MITOS2_conflict=sum(c["source"]=="MITOS2" for c in conflicts),
                n_tRNAscan_conflict=sum(c["source"]=="tRNAscan" for c in conflicts),
                n_OriH_noncontained_bp=oh_bp-best_intersection)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--failures",default="strict_orih_gene_desert_qc/failure_diagnostics/OriH_predicted_but_not_selected.tsv")
    p.add_argument("--fasta-dir",default="data/unrot_fasta")
    p.add_argument("--genbank-annotations",default="data/kingdom_genes_v2.csv")
    p.add_argument("--gene-aliases",default="data/genefiltfinal5_2026_09_24.csv")
    p.add_argument("--trnascan-csv",default="data/consolidated_trnascan-2.csv")
    p.add_argument("--mitos2-summary",default="data/mitos2_species_summary.tsv")
    p.add_argument("--ranks",nargs="+",default=["kingdom"])
    p.add_argument("--out-dir",default="strict_orih_gene_desert_qc/coordinate_audit")
    a=p.parse_args()
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    failure=collections.defaultdict(set)
    with open(a.failures,newline="") as fh:
        for row in csv.DictReader(fh,delimiter="\t"):
            if row["rank"] in a.ranks and row["predicted_but_not_selected"].lower()=="true":
                failure[row["rank"],norm(row["species"])].add(row["mode"])
    if not failure:raise ValueError("No requested failed species in input")
    aliases=pipeline.load_aliases(a.gene_aliases)
    gb=read_genbank_labeled(a.genbank_annotations,aliases)
    tr=read_trnas_labeled(a.trnascan_csv)
    gff_map=gff_paths(a.mitos2_summary)
    main_rows=[];oh_rows=[];conflict_rows=[]
    cache={}
    for rank in a.ranks:
        for species,seq in fasta(Path(a.fasta_dir)/(rank+"_chrm.faa")):
            key=norm(species)
            modes=failure.get((rank,key))
            if not modes:continue
            L=len(seq)
            path=gff_map.get(key)
            if not path or not Path(path).is_file():
                print(f"[SKIP missing GFF] {rank} {species}: {path}",flush=True)
                continue
            if (key,L) not in cache:
                cache[key,L]=annotated_features(path,L)
            m_all,m_canon,oh,raw_oh,gff_genes=cache[key,L]
            for g in raw_oh:
                oh_rows.append(dict(rank=rank,species=species,genome_bp=L,
                                    gff_path=path,**g))
            def features_for(mode):
                gs=[];occupied=[]
                for g,loci,is_canon,ftype,start,end in gff_genes:
                    if mode=="no_genes_at_all" or is_canon:
                        gs.append((g,loci,"MITOS2",is_canon));occupied.extend(loci)
                for first,last,is_canon,g in gb.get(key,[]):
                    if mode=="no_genes_at_all" or is_canon:
                        loci=spans(first,last,L)
                        if loci:
                            occupied.extend(loci);gs.append((g,loci,"GenBank",is_canon))
                for first,last,name in tr.get(key,[]):
                    loci=spans(first,last,L)
                    if loci:
                        occupied.extend(loci);gs.append((name,loci,"tRNAscan",True))
                return occupied,gs
            for mode in sorted(modes):
                occupied,genes=features_for(mode)
                stats=audit(occupied,oh,L,genes)
                valid_gff=all(r["coordinate_valid"] for r in raw_oh)
                main_rows.append(dict(rank=rank,species=species,mode=mode,genome_bp=L,
                    gff_path=path,n_oriH_gff_features=len(raw_oh),
                    all_OriH_GFF_coordinates_valid=int(valid_gff),
                    n_OriH_coordinate_pieces_parsed=len(oh),
                    OriH_fraction_genome=round(stats["OriH_union_bp"]/L,6),
                    **stats))
                oh_union=union(oh)
                for name,loci,source,is_canon in genes:
                    bp=pieces_overlap(oh_union,loci,L)
                    if bp:
                        conflict_rows.append(dict(rank=rank,species=species,mode=mode,
                            source=source,gene=name,overlap_bp=bp,
                            overlaps_OriH_fraction=round(bp/max(stats["OriH_union_bp"],1),6)))
    if not main_rows:raise ValueError("No FASTA species matched failed cases")
    write(out/"failed_species_OriH_coordinates_and_gene_conflicts.tsv",main_rows)
    if oh_rows:write(out/"raw_MITOS2_OriH_GFF_features.tsv",oh_rows)
    if conflict_rows:write(out/"genes_overlapping_predicted_OriH.tsv",conflict_rows)
    print("Wrote coordinate audit to",out)
    print("Inspect: raw_MITOS2_OriH_GFF_features.tsv (actual start/end), "
          "failed_species_OriH_coordinates_and_gene_conflicts.tsv, "
          "genes_overlapping_predicted_OriH.tsv")
    print("Caution: broad GFF intervals are not alone proof of an annotation bug;"
          " verify original GFF coordinate system and FASTA sequence identity.")

if __name__=="__main__":
    main()
