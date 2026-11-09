#!/usr/bin/env python3
"""Read-only GenBank coordinate audit for failed strict OriH gene-desert cases.

Replicates the gene-footprint union of analyze_strict_orih_gene_deserts.py.
Focuses on zero-gap cases and flags GenBank start>end, very long interpreted
spans, and out-of-range coordinates. Also reports a COUNTERFACTUAL union in
which start>end GenBank intervals are treated as [end,start] (strand is NOT
sufficient evidence that this alternative is biologically correct).

Output:
  genbank_coordinate_rows.tsv        per source GenBank annotation
  gene_coverage_by_species.tsv       observed vs hypothetical footprint coverage
  suspicious_genbank_rows.tsv        high-priority records to inspect
All coordinates are from the original unrotated FASTA/GenBank CSV.
"""
import argparse, collections, csv, re
from pathlib import Path

import chrM_block_align as base
import resolve_gene_order_annotations as go
import run_genbank_mitos_blocks as pipeline
from analyze_longest_all_gene_desert import (
    RANKS, TYPES, norm, number, spans, union, fasta, trnascan_genes,
    gff_paths, write,
)
from analyze_strict_orih_gene_deserts import (
    CANONICAL, SKIP_LABELS, cleaned, mitos2_features, circular_gaps,
)

MODES=("no_human_homologous_genes","no_genes_at_all")

def union_bp(pieces):
    return sum(b-a for a,b in union(pieces))

def read_raw_genbank(path,aliases):
    by=collections.defaultdict(list)
    with open(path,newline="",encoding="utf-8-sig") as f:
        for n,row in enumerate(csv.DictReader(f),2):
            key=norm(row.get("Organism") or row.get("organism") or row.get("species") or "")
            if not key:continue
            typ=str(row.get("feature") or row.get("Feature") or row.get("type") or "").strip().lower()
            labels=[row.get(x,"") for x in ("Gene","G_Symbol","G_locus_tag","alt_name","alt_Product")]
            # Match the strict pipeline's GenBank source filtering exactly.
            if typ and typ not in TYPES:continue
            if any(cleaned(v) in SKIP_LABELS for v in labels):continue
            gene=next((str(x) for x in labels if x),"unknown")
            by[key].append(dict(line=n,acc=row.get("acc",""),gene=gene,
                feature_type=typ,start=number(row.get("start")),
                end=number(row.get("end")),strand=str(row.get("strand","")).strip(),
                canonical=bool(go.candidates(row,aliases)&CANONICAL)))
    return by

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--species-tsv",
                   default="strict_orih_gene_desert_qc/all_ranks_species_strict_OriH.tsv")
    p.add_argument("--genbank-annotations",default="data/kingdom_genes_v2.csv")
    p.add_argument("--gene-aliases",default="data/genefiltfinal5_2026_09_24.csv")
    p.add_argument("--fasta-dir",default="data/unrot_fasta")
    p.add_argument("--trnascan-csv",default="data/consolidated_trnascan-2.csv")
    p.add_argument("--mitos2-summary",default="data/mitos2_species_summary.tsv")
    p.add_argument("--ranks",nargs="+",default=["kingdom"],choices=RANKS)
    p.add_argument("--only-zero-gap",action="store_true",
                   help="Only species/modes with n_gene_free_gaps == 0")
    p.add_argument("--min-gene-fraction",type=float,default=0.5,
                   help="Flag single GenBank gene intervals occupying this fraction of genome")
    p.add_argument("--out-dir",
                   default="strict_orih_gene_desert_qc/genbank_coordinate_audit")
    a=p.parse_args()
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    targets=collections.defaultdict(set)
    with open(a.species_tsv,newline="",encoding="utf-8-sig") as f:
        for row in csv.DictReader(f,delimiter="\t"):
            if row["rank"] not in a.ranks:continue
            mode=row["mode"]
            if mode not in MODES:continue
            if a.only_zero_gap and number(row.get("n_gene_free_gaps"))!=0:continue
            # Investigate failed OriH-containing selections only.
            if number(row.get("n_OriH_features")) is None or number(row["n_OriH_features"])<1:continue
            if number(row.get("selected"))!=0:continue
            targets[row["rank"],norm(row["species"])].add(mode)
    if not targets:raise ValueError("No matching failed species in QC TSV")
    aliases=pipeline.load_aliases(a.gene_aliases)
    gb=read_raw_genbank(a.genbank_annotations,aliases)
    tr=trnascan_genes(a.trnascan_csv)
    mitos_files=gff_paths(a.mitos2_summary)
    rows=[];coverage=[]
    cache={}
    seen=set()
    for rank in a.ranks:
        for name,seq in fasta(Path(a.fasta_dir)/(rank+"_chrm.faa")):
            key=norm(name)
            modes=targets.get((rank,key))
            if not modes:continue
            seen.add((rank,key));L=len(seq)
            gff=mitos_files.get(key)
            if (key,L) not in cache:
                cache[key,L]=mitos2_features(gff,L) if gff and Path(gff).is_file() else ([],[],[])
            m_all,m_canon,origins=cache[key,L]
            per_gene=[]
            for entry in gb.get(key,[]):
                start,end=entry["start"],entry["end"]
                valid=(start is not None and end is not None
                       and 1<=start<=L and 1<=end<=L)
                span=spans(start,end,L) if valid else []
                # Diagnostic-only alternative; do not treat "-" as evidence.
                alt=spans(min(start,end),max(start,end),L) if valid else []
                current_bp=sum(b-a for a,b in span)
                alt_bp=sum(b-a for a,b in alt)
                reasons=[]
                if not valid:reasons.append("invalid_or_out_of_range")
                if valid and start>end:reasons.append("start_gt_end")
                if valid and current_bp/L>=a.min_gene_fraction:
                    reasons.append("large_circular_span")
                if valid and start>end and entry["strand"]=="-":
                    reasons.append("reverse_strand_start_gt_end_review")
                if valid and start>end and current_bp>alt_bp*3:
                    reasons.append("circular_vs_sorted_size_discrepancy")
                data=dict(rank=rank,species=name,genome_bp=L,
                          csv_line=entry["line"],accession=entry["acc"],
                          gene=entry["gene"],feature_type=entry["feature_type"],
                          strand=entry["strand"],start=start,end=end,
                          is_human_homolog=int(entry["canonical"]),
                          valid_coordinates=int(valid),
                          start_gt_end=int(valid and start>end),
                          circular_interpretation_bp=current_bp if valid else "",
                          sorted_endpoint_interpretation_bp=alt_bp if valid else "",
                          circular_fraction_genome=round(current_bp/L,5) if valid else "",
                          review_flags=";".join(reasons),
                          suspicious=int(bool(reasons)))
                rows.append(data)
                per_gene.append((entry["canonical"],span,alt))
            tr_parts=[]
            for f in tr.get(key,[]):
                tr_parts.extend(spans(f.start+1,f.end,L))
            for mode in modes:
                canon=(mode=="no_human_homologous_genes")
                gb_original=[p for is_c,part,alt in per_gene if not canon or is_c for p in part]
                gb_alternative=[p for is_c,part,alt in per_gene if not canon or is_c for p in alt]
                non_gb=(m_canon if canon else m_all)+tr_parts
                current=non_gb+gb_original
                alternative=non_gb+gb_alternative
                current_gaps=circular_gaps(current,L)
                alternative_gaps=circular_gaps(alternative,L)
                coverage.append(dict(rank=rank,species=name,mode=mode,genome_bp=L,
                    n_GenBank_rows=len(gb.get(key,[])),
                    n_GenBank_start_gt_end=sum(
                        bool(e["start"] is not None and e["end"] is not None
                        and e["start"]>e["end"]) for e in gb.get(key,[])),
                    GenBank_union_bp=union_bp(gb_original),
                    GenBank_union_fraction=round(union_bp(gb_original)/L,5),
                    other_sources_union_bp=union_bp(non_gb),
                    original_all_sources_union_bp=union_bp(current),
                    original_all_sources_union_fraction=round(union_bp(current)/L,5),
                    original_gene_free_gaps=len(current_gaps),
                    original_longest_gap_bp=max((b-a for a,b in current_gaps),default=0),
                    counterfactual_all_sources_union_bp=union_bp(alternative),
                    counterfactual_gene_free_gaps=len(alternative_gaps),
                    counterfactual_longest_gap_bp=max((b-a for a,b in alternative_gaps),default=0),
                    counterfactual_recovers_gap=int(not current_gaps and bool(alternative_gaps)),
                    note="counterfactual is exploratory: does not establish correct strand semantics"))
    missing=set(targets)-seen
    if missing:raise ValueError(f"QC species missing from FASTA (first five): {list(missing)[:5]}")
    if not rows:raise ValueError("No GenBank annotation records for target species")
    write(out/"genbank_coordinate_rows.tsv",rows)
    write(out/"gene_coverage_by_species.tsv",coverage)
    write(out/"suspicious_genbank_rows.tsv",[r for r in rows if r["suspicious"]]
          or [dict.fromkeys(rows[0],"")])
    print("Wrote",out)
    print("Species/mode cases:",len(coverage))
    print("GenBank rows flagged:",sum(r["suspicious"] for r in rows))
    print("Counterfactual zero-gap cases recovering a gap:",
          sum(r["counterfactual_recovers_gap"] for r in coverage))
    print("IMPORTANT: start>end need not be an error; inspect original GenBank")
    print("feature location and match the exact FASTA accession before changes.")

if __name__=="__main__":
    main()
