#!/usr/bin/env python3
"""Summarize why strict OriH-contained gene-free region selection fails.

Input produced by analyze_strict_orih_gene_deserts.py. One record per
rank/species/mode; all-species denominator, human excluded upstream.

Full containment is checked BEFORE selecting the longest qualifying gap.
Therefore predicted-OriH failures are NOT caused by choosing the longest
unfiltered gap ahead of time.

We can categorize definite size impossibility from the QC TSV, but when
OriH could fit in a gap by length alone, positional conflict versus an
annotated gene requires inspecting original coordinates/GFF.
"""
import argparse
from pathlib import Path
import pandas as pd

RANKS=("kingdom","phylum","subphylum","superclass","class","superorder",
       "order","suborder","infraorder","parvorder")
MODES=("no_human_homologous_genes","no_genes_at_all")

def classify(d):
    d=d.copy()
    for field in ("selected","mitos2_gff_available","n_OriH_features",
                  "n_gene_free_gaps","OriH_annotated_bp","largest_unfiltered_gap_bp",
                  "n_fully_containing_gaps"):
        d[field]=pd.to_numeric(d[field],errors="coerce").fillna(0)
    def group(row):
        if row["selected"]==1:
            return "selected_full_containment"
        if row["mitos2_gff_available"]==0:
            return "missing_mitos2_gff"
        if row["n_OriH_features"]==0:
            return "no_OriH_prediction"
        if row["n_gene_free_gaps"]==0:
            return "OriH_predicted_no_gene_free_gap"
        if row["OriH_annotated_bp"]>row["largest_unfiltered_gap_bp"]:
            return "OriH_predicted_too_long_for_any_gap"
        return "OriH_predicted_no_single_gap_fully_contains"
    d["failure_group"]=d.apply(group,axis=1)
    d["OriH_predicted"]=d["n_OriH_features"]>0
    d["predicted_but_not_selected"]=d["OriH_predicted"] & (d["selected"]!=1)
    d["multiple_OriH_features"]=d["n_OriH_features"]>1
    return d

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input",default="strict_orih_gene_desert_qc/all_ranks_species_strict_OriH.tsv")
    ap.add_argument("--out-dir",default="strict_orih_gene_desert_qc/failure_diagnostics")
    args=ap.parse_args()
    out=Path(args.out_dir);out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(args.input,sep="\t",low_memory=False)
    required={"rank","species","mode","selected","mitos2_gff_available",
              "n_OriH_features","n_gene_free_gaps","OriH_annotated_bp",
              "largest_unfiltered_gap_bp","n_fully_containing_gaps"}
    if required-set(d.columns):
        raise ValueError("Missing expected QC columns: "+str(required-set(d.columns)))
    if d.duplicated(["rank","species","mode"]).any():
        raise ValueError("Repeated rank/species/mode records")
    d=classify(d)
    d.to_csv(out/"all_species_failure_categories.tsv",sep="\t",index=False)
    failures=d[d["predicted_but_not_selected"]].copy()
    failures.to_csv(out/"OriH_predicted_but_not_selected.tsv",sep="\t",index=False)

    summary=[]
    detail=[]
    for mode in MODES:
        for rank in RANKS:
            x=d[(d["rank"]==rank)&(d["mode"]==mode)]
            if x.empty:continue
            counts=x.failure_group.value_counts()
            n=len(x)
            n_pred=int(x.OriH_predicted.sum())
            n_chosen=int((x.selected==1).sum())
            n_fail=int(x.predicted_but_not_selected.sum())
            if n_pred!=n_chosen+n_fail:
                raise ValueError(f"Inconsistent OriH counts for {rank}/{mode}")
            p=x[x.predicted_but_not_selected]
            summary.append(dict(
                mode=mode,rank=rank,n_all_species=n,
                n_mitos2_gff_missing=int(counts.get("missing_mitos2_gff",0)),
                n_no_OriH_prediction=int(counts.get("no_OriH_prediction",0)),
                n_OriH_predicted=n_pred,
                n_OriH_selected_full_containment=n_chosen,
                n_OriH_predicted_but_missed=n_fail,
                pct_predicted_but_missed_among_all_species=100*n_fail/n,
                pct_predicted_but_missed_among_predicted=100*n_fail/n_pred if n_pred else float("nan"),
                n_predicted_no_gene_free_gap=int(counts.get("OriH_predicted_no_gene_free_gap",0)),
                n_predicted_OriH_too_long_for_any_gap=int(counts.get("OriH_predicted_too_long_for_any_gap",0)),
                n_predicted_no_single_gap_fully_contains=int(counts.get("OriH_predicted_no_single_gap_fully_contains",0)),
                n_failed_with_multiple_OriH_features=int(p.multiple_OriH_features.sum()),
            ))
            for category,count in counts.items():
                detail.append(dict(mode=mode,rank=rank,failure_group=category,
                                   n_species=int(count),pct_all_species=100*count/n))
    pd.DataFrame(summary).to_csv(out/"OriH_failure_summary_by_rank.tsv",sep="\t",index=False)
    pd.DataFrame(detail).to_csv(out/"failure_categories_by_rank.tsv",sep="\t",index=False)
    with pd.option_context("display.max_rows",None,"display.width",240):
        print(pd.DataFrame(summary)[["mode","rank","n_all_species",
            "n_OriH_predicted","n_OriH_selected_full_containment",
            "n_OriH_predicted_but_missed",
            "n_predicted_no_gene_free_gap",
            "n_predicted_OriH_too_long_for_any_gap",
            "n_predicted_no_single_gap_fully_contains",
            "n_failed_with_multiple_OriH_features"]].to_string(index=False))
    print("\nWrote:",out)
    print("NOTE: 'no_single_gap_fully_contains' is a positional failure, "
          "not proof of gene overlap; inspect GFF/coordinates to distinguish "
          "gene overlap, spanning two gaps, and possible annotation errors.")

if __name__=="__main__":
    main()
