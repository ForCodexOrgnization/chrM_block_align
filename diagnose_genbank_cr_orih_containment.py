#!/usr/bin/env python3
"""Diagnose apparent low 100% GenBank CR OriH containment from species TSV.

No coordinate transformation is applied. This report distinguishes missing CR,
zero overlap, partial overlap, and full overlap, using the exact single-region
maximal OriH intersection already produced by the comparison script.
"""
import argparse
import collections
import csv
from pathlib import Path

from compare_pairwise_cr_to_longest_intervals import RANKS


def label(r):
    if r["ori_status"]!="predicted":
        return "OriH_unavailable"
    length=int(r["ori_bp"])
    n=int(r["n_genbank_control_regions"])
    if n==0:
        return "no_full_genbank_CR_annotation"
    if r["coordinate_qc"]=="multiple_genbank_accessions":
        return "ambiguous_genbank_accessions"
    bp=int(r["GenBank_CR_max_single_region_overlap_bp"] or 0)
    if bp==0:
        return "zero_overlap"
    if bp==length:
        return "100pct_contained"
    pct=100*bp/length
    if pct>=99:
        return "99_to_under100pct"
    if pct>=95:
        return "95_to_under99pct"
    if pct>=80:
        return "80_to_under95pct"
    return "gt0_to_under80pct"


def make_qc(path,out_dir):
    out=Path(out_dir)
    out.mkdir(parents=True,exist_ok=True)
    with open(path,newline="",encoding="utf-8-sig") as fh:
        rows=list(csv.DictReader(fh,delimiter="\t"))
    for r in rows:
        r["genbank_failure_class"]=label(r)
        if r["ori_status"]=="predicted":
            ori=int(r["ori_bp"])
            crbp=int(r["GenBank_CR_max_single_region_overlap_bp"] or 0)
            r["max_single_CR_coverage_pct"]=round(crbp*100/ori,4) if ori else ""
        else:
            r["max_single_CR_coverage_pct"]=""
    with (out/"genbank_CR_containment_species_diagnostics.tsv").open("w",newline="") as fh:
        writer=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    counts=[]
    for rank in RANKS:
        subset=[r for r in rows if r["rank"]==rank]
        pred=[r for r in subset if r["ori_status"]=="predicted"]
        categories=collections.Counter(r["genbank_failure_class"] for r in pred)
        all_n=len(subset)
        pred_n=len(pred)
        cr_present=[r for r in pred if int(r["n_genbank_control_regions"])>0
                    and r["coordinate_qc"]!="multiple_genbank_accessions"]
        if pred_n:
            counts.append(dict(
                rank=rank,n_all_species=all_n,n_orih_predicted=pred_n,
                n_orih_with_genbank_CR_annotation=len(cr_present),
                n_genbank_full=sum(r["GenBank_CR_full"]=="1" for r in pred),
                n_longest_full=sum(r["Longest_full"]=="1" for r in pred),
                genbank_full_pct_of_predicted=round(100*sum(r["GenBank_CR_full"]=="1" for r in pred)/pred_n,3),
                genbank_full_pct_given_CR=round(100*sum(r["GenBank_CR_full"]=="1" for r in cr_present)/len(cr_present),3) if cr_present else "",
                n_genbank_zero_overlap=categories["zero_overlap"],
                n_genbank_no_annotation=categories["no_full_genbank_CR_annotation"],
                n_genbank_ambiguous_accessions=categories["ambiguous_genbank_accessions"],
                n_genbank_0_80=categories["gt0_to_under80pct"],
                n_genbank_80_95=categories["80_to_under95pct"],
                n_genbank_95_99=categories["95_to_under99pct"],
                n_genbank_99_100=categories["99_to_under100pct"],
                n_genbank_full_but_longest_not=sum(
                    r["GenBank_CR_full"]=="1" and r["Longest_full"]!="1" for r in pred),
                n_longest_full_but_genbank_not=sum(
                    r["Longest_full"]=="1" and r["GenBank_CR_full"]!="1" for r in pred)))
    if counts:
        with (out/"genbank_CR_containment_rank_diagnostics.tsv").open("w",newline="") as fh:
            writer=csv.DictWriter(fh,fieldnames=list(counts[0]),delimiter="\t")
            writer.writeheader()
            writer.writerows(counts)
    print("Wrote failure-class diagnostics to",out)
    return rows,counts


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-tsv",
        default="longest_all_gene_desert_qc/genbank_CR_comparison/species_100pct_containment.tsv")
    parser.add_argument("--out-dir",
        default="longest_all_gene_desert_qc/genbank_CR_comparison")
    args=parser.parse_args()
    make_qc(args.species_tsv,args.out_dir)


if __name__=="__main__":
    main()
