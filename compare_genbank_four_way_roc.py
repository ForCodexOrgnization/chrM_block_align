#!/usr/bin/env python3
"""GenBank+MITOS2 four-way genome-wide ROC across all ten ranks.

Reuses the existing four-way ROC implementation and exactly the same
ClinVar/MITOMAP labeling, full-case variant cohort, position-clustered
bootstrap, and plotting logic. Only the result directories and filenames
differ. Existing results are read-only.
"""
import argparse
from pathlib import Path
import sys
from compare_all_ranks_padding_roc import main as existing_main

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--old-csv",default="data/2026_10_10_mito_cons_anns_main.csv")
    p.add_argument("--results-dir",default="block_results_genbank_mitos_v1")
    p.add_argument("--out-dir",default="block_results_genbank_mitos_v1/roc_comparison")
    p.add_argument("--metric",choices=["CI_MSA","CI_callable"],default="CI_MSA")
    p.add_argument("--bootstrap-reps",type=int,default=1000)
    p.add_argument("--seed",type=int,default=20261030)
    a=p.parse_args()
    root=Path(a.results_dir)
    # The existing runner expects --padded-dir/{pad50,pad100} and
    # --gene-only-dir for the three new CI sources.
    # Namespace is kept separate to avoid overwriting the old figures.
    sys.argv=[
        "compare_all_ranks_padding_roc.py",
        "--old-csv",a.old_csv,
        "--padded-dir",str(root),
        "--gene-only-dir",str(root/"gene_only"),
        "--out-dir",a.out_dir,
        "--metric",a.metric,
        "--bootstrap-reps",str(a.bootstrap_reps),
        "--seed",str(a.seed),
    ]
    existing_main()
    out=Path(a.out_dir)
    for stem in ("all_ranks_four_way_roc.png",
                 "all_ranks_four_way_roc_summary.csv",
                 "all_ranks_four_way_roc_bootstrap.csv"):
        old=out/stem
        if old.is_file():
            old.rename(out/("genbank_"+stem))
    print("GenBank+MITOS2 four-way ROC outputs:",out)

if __name__=="__main__":main()
