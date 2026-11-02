#!/usr/bin/env python3
"""Scatter old vs new chrM CI for gene-only, pad50, pad100 across ten ranks.

Works for MITOS2 and GenBank. Missing rank/method CI files are marked unavailable.
Plots one scatter panel per rank, reporting position-level Pearson and Spearman.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RANKS=("kingdom","phylum","subphylum","superclass","class","superorder",
       "order","suborder","infraorder","parvorder")
METHODS=("gene_only","pad50","pad100")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--old-csv",default="data/2026_10_10_mito_cons_anns_main.csv")
    p.add_argument("--pipeline",choices=["mitos2","genbank"],required=True)
    p.add_argument("--mitos2-gene-only-root",default="block_results_v6")
    p.add_argument("--mitos2-padded-root",default="block_results_padded_v1")
    p.add_argument("--genbank-root",default="block_results_genbank_mitos_v1")
    p.add_argument("--metric",choices=["CI_MSA","CI_callable"],default="CI_MSA")
    p.add_argument("--out-dir",default=None)
    a=p.parse_args()
    output=Path(a.out_dir or f"old_vs_new_ci_{a.pipeline}")
    output.mkdir(parents=True,exist_ok=True)
    cols=["POS"]+[f"{r}_rot_ref" for r in RANKS]
    old=pd.read_csv(a.old_csv,usecols=cols,low_memory=False)
    old.POS=pd.to_numeric(old.POS,errors="raise").astype(int)
    if not old.POS.between(1,16569).all():raise ValueError("Old POS must be 1..16569")
    # One observation per human position; avoid giving sites with many alleles extra weight.
    for c in cols[1:]:old[c]=pd.to_numeric(old[c],errors="coerce")
    old=old.groupby("POS",sort=True)[cols[1:]].median()
    rows=[]
    for method in METHODS:
        fig,axes=plt.subplots(2,5,figsize=(20,8),sharex=True,sharey=True)
        for rank,ax in zip(RANKS,axes.flat):
            if a.pipeline=="genbank":
                path=Path(a.genbank_root)/method/rank/"conservation"/"conservation_ci.csv"
            elif method=="gene_only":
                path=Path(a.mitos2_gene_only_root)/rank/"conservation"/"conservation_ci.csv"
            else:
                path=Path(a.mitos2_padded_root)/method/rank/"conservation"/"conservation_ci.csv"
            ax.set_title(rank)
            ax.set_xlim(0,100);ax.set_ylim(0,100)
            ax.plot([0,100],[0,100],"--",color="#777777",lw=.7)
            if not path.is_file():
                print(f"[SKIP] {method}/{rank}: {path}")
                ax.text(.5,.5,"Not available",ha="center",va="center",transform=ax.transAxes)
                rows.append(dict(pipeline=a.pipeline,method=method,rank=rank,
                    n_positions=0,pearson_r=np.nan,spearman_rho=np.nan,status="missing"))
                continue
            new=pd.read_csv(path,usecols=["Position",a.metric])
            new.Position=pd.to_numeric(new.Position,errors="raise").astype(int)
            if len(new)!=16569 or new.Position.tolist()!=list(range(1,16570)):
                raise ValueError(f"Invalid human coordinate map: {path}")
            x=old[f"{rank}_rot_ref"]
            y=new.set_index("Position")[a.metric].pipe(pd.to_numeric,errors="coerce")
            xy=pd.concat([x.rename("old"),y.rename("new")],axis=1).dropna()
            xy=xy.loc[np.isfinite(xy.old)&np.isfinite(xy.new)]
            if not xy.empty and not xy[["old","new"]].stack().between(0,100).all():
                raise ValueError(f"CI must be 0..100: {path}")
            pearson=xy.old.corr(xy.new,method="pearson") if len(xy)>2 else np.nan
            spearman=xy.old.corr(xy.new,method="spearman") if len(xy)>2 else np.nan
            ax.scatter(xy.old,xy.new,s=4,alpha=.14,rasterized=True,color="#286BA1",linewidths=0)
            ax.text(.04,.96,f"Pearson r = {pearson:.3f}\nSpearman ρ = {spearman:.3f}\nn = {len(xy):,}",
                    transform=ax.transAxes,ha="left",va="top",fontsize=9,
                    bbox=dict(facecolor="white",edgecolor="none",alpha=.85))
            rows.append(dict(pipeline=a.pipeline,method=method,rank=rank,
                n_positions=len(xy),pearson_r=pearson,spearman_rho=spearman,status="available"))
        for ax in axes[-1]:ax.set_xlabel("Old chrM CI (0–100)")
        for row in axes:row[0].set_ylabel("New chrM CI (0–100)")
        fig.suptitle(f"Old vs {a.pipeline.upper()} {method} chrM CI — by taxonomic rank")
        fig.tight_layout()
        for ext in ("png","pdf"):
            fig.savefig(output/f"old_vs_{a.pipeline}_{method}_dotplots.{ext}",dpi=200)
        plt.close(fig)
    pd.DataFrame(rows).to_csv(output/"old_vs_new_ci_correlation_summary.csv",index=False)
    print("Saved figures and summary to",output)

if __name__=="__main__":
    main()
