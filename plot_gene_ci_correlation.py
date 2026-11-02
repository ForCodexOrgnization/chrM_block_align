#!/usr/bin/env python3
"""Compare old chrM CI and newly aligned chrM CI to rank-specific gene CI.

The source annotation CSV calls the rank-specific gene score <rank>_ref_gene.
Some duplicate POS/REF/ALT rows have conflicting scores; report such rows
and aggregate by per-variant median before correlation to avoid extra votes.
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

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--annotation-csv",default="data/2026_10_10_mito_cons_anns_main.csv")
    p.add_argument("--new-root",default="block_results_genbank_mitos_v1")
    p.add_argument("--methods",nargs="+",default=["gene_only","pad50","pad100"])
    p.add_argument("--new-metric",choices=["CI_MSA","CI_callable"],default="CI_MSA")
    p.add_argument("--out-dir",default="block_results_genbank_mitos_v1/gene_ci_correlations")
    p.add_argument("--old-only",action="store_true")
    args=p.parse_args();out=Path(args.out_dir);out.mkdir(parents=True,exist_ok=True)
    cols=["POS","REF","ALT"]+[f"{r}_{suffix}" for r in RANKS for suffix in ("ref_gene","rot_ref")]
    df=pd.read_csv(args.annotation_csv,usecols=cols,low_memory=False)
    df.POS=pd.to_numeric(df.POS,errors="raise").astype(int)
    grouped=df.groupby(["POS","REF","ALT"],sort=True,dropna=False)
    conflicted=grouped[cols[3:]].nunique(dropna=True).gt(1)
    flags=conflicted.stack().loc[lambda x:x].reset_index()
    flags.columns=["POS","REF","ALT","score_column","conflict"]
    flags.drop(columns="conflict").to_csv(out/"duplicate_variant_score_conflicts.csv",index=False)
    variants=grouped[cols[3:]].median().reset_index()
    methods=["old"] if args.old_only else ["old"]+args.methods
    summaries=[]
    for method in methods:
        fig,axs=plt.subplots(2,5,figsize=(22,9),sharex=True,sharey=True)
        for rank,ax in zip(RANKS,axs.flat):
            x=variants[f"{rank}_ref_gene"]
            if method=="old":y=variants[f"{rank}_rot_ref"]
            else:
                path=Path(args.new_root)/method/rank/"conservation"/"conservation_ci.csv"
                if not path.exists():raise FileNotFoundError(path)
                table=pd.read_csv(path,usecols=["Position",args.new_metric])
                if len(table)!=16569 or table.Position.nunique()!=16569:
                    raise ValueError(f"Expected 16,569 unique positions: {path}")
                table=table.rename(columns={"Position":"POS",args.new_metric:"score"})
                y=variants[["POS"]].merge(table,on="POS",how="left",validate="many_to_one").score
            good=x.notna()&y.notna()&np.isfinite(x)&np.isfinite(y)
            xv=x[good].to_numpy(float);yv=y[good].to_numpy(float)
            pearson=float(np.corrcoef(xv,yv)[0,1]) if len(xv)>2 and np.std(xv)>0 and np.std(yv)>0 else float("nan")
            spearman=float(pd.Series(xv).corr(pd.Series(yv),method="spearman")) if len(xv)>2 else float("nan")
            ax.hexbin(xv,yv,gridsize=45,mincnt=1,bins="log",cmap="Blues",extent=(0,100,0,100))
            ax.plot([0,100],[0,100],"--",color="gray",lw=.8)
            ax.text(.04,.96,f"Pearson r = {pearson:.3f}\nSpearman = {spearman:.3f}\nn = {len(xv):,}",
                    transform=ax.transAxes,va="top",fontsize=9,
                    bbox=dict(facecolor="white",edgecolor="none",alpha=.85))
            ax.set_title(rank);ax.set_xlim(0,100);ax.set_ylim(0,100)
            summaries.append(dict(method=method,rank=rank,n_variants=len(xv),pearson_r=pearson,spearman_rho=spearman))
        for ax in axs[-1]:ax.set_xlabel("Gene CI (0-100)")
        for row in axs:row[0].set_ylabel("chrM CI (0-100)")
        fig.suptitle(f"{method}: chrM CI vs gene CI (unique POS/REF/ALT variants)")
        fig.tight_layout()
        fig.savefig(out/f"{method}_vs_gene_ci.png",dpi=180)
        fig.savefig(out/f"{method}_vs_gene_ci.pdf")
        plt.close(fig)
    result=pd.DataFrame(summaries)
    result.to_csv(out/"gene_ci_correlation_summary.csv",index=False)
    print(result.to_string(index=False))
if __name__=="__main__":main()
