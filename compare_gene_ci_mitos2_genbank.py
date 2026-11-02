#!/usr/bin/env python3
"""Compare gene CI to old, MITOS2, and GenBank chrM CI at ten ranks.

Methods:
  old: <rank>_rot_ref from the annotation CSV
  mitos2_gene_only: block_results_v6/<rank>/conservation/conservation_ci.csv
  mitos2_pad50/100: block_results_padded_v1/pad{50,100}/<rank>/...
  genbank_gene_only/pad50/pad100: block_results_genbank_mitos_v1/<method>/<rank>/...

Gene CI uses <rank>_ref_gene (actual column names in the annotation CSV).
Each point is a unique (POS, REF, ALT) variant. Duplicate contradictory
score records are flagged and resolved by median. Comparisons use the same
finite-value cohort PER RANK across all methods, so r differences are paired.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RANKS=("kingdom","phylum","subphylum","superclass","class",
       "superorder","order","suborder","infraorder","parvorder")
NEW_METHODS=("mitos2_gene_only","mitos2_pad50","mitos2_pad100",
             "genbank_gene_only","genbank_pad50","genbank_pad100")
METHODS=("old",)+NEW_METHODS

def new_ci_path(a,method,rank):
    parts=method.split("_",1)
    if parts[0]=="mitos2":
        root=Path(a.mitos2_gene_only_root) if parts[1]=="gene_only" else Path(a.mitos2_padded_root)/parts[1]
    else:
        root=Path(a.genbank_root)/parts[1]
    return root/rank/"conservation"/"conservation_ci.csv"

def load_ci(path,metric):
    if not path.exists():raise FileNotFoundError(f"Missing CI table: {path}")
    d=pd.read_csv(path,usecols=["Position",metric])
    d["Position"]=pd.to_numeric(d["Position"],errors="raise").astype(int)
    if len(d)!=16569 or d.Position.tolist()!=list(range(1,16570)):
        raise ValueError(f"Expected ordered 1..16569 human coordinates: {path}")
    d[metric]=pd.to_numeric(d[metric],errors="coerce")
    if not d[metric].dropna().between(0,100).all():
        raise ValueError(f"CI scores outside 0-100: {path}")
    return d.set_index("Position")[metric]

def read_variants(path):
    cols=["POS","REF","ALT"]+[f"{r}_{suffix}" for r in RANKS for suffix in ("ref_gene","rot_ref")]
    d=pd.read_csv(path,usecols=cols,low_memory=False)
    d.POS=pd.to_numeric(d.POS,errors="raise").astype(int)
    if not d.POS.between(1,16569).all():raise ValueError("Annotation POS outside human chrM")
    for col in cols[3:]:
        d[col]=pd.to_numeric(d[col],errors="coerce")
        if not d[col].dropna().between(0,100).all():
            raise ValueError(f"{col} has scores outside 0-100")
    grouped=d.groupby(["POS","REF","ALT"],sort=True,dropna=False)
    conflicted=grouped[cols[3:]].nunique(dropna=True).gt(1)
    flags=conflicted.stack().loc[lambda x:x].reset_index()
    flags.columns=["POS","REF","ALT","score_column","conflict"]
    return grouped[cols[3:]].median().reset_index(),flags.drop(columns="conflict")

def correlations(x,y):
    n=len(x)
    if n<3 or np.std(x)==0 or np.std(y)==0:return float("nan"),float("nan")
    return float(np.corrcoef(x,y)[0,1]),float(pd.Series(x).corr(pd.Series(y),method="spearman"))

def figure(d,method,out,metric):
    fig,axes=plt.subplots(2,5,figsize=(21,9),sharex=True,sharey=True)
    for rank,ax in zip(RANKS,axes.flat):
        x=d[f"{rank}_gene"].to_numpy(float)
        y=d[f"{rank}_{method}"].to_numpy(float)
        good=np.isfinite(x)&np.isfinite(y)
        x,y=x[good],y[good]
        r,rho=correlations(x,y)
        ax.hexbin(x,y,gridsize=45,mincnt=1,bins="log",cmap="Blues",extent=(0,100,0,100))
        ax.plot([0,100],[0,100],"--",color="#777777",linewidth=.8)
        ax.text(.04,.96,f"Pearson r = {r:.3f}\nSpearman \u03c1 = {rho:.3f}\nn = {len(x):,}",
                transform=ax.transAxes,va="top",fontsize=9,
                bbox={"facecolor":"white","alpha":.84,"edgecolor":"none"})
        ax.set_title(rank);ax.set_xlim(0,100);ax.set_ylim(0,100)
    for ax in axes[-1]:ax.set_xlabel("Gene CI (0-100)")
    for row in axes:row[0].set_ylabel("chrM CI (0-100)")
    fig.suptitle(f"{method.replace('_',' ')} chrM CI vs gene CI | {metric if method!='old' else 'original'}",
                 fontsize=15)
    fig.tight_layout()
    for ext in ("png","pdf"):
        fig.savefig(out/f"{method}_vs_gene_ci.{ext}",dpi=200)
    plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--annotation-csv",default="data/2026_10_10_mito_cons_anns_main.csv")
    p.add_argument("--mitos2-gene-only-root",default="block_results_v6")
    p.add_argument("--mitos2-padded-root",default="block_results_padded_v1")
    p.add_argument("--genbank-root",default="block_results_genbank_mitos_v1")
    p.add_argument("--metric",choices=["CI_MSA","CI_callable"],default="CI_MSA")
    p.add_argument("--out-dir",default="gene_ci_comparison_all_pipelines")
    a=p.parse_args()
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    d,conflicts=read_variants(a.annotation_csv)
    conflicts.to_csv(out/"duplicate_variant_score_conflicts.csv",index=False)
    scores={}
    for rank in RANKS:
        scores[f"{rank}_gene"]=d[f"{rank}_ref_gene"].to_numpy(float)
        scores[f"{rank}_old"]=d[f"{rank}_rot_ref"].to_numpy(float)
        for method in NEW_METHODS:
            s=load_ci(new_ci_path(a,method,rank),a.metric)
            scores[f"{rank}_{method}"]=d.POS.map(s).to_numpy(float)
    # Same variants for all 7 methods within a rank; no misleading differing N.
    shared={}
    for rank in RANKS:
        rank_cols=[f"{rank}_gene"]+[f"{rank}_{m}" for m in METHODS]
        ok=np.ones(len(d),dtype=bool)
        for col in rank_cols:ok&=np.isfinite(scores[col])
        shared[rank]=ok
        print(f"{rank}: shared cohort {ok.sum()} / {len(d)} unique variants")
        if ok.sum()<3:raise ValueError(f"{rank}: too few complete-case variants")
    # Reuse a single table containing NaN outside the rank-specific paired cohort.
    paired={}
    for rank in RANKS:
        for suffix in ("gene",)+METHODS:
            col=f"{rank}_{suffix}"
            paired[col]=np.where(shared[rank],scores[col],np.nan)
    plotted=pd.DataFrame(paired)
    summary=[]
    for method in METHODS:
        figure(plotted,method,out,a.metric)
        for rank in RANKS:
            ok=shared[rank]
            x=scores[f"{rank}_gene"][ok]
            y=scores[f"{rank}_{method}"][ok]
            r,rho=correlations(x,y)
            summary.append(dict(rank=rank,method=method,n_variants=int(ok.sum()),
                                pearson_r=r,spearman_rho=rho,
                                method_metric="original" if method=="old" else a.metric))
    table=pd.DataFrame(summary)
    old=table[table.method=="old"].set_index("rank")
    table["delta_pearson_vs_old"]=table.apply(lambda row:row.pearson_r-old.loc[row["rank"],"pearson_r"],axis=1)
    table["delta_spearman_vs_old"]=table.apply(lambda row:row.spearman_rho-old.loc[row["rank"],"spearman_rho"],axis=1)
    table.to_csv(out/"all_pipelines_gene_ci_correlation_summary.csv",index=False)
    pivot=table.pivot(index="rank",columns="method",values="pearson_r").reindex(RANKS)
    pivot.to_csv(out/"pearson_comparison_matrix.csv")
    fig,ax=plt.subplots(figsize=(14,5))
    xs=np.arange(len(RANKS))
    for method in METHODS:
        row=pivot[method].to_numpy()
        ax.plot(xs,row,marker="o",markersize=4,label=method.replace("_"," "))
    ax.set_xticks(xs,RANKS,rotation=30,ha="right")
    ax.set_ylabel("Pearson r with gene CI")
    ax.set_title("chrM CI vs gene CI: original, MITOS2 and GenBank methods")
    ax.grid(axis="y",alpha=.25)
    ax.legend(ncol=4,fontsize=8,loc="lower center",bbox_to_anchor=(.5,-.35))
    fig.tight_layout()
    fig.savefig(out/"pearson_comparison_all_methods.png",dpi=200,bbox_inches="tight")
    fig.savefig(out/"pearson_comparison_all_methods.pdf",bbox_inches="tight")
    plt.close(fig)
    print(f"Saved seven 10-panel figures and paired correlation summaries to {out}")

if __name__=="__main__":main()
