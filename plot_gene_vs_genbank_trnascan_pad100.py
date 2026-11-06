#!/usr/bin/env python3
"""Gene CI vs original and GenBank+tRNAscan pad100 chrM CI, by rank.

One data point per human reference position (1..16569), rather than per variant.
The gene-CI baseline is <rank>_ref_gene from the original annotation CSV.
Feature categories are mapped from the 37 genes in human GenBank. Where genes
overlap, color precedence is tRNA > rRNA > protein > non-coding.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RANKS=("kingdom","phylum","subphylum","superclass","class","superorder",
       "order","suborder","infraorder","parvorder")
COLOR={"Non-coding":"#F5C400","Protein":"#377EB8","rRNA":"#984EA3","tRNA":"#FF7F00"}

def categories_from_genbank(path):
    from Bio import SeqIO
    rec=SeqIO.read(path,"genbank")
    if len(rec.seq)!=16569:raise ValueError("Expected human NC_012920 length 16569")
    labels=np.array(["Non-coding"]*16569,dtype=object)
    for feature_type,label in (("CDS","Protein"),("rRNA","rRNA"),("tRNA","tRNA")):
        for f in rec.features:
            if f.type!=feature_type:continue
            for pos in f.location:
                labels[int(pos)]=label
    return pd.Series(labels,index=np.arange(1,16570),name="category")

def old_scores(path):
    usecols=["POS"]+[f"{r}_{s}" for r in RANKS for s in ("ref_gene","rot_ref")]
    d=pd.read_csv(path,usecols=usecols,low_memory=False)
    d["POS"]=pd.to_numeric(d["POS"],errors="raise").astype(int)
    if not d.POS.between(1,16569).all():raise ValueError("POS out of bounds")
    for col in usecols[1:]: d[col]=pd.to_numeric(d[col],errors="coerce")
    # A position may have several REF/ALT records; use the per-position median.
    return d.groupby("POS",sort=True)[usecols[1:]].median()

def new_scores(path,metric):
    d=pd.read_csv(path,usecols=["Position",metric])
    d["Position"]=pd.to_numeric(d["Position"],errors="raise").astype(int)
    if len(d)!=16569 or not d.Position.tolist()==list(range(1,16570)):
        raise ValueError(f"Expected 16569 ordered positions: {path}")
    score=pd.to_numeric(d[metric],errors="coerce")
    if not score.dropna().between(0,100).all():
        raise ValueError(f"CI outside 0..100 in {path}")
    return pd.Series(score.to_numpy(),index=d.Position.to_numpy(),name="new")

def scatter(ax,x,y,categories,name):
    mask=np.isfinite(x)&np.isfinite(y)
    x=x[mask];y=y[mask];categories=categories[mask]
    for label in ("Non-coding","Protein","rRNA","tRNA"):
        subset=categories==label
        if subset.any():
            ax.scatter(x[subset],y[subset],s=5,alpha=.38,color=COLOR[label],
                       rasterized=True,lw=0)
    ax.plot([0,100],[0,100],ls="--",lw=.7,color="black",alpha=.7)
    ax.set_xlim(0,100);ax.set_ylim(0,100)
    ax.set_title(name,fontsize=9,fontweight="bold")
    if len(x)>2 and np.std(x)>0 and np.std(y)>0:
        pr=pd.Series(x).corr(pd.Series(y),method="pearson")
        sr=pd.Series(x).corr(pd.Series(y),method="spearman")
    else:pr=sr=float("nan")
    ax.text(.035,.965,f"Pearson r = {pr:.3f}\nSpearman ρ = {sr:.3f}\nn = {len(x):,}",
            transform=ax.transAxes,ha="left",va="top",fontsize=7,
            bbox={"facecolor":"white","edgecolor":"none","alpha":.78})
    return len(x),pr,sr

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--annotation-csv",default="data/2026_10_10_mito_cons_anns_main.csv")
    ap.add_argument("--human-genbank",default="data/sequence.gb")
    ap.add_argument("--genbank-root",default="block_results_genbank_trnascan_v1")
    ap.add_argument("--metric",choices=["CI_MSA","CI_callable"],default="CI_MSA")
    ap.add_argument("--out-dir",default="gene_vs_genbank_trnascan_pad100")
    ap.add_argument("--new-only",action="store_true",help="Show 10 new-only panels (2x5), rather than 10 old/new pairs (5x4)")
    a=ap.parse_args()
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    old=old_scores(a.annotation_csv)
    cat=categories_from_genbank(a.human_genbank)
    if a.new_only:
        fig,axes=plt.subplots(2,5,figsize=(19,8),sharex=True,sharey=True)
    else:
        fig,axes=plt.subplots(5,4,figsize=(21,17),sharex=True,sharey=True)
    rows=[]
    for i,rank in enumerate(RANKS):
        path=Path(a.genbank_root)/"pad100"/rank/"conservation"/"conservation_ci.csv"
        if not path.is_file():raise FileNotFoundError(path)
        new=new_scores(path,a.metric)
        gene=old[f"{rank}_ref_gene"].reindex(cat.index)
        orig=old[f"{rank}_rot_ref"].reindex(cat.index)
        # Identical paired cohort for old and new panels.
        paired=np.isfinite(gene.to_numpy()) & np.isfinite(orig.to_numpy()) & np.isfinite(new.to_numpy())
        c=cat.to_numpy()[paired]
        y=gene.to_numpy()[paired]
        if a.new_only:
            ax=axes.flat[i]
            n,r,rho=scatter(ax,new.to_numpy()[paired],y,c,rank+" (new pad100)")
            rows.append(dict(rank=rank,method="genbank_trnascan_pad100",
                             n_positions=n,pearson_r=r,spearman_rho=rho))
        else:
            ax_old,ax_new=axes.flat[2*i:2*i+2]
            n,r,rho=scatter(ax_old,orig.to_numpy()[paired],y,c,rank+" (old)")
            rows.append(dict(rank=rank,method="old",n_positions=n,
                             pearson_r=r,spearman_rho=rho))
            n,r,rho=scatter(ax_new,new.to_numpy()[paired],y,c,rank+" (new pad100)")
            rows.append(dict(rank=rank,method="genbank_trnascan_pad100",
                             n_positions=n,pearson_r=r,spearman_rho=rho))
    for ax in axes.flat:
        ax.set_xlabel("chrM CI",fontsize=8)
        ax.set_ylabel("gene CI",fontsize=8)
        ax.tick_params(labelsize=7)
    handles=[Line2D([0],[0],marker="o",linestyle="",markersize=6,
                    color=COLOR[k],label=k) for k in COLOR]
    fig.legend(handles=handles,loc="upper center",ncol=4,frameon=False,
               bbox_to_anchor=(.5,.966),fontsize=9)
    fig.suptitle("Gene vs chrM CI scores: "+("GENBANK+tRNAscan pad100" if a.new_only
                  else "OLD vs GENBANK+tRNAscan NEW (pad100)"),
                 fontsize=13,fontweight="bold",y=.995)
    fig.tight_layout(rect=[0,0,1,.944])
    stem="gene_vs_genbank_trnascan_pad100"+("_new_only" if a.new_only else "_old_vs_new")
    for suffix in ("png","pdf"):
        fig.savefig(out/f"{stem}.{suffix}",dpi=200)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(out/f"{stem}_correlations.csv",index=False)
    print("Saved",out/stem, "(PNG, PDF and correlations CSV)")

if __name__=="__main__":
    main()
