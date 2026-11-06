#!/usr/bin/env python3
"""Plot longest all-gene-free candidate interval (one per species) for 10 ranks.
No length cutoff in selection, and no vertical threshold reference lines.
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
EDGES=[-0.5,0,10,50,100,200,300,500,1000,1500,2500,5000,float("inf")]
LABELS=["0","1-10","11-50","51-100","101-200","201-300","301-500",
        "501-1000","1001-1500","1501-2500","2501-5000",">5000"]

def save(fig,out,stem):
    fig.tight_layout(rect=(0,0,1,.95))
    for ext in ("png","pdf"):fig.savefig(out/(stem+"."+ext),dpi=220)
    plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",default="longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv")
    p.add_argument("--out-dir",default="longest_all_gene_desert_qc/figures")
    a=p.parse_args();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(a.input,sep="\t")
    # Fixed CDF denominators are ALL species in each rank, including species
    # without an OriH call and species without a definable gene-free interval.
    if d.duplicated(["rank","species"]).any():
        raise ValueError("More than one input record per rank/species")
    rank_totals=d.groupby("rank").size().to_dict()
    d=d[d.candidate_defined==1].copy()
    if d.empty:raise ValueError("No defined intervals")
    d["longest_gap_bp"]=pd.to_numeric(d.longest_gap_bp,errors="coerce")
    d["OriH_overlaps"]=pd.to_numeric(d.OriH_overlaps,errors="coerce")
    d["OriH_fully_contained"]=pd.to_numeric(d.OriH_fully_contained,errors="coerce")
    d=d.dropna(subset=["longest_gap_bp"])
    if d.duplicated(["rank","species"]).any():raise ValueError("More than one candidate per rank/species")
    d["bin"]=pd.cut(d.longest_gap_bp,bins=EDGES,labels=LABELS)
    bins=[]
    for rank in RANKS:
        r=d[d["rank"]==rank]
        for label in LABELS:
            x=r[r["bin"]==label]
            known=x[x["n_OriH_features"].fillna(0)>0]
            n=len(known)
            overlapping=int((known.OriH_overlaps==1).sum())
            full=int((known.OriH_fully_contained==1).sum())
            bins.append(dict(rank=rank,length_bin=label,n_species=len(x),
                 n_with_OriH_prediction=n,n_overlapping=overlapping,
                 n_fully_contained=full,
                 pct_overlapping=100*overlapping/n if n else np.nan,
                 pct_fully_contained=100*full/n if n else np.nan))
    data=pd.DataFrame(bins)
    data.to_csv(out/"longest_gene_free_interval_length_bins.tsv",sep="\t",index=False)
    fig,axes=plt.subplots(2,5,figsize=(23,9))
    for ax,rank in zip(axes.flat,RANKS):
        r=data[data["rank"]==rank]
        x=np.arange(len(LABELS))
        ax.bar(x,r.n_species,color="#658EAD",width=.9)
        ax.set_title(f"{rank} (n={r.n_species.sum():,})")
        ax.set_xticks(x,LABELS,rotation=60,ha="right",fontsize=7)
        ax.set_yscale("symlog",linthresh=1)
        ax.grid(axis="y",alpha=.15);ax.set_axisbelow(True)
    fig.suptitle("Longest annotated gene-free interval per species: length distribution")
    save(fig,out,"longest_gap_length_distribution")
    fig,axes=plt.subplots(2,5,figsize=(23,10),sharey=True)
    # Show representative length bins in each panel; data points for ALL bins
    # are still plotted. Avoid dense x labels and overlapping shared y labels.
    tick_idx=[0,2,4,6,8,11]
    for i,(ax,rank) in enumerate(zip(axes.flat,RANKS)):
        r=data[data["rank"]==rank];x=np.arange(len(LABELS))
        ax.plot(x,r.pct_overlapping,marker="o",ms=3,lw=1.6,color="#2D759E",label=">=1 bp overlap")
        ax.plot(x,r.pct_fully_contained,marker="s",ms=3,lw=1.6,color="#B77B34",label="100% containment")
        ax.set_title(rank,fontsize=12,pad=8)
        ax.set_xlim(-.3,len(LABELS)-.7);ax.set_ylim(0,103)
        ax.set_yticks([0,20,40,60,80,100])
        ax.set_xticks(tick_idx)
        ax.set_xticklabels([LABELS[j] for j in tick_idx],
                           rotation=45,ha="right",fontsize=8)
        ax.tick_params(axis="y",labelsize=9)
        ax.grid(axis="y",alpha=.15)
        ax.set_axisbelow(True)
        if i%5==0:
            ax.set_ylabel("Species with OriH support (%)",fontsize=11,labelpad=8)
    # Use explicit margins instead of tight_layout, which pushed the long
    # figure-level y label into the middle of the two panel rows.
    handles,legend_labels=axes.flat[0].get_legend_handles_labels()
    fig.suptitle("OriH validation of the longest gene-free interval (per rank)",
                 fontsize=15,y=.985)
    fig.legend(handles,legend_labels,ncol=2,loc="upper center",
               bbox_to_anchor=(.5,.952),frameon=False,fontsize=10)
    fig.subplots_adjust(left=.06,right=.99,bottom=.10,top=.88,
                        wspace=.18,hspace=.47)
    for ext in ("png","pdf"):
        fig.savefig(out/("longest_gap_OriH_overlap_by_length."+ext),
                    dpi=220)
    plt.close(fig)
    # Genuine cumulative distribution, with a fixed denominator equal to
    # EVERY species in the rank FASTA, including those lacking an OriH call
    # or a valid longest gene-free interval.
    #
    # Gray: cumulative fraction with a definable longest gap of length <= L.
    # Blue: cumulative fraction of ALL species with >=1 bp OriH overlap.
    # Orange: cumulative fraction of ALL species with 100% OriH containment.
    cumulative=[]
    for rank in RANKS:
        r=data[data["rank"]==rank].reset_index(drop=True)
        denom=int(rank_totals.get(rank,0))
        if denom==0:continue
        n_candidates=r["n_species"].cumsum()
        n_with_oriH=r["n_with_OriH_prediction"].cumsum()
        hits=r["n_overlapping"].cumsum()
        full=r["n_fully_contained"].cumsum()
        for j,label in enumerate(LABELS):
            cumulative.append(dict(
                rank=rank,length_bin=label,n_total_species=denom,
                cumulative_n_candidates=int(n_candidates.iloc[j]),
                cumulative_n_with_OriH_prediction=int(n_with_oriH.iloc[j]),
                cumulative_n_overlapping=int(hits.iloc[j]),
                cumulative_n_fully_contained=int(full.iloc[j]),
                cumulative_pct_candidates=100*n_candidates.iloc[j]/denom,
                cumulative_pct_overlap_ge1bp=100*hits.iloc[j]/denom,
                cumulative_pct_fully_contained=100*full.iloc[j]/denom))
    cumulative=pd.DataFrame(cumulative)
    cumulative.to_csv(out/"longest_gap_OriH_cumulative_by_length.tsv",
                      sep="\t",index=False)
    fig,axes=plt.subplots(2,5,figsize=(23,10),sharey=True)
    for i,(ax,rank) in enumerate(zip(axes.flat,RANKS)):
        r=cumulative[cumulative["rank"]==rank]
        x=np.arange(len(r))
        ax.plot(x,r.cumulative_pct_candidates,marker=".",ms=4,lw=1.4,
                color="#888888",ls="--",label="Species with a defined interval")
        ax.plot(x,r.cumulative_pct_overlap_ge1bp,marker="o",ms=3,lw=1.6,
                color="#2D759E",label=">=1 bp OriH overlap")
        ax.plot(x,r.cumulative_pct_fully_contained,marker="s",ms=3,lw=1.6,
                color="#B77B34",label="100% OriH containment")
        ax.set_title(f"{rank} (N={rank_totals.get(rank,0):,})",
                     fontsize=12,pad=8)
        ax.set_xlim(-.3,len(LABELS)-.7);ax.set_ylim(0,103)
        ax.set_yticks([0,20,40,60,80,100])
        ax.set_xticks(tick_idx)
        ax.set_xticklabels([LABELS[j] for j in tick_idx],
                           rotation=45,ha="right",fontsize=8)
        ax.tick_params(axis="y",labelsize=9)
        ax.grid(axis="y",alpha=.15)
        ax.set_axisbelow(True)
        if i%5==0:
            ax.set_ylabel("Cumulative % of all species in rank",
                          fontsize=11,labelpad=8)
    handles,legend_labels=axes.flat[0].get_legend_handles_labels()
    fig.suptitle("Cumulative longest gene-free interval and OriH support (all species)",
                 fontsize=15,y=.985)
    fig.legend(handles,legend_labels,ncol=3,loc="upper center",
               bbox_to_anchor=(.5,.952),frameon=False,fontsize=10)
    fig.subplots_adjust(left=.06,right=.99,bottom=.10,top=.88,
                        wspace=.18,hspace=.47)
    for ext in ("png","pdf"):
        fig.savefig(out/("longest_gap_OriH_cumulative_by_length."+ext),
                    dpi=220)
    plt.close(fig)
    # Species-level threshold sensitivity: one longest candidate per species,
    # denominators require an available OriH prediction.
    curves=[]
    for rank in RANKS:
        r=d[(d["rank"]==rank)&(d["n_OriH_features"].fillna(0)>0)]
        for threshold in np.sort(r.longest_gap_bp.unique()):
            s=r[r.longest_gap_bp>=threshold]
            curves.append(dict(rank=rank,min_length_bp=threshold,n_species=len(s),
                pct_OriH_overlap=100*(s.OriH_overlaps==1).mean(),
                pct_OriH_contained=100*(s.OriH_fully_contained==1).mean()))
    curve=pd.DataFrame(curves)
    curve.to_csv(out/"minimum_length_retention_curves.tsv",sep="\t",index=False)
    fig,axes=plt.subplots(2,5,figsize=(22,9),sharey=True)
    for ax,rank in zip(axes.flat,RANKS):
        r=curve[curve["rank"]==rank]
        if not r.empty:
            ax.plot(r.min_length_bp,r.pct_OriH_overlap,color="#2D759E",lw=1.6)
            ax.plot(r.min_length_bp,r.pct_OriH_contained,color="#B77B34",lw=1.6)
        ax.set_xscale("symlog",linthresh=10);ax.set_ylim(0,103)
        ax.set_title(rank);ax.set_xlabel("Minimum longest-gap length (bp)",fontsize=9)
        ax.grid(alpha=.15)
    fig.suptitle("OriH validation rate when retaining longer gene-free candidates")
    save(fig,out,"minimum_length_OriH_rate")
    print("Saved four figures (PNG and PDF) and TSV summaries to",out)

if __name__=="__main__":main()
