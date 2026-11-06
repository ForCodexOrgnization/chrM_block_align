#!/usr/bin/env python3
"""Four figures: strict OH containment by longest eligible interval, 10 ranks.

For each of two gene-free definitions: (1) binned bar chart and (2) CDF.
Each species contributes AT MOST ONE interval (the longest meeting the full
OriH-containment criterion). Each percentage denominator is ALL species in
that rank, including missing-OriH and no-candidate species. No threshold lines.
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
MODES=("no_human_homologous_genes","no_genes_at_all")
BINS=[0,10,50,100,200,300,500,1000,1500,2500,5000,np.inf]
LABELS=["1-10","11-50","51-100","101-200","201-300","301-500",
        "501-1000","1001-1500","1501-2500","2501-5000",">5000"]
TICKS=[0,2,4,6,8,10]
BLUE="#327FA6"
ORANGE="#C78637"

def figures(df,out,mode):
    part=df[df["mode"]==mode].copy()
    if part.empty:raise ValueError("No records for "+mode)
    part["selected_gap_bp"]=pd.to_numeric(part["selected_gap_bp"],errors="coerce")
    selected=part[part["selected"]==1].copy()
    if selected.selected_gap_bp.isna().any() or (selected.selected_gap_bp<=0).any():
        raise ValueError("Selected intervals must have positive lengths")
    selected["bin"]=pd.cut(selected.selected_gap_bp,bins=BINS,labels=LABELS,
                           include_lowest=True)
    if selected["bin"].isna().any():
        raise ValueError("Unbinned selected intervals")
    bins=[];overall=[]
    for rank in RANKS:
        group=part[part["rank"]==rank]
        n=len(group)
        chosen=selected[selected["rank"]==rank]
        counts=chosen["bin"].value_counts().reindex(LABELS,fill_value=0)
        running=0
        for label,count in counts.items():
            running+=int(count)
            bins.append(dict(rank=rank,mode=mode,length_bin=label,
                n_total_species=n,n_selected_in_bin=int(count),
                pct_all_species_in_bin=100*count/n if n else np.nan,
                cumulative_n_selected=running,
                cumulative_pct_all_species=100*running/n if n else np.nan))
        overall.append(dict(rank=rank,mode=mode,n_total_species=n,
            n_selected_full_containment=len(chosen),
            pct_all_species_selected=100*len(chosen)/n if n else np.nan,
            n_missing_OriH=int((group["n_OriH_features"]==0).sum())))
    table=pd.DataFrame(bins)
    table.to_csv(out/(mode+"_length_bin_counts_and_cdf.tsv"),sep="\t",index=False)
    pd.DataFrame(overall).to_csv(out/(mode+"_rank_summary.tsv"),sep="\t",index=False)

    # Bin percentage: denominator remains all rank species, NOT bin total.
    fig,axes=plt.subplots(2,5,figsize=(23,10),sharey=True)
    for i,(ax,rank) in enumerate(zip(axes.flat,RANKS)):
        r=table[table["rank"]==rank]
        x=np.arange(len(LABELS))
        ax.bar(x,r["pct_all_species_in_bin"].to_numpy(),color=BLUE,width=.82)
        ax.set_title(f"{rank} (N={r.n_total_species.iloc[0]:,})",fontsize=11,pad=8)
        ax.set_ylim(0,100);ax.set_yticks([0,20,40,60,80,100])
        ax.set_xticks(TICKS,[LABELS[j] for j in TICKS],
                      rotation=45,ha="right",fontsize=8)
        ax.grid(axis="y",alpha=.16);ax.set_axisbelow(True)
        if i%5==0:ax.set_ylabel("Species in length bin (% of all species)",fontsize=10)
    fig.suptitle(f"Fully OriH-containing intervals by length: {mode.replace('_',' ')}",
                 fontsize=15,y=.98)
    fig.subplots_adjust(left=.055,right=.99,bottom=.11,top=.91,wspace=.20,hspace=.48)
    for ext in ("png","pdf"):
        fig.savefig(out/(mode+"_length_bin_barplot."+ext),dpi=220)
    plt.close(fig)

    # CDF uses exact same fixed rank denominator across every length bin.
    fig,axes=plt.subplots(2,5,figsize=(23,10),sharey=True)
    for i,(ax,rank) in enumerate(zip(axes.flat,RANKS)):
        r=table[table["rank"]==rank]
        x=np.arange(len(LABELS))
        y=r["cumulative_pct_all_species"].to_numpy()
        if np.any(np.diff(y)<-1e-9):
            raise AssertionError(f"CDF must be nondecreasing: {rank}")
        ax.step(x,y,where="post",lw=1.8,color=ORANGE)
        ax.scatter(x,y,s=14,color=ORANGE,zorder=3)
        ax.set_title(f"{rank} (N={r.n_total_species.iloc[0]:,})",fontsize=11,pad=8)
        ax.set_ylim(0,100);ax.set_yticks([0,20,40,60,80,100])
        ax.set_xlim(-.15,len(LABELS)-.85)
        ax.set_xticks(TICKS,[LABELS[j] for j in TICKS],
                      rotation=45,ha="right",fontsize=8)
        ax.grid(axis="y",alpha=.16);ax.set_axisbelow(True)
        if i%5==0:ax.set_ylabel("Cumulative % of all species",fontsize=11)
    fig.suptitle(f"Cumulative full OriH containment: {mode.replace('_',' ')}",
                 fontsize=15,y=.98)
    fig.subplots_adjust(left=.055,right=.99,bottom=.11,top=.91,wspace=.20,hspace=.48)
    for ext in ("png","pdf"):
        fig.savefig(out/(mode+"_cumulative_CDF."+ext),dpi=220)
    plt.close(fig)
    return pd.DataFrame(overall)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",default="strict_orih_gene_desert_qc/all_ranks_species_strict_OriH.tsv")
    p.add_argument("--out-dir",default="strict_orih_gene_desert_qc/figures")
    a=p.parse_args();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(a.input,sep="\t")
    required={"rank","species","mode","selected","selected_gap_bp","n_OriH_features"}
    if not required.issubset(df.columns):
        raise ValueError(f"Missing columns: {required-set(df.columns)}")
    if df.duplicated(["rank","species","mode"]).any():
        raise ValueError("Expected at most one record per species/rank/mode")
    summaries=[]
    for mode in MODES:
        summaries.append(figures(df,out,mode))
    summary=pd.concat(summaries,ignore_index=True)
    # A no-gene-at-all gap is also a no-human-homologous-gene gap.
    wide=summary.pivot(index="rank",columns="mode",values="n_selected_full_containment")
    if (wide["no_genes_at_all"]>wide["no_human_homologous_genes"]).any():
        raise AssertionError("All-gene exclusion cannot select more species than canonical-only exclusion")
    summary.to_csv(out/"both_modes_rank_summary.tsv",sep="\t",index=False)
    print("Four figures (PNG/PDF), binned counts, and summaries saved to",out)

if __name__=="__main__":main()
