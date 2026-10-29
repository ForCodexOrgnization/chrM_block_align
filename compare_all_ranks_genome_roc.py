#!/usr/bin/env python3
"""Genome-wide ROC across all ranks, including kingdom, original vs block v6.

Requires all ten conservation CSVs. Uses an identical labeled variant cohort
for all ranks and scores, paired position-cluster bootstrap, and 0-100 scores.
Does not alter or rerun any alignments.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,roc_curve
from compare_old_new_ci import RANKS,read_old,read_new_rank
from analyze_ci_diagnostics import POSITIVE

def load(old_path,root,metric):
    d=read_old(old_path)
    clin=d.clinvar_interp_new.fillna("").astype(str).str.strip()
    mito=d.mitomap_status_new.fillna("").astype(str)
    pos=clin.isin(POSITIVE)|mito.str.contains("Cfrm",regex=False)
    neg=clin.eq("Benign")&~pos
    d=d.loc[pos|neg].copy()
    d["label"]=pos.loc[d.index].astype(np.int8)
    d=d.drop_duplicates(subset=["POS","REF","ALT","label"])
    for rank in RANKS:
        new,path=read_new_rank(root,rank,metric)
        if new is None:raise FileNotFoundError(f"Required rank absent: {path}")
        d=d.merge(new[["POS",f"{rank}_new_ci"]],on="POS",how="left",validate="many_to_one")
    score_columns=[f"{r}_{suffix}" for r in RANKS for suffix in ("rot_ref","new_ci")]
    for c in score_columns:
        d[c]=pd.to_numeric(d[c],errors="coerce")
    d=d.dropna(subset=score_columns).reset_index(drop=True)
    if not d[score_columns].stack().between(0,100).all():
        raise ValueError("CI outside 0..100")
    if d.label.nunique()!=2:raise ValueError("Shared cohort has fewer than two labels")
    return d

def run(d,nboot,seed,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    y=d.label.to_numpy(dtype=int)
    cluster,inverse=np.unique(d.POS.to_numpy(),return_inverse=True)
    rng=np.random.default_rng(seed)
    rows=[];draws=[]
    fig,axs=plt.subplots(2,5,figsize=(19,8),sharex=True,sharey=True)
    for rank,ax in zip(RANKS,axs.flat):
        old=d[f"{rank}_rot_ref"].to_numpy(dtype=float)
        new=d[f"{rank}_new_ci"].to_numpy(dtype=float)
        a=roc_auc_score(y,old);b=roc_auc_score(y,new)
        for scores,tag,aucv in ((old,"Original",a),(new,"Block v6",b)):
            fpr,tpr,_=roc_curve(y,scores)
            ax.plot(fpr,tpr,label=f"{tag}: {aucv:.3f}")
        ax.plot([0,1],[0,1],"k--",alpha=.3)
        ax.set_title(rank);ax.legend(fontsize=8)
        diff=[]
        for iteration in range(nboot):
            ids=rng.integers(len(cluster),size=len(cluster))
            w=np.bincount(ids,minlength=len(cluster))[inverse]
            if w[y==1].sum()==0 or w[y==0].sum()==0:continue
            oa=roc_auc_score(y,old,sample_weight=w)
            nb=roc_auc_score(y,new,sample_weight=w)
            diff.append(nb-oa)
            draws.append((rank,iteration,oa,nb,nb-oa))
        lo,hi=np.quantile(diff,[.025,.975]) if diff else (np.nan,np.nan)
        rows.append((rank,len(d),len(cluster),int(y.sum()),int((1-y).sum()),a,b,b-a,lo,hi,len(diff)))
    for i in range(5):axs[1,i].set_xlabel("False positive rate")
    for row in axs:row[0].set_ylabel("True positive rate")
    fig.suptitle("Genome-wide ROC, shared variants across all ten ranks")
    fig.tight_layout()
    fig.savefig(out/"all_ranks_roc.png",dpi=170);plt.close(fig)
    columns=["rank","n_variants","n_positions","n_pathogenic","n_benign",
             "auc_original","auc_block_v6","delta_auc","ci_lower","ci_upper","valid_bootstraps"]
    summary=pd.DataFrame(rows,columns=columns)
    summary.to_csv(out/"all_ranks_roc_summary.csv",index=False)
    pd.DataFrame(draws,columns=["rank","replicate","auc_original","auc_block_v6","delta_auc"]).to_csv(
        out/"all_ranks_roc_bootstrap.csv",index=False)
    print(summary.to_string(index=False))

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--old-csv",default="data/2026_10_10_mito_cons_anns_main.csv")
    p.add_argument("--results-dir",default="block_results_v6")
    p.add_argument("--metric",choices=("CI_MSA","CI_callable"),default="CI_MSA")
    p.add_argument("--out-dir",default="block_results_v6/all_ranks_genome_roc")
    p.add_argument("--bootstrap-reps",type=int,default=1000)
    p.add_argument("--seed",type=int,default=20261029)
    args=p.parse_args()
    if args.bootstrap_reps<1:p.error("bootstrap-reps must be positive")
    out=Path(args.out_dir);out.mkdir(parents=True,exist_ok=True)
    d=load(args.old_csv,args.results_dir,args.metric)
    d.to_csv(out/"shared_labeled_variants.csv",index=False)
    run(d,args.bootstrap_reps,args.seed,out)
if __name__=="__main__":main()
