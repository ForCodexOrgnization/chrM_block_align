#!/usr/bin/env python3
"""Genome-wide ROC for original, pad50 and pad100 CI across all ten ranks.

Shared complete-case ClinVar/MITOMAP labeled variants; paired position-cluster
bootstrap for differences. No alignment or conservation files are changed.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,roc_curve
from compare_old_new_ci import RANKS,read_old,read_new_rank
from analyze_ci_diagnostics import POSITIVE

METHODS=("original","pad50","pad100")
COMPARISONS=(("pad50","original"),("pad100","original"),("pad100","pad50"))

def load(old_csv,padded_dir,metric):
    d=read_old(old_csv)
    clin=d.clinvar_interp_new.fillna("").astype(str).str.strip()
    mito=d.mitomap_status_new.fillna("").astype(str)
    pos=clin.isin(POSITIVE)|mito.str.contains("Cfrm",regex=False)
    neg=clin.eq("Benign")&~pos
    d=d.loc[pos|neg].copy()
    d["label"]=pos.loc[d.index].astype(int)
    d=d.drop_duplicates(subset=["POS","REF","ALT","label"])
    for rank in RANKS:
        d=d.rename(columns={f"{rank}_rot_ref":f"{rank}_original"})
        for pad in ("pad50","pad100"):
            path=Path(padded_dir)/pad/rank/"conservation"/"conservation_ci.csv"
            if not path.exists():raise FileNotFoundError(path)
            fresh=pd.read_csv(path,usecols=["Position",metric])
            if len(fresh)!=16569 or fresh.Position.nunique()!=16569 or not fresh.Position.between(1,16569).all():
                raise ValueError(f"Invalid 16569-position CI table: {path}")
            fresh=fresh.rename(columns={"Position":"POS",metric:f"{rank}_{pad}"})
            d=d.merge(fresh,on="POS",how="left",validate="many_to_one")
    cols=[f"{rank}_{method}" for rank in RANKS for method in METHODS]
    for col in cols:d[col]=pd.to_numeric(d[col],errors="coerce")
    d=d.dropna(subset=cols).reset_index(drop=True)
    if not d[cols].stack().between(0,100).all():raise ValueError("CI must be 0-100")
    if d.label.nunique()!=2:raise ValueError("No positive and negative labels after filtering")
    return d

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--old-csv",default="data/2026_10_10_mito_cons_anns_main.csv")
    p.add_argument("--padded-dir",default="block_results_padded_v1")
    p.add_argument("--metric",choices=["CI_MSA","CI_callable"],default="CI_MSA")
    p.add_argument("--bootstrap-reps",type=int,default=1000)
    p.add_argument("--seed",type=int,default=20261030)
    p.add_argument("--out-dir",default="block_results_padded_v1/all_ranks_padding_roc")
    a=p.parse_args()
    if a.bootstrap_reps<1:p.error("bootstrap-reps must be positive")
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    d=load(a.old_csv,a.padded_dir,a.metric)
    d.to_csv(out/"shared_labeled_variants.csv",index=False)
    y=d.label.to_numpy()
    pos,inv=np.unique(d.POS.to_numpy(),return_inverse=True)
    rng=np.random.default_rng(a.seed)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,5,figsize=(20,8),sharex=True,sharey=True)
    rows=[];boot=[]
    for rank,ax in zip(RANKS,axes.flat):
        scores={m:d[f"{rank}_{m}"].to_numpy(dtype=float) for m in METHODS}
        point={m:roc_auc_score(y,scores[m]) for m in METHODS}
        for m in METHODS:
            fpr,tpr,_=roc_curve(y,scores[m])
            ax.plot(fpr,tpr,label=f"{m}: {point[m]:.3f}")
        ax.plot([0,1],[0,1],"k--",alpha=.3)
        ax.set_title(rank);ax.legend(fontsize=8)
        diffs={pair:[] for pair in COMPARISONS}
        for iteration in range(a.bootstrap_reps):
            sample=rng.integers(len(pos),size=len(pos))
            weights=np.bincount(sample,minlength=len(pos))[inv]
            if weights[y==0].sum()==0 or weights[y==1].sum()==0:continue
            values={m:roc_auc_score(y,scores[m],sample_weight=weights) for m in METHODS}
            boot.append({"rank":rank,"replicate":iteration,**{f"auc_{m}":values[m] for m in METHODS}})
            for pair in COMPARISONS:
                b,base=pair;diffs[pair].append(values[b]-values[base])
        row={"rank":rank,"n_variants":len(d),"n_positions":len(pos),
             "n_pathogenic":int(y.sum()),"n_benign":int((1-y).sum()),
             **{f"auc_{m}":point[m] for m in METHODS}}
        for pair in COMPARISONS:
            b,base=pair;tag=f"{b}_minus_{base}"
            diff=np.asarray(diffs[pair])
            row[f"delta_{tag}"]=point[b]-point[base]
            row[f"ci_lower_{tag}"]=float(np.quantile(diff,.025)) if len(diff) else np.nan
            row[f"ci_upper_{tag}"]=float(np.quantile(diff,.975)) if len(diff) else np.nan
        row["valid_bootstrap_replicates"]=len(diffs[COMPARISONS[0]])
        rows.append(row)
    for i in range(5):axes[1,i].set_xlabel("False positive rate")
    for i in range(2):axes[i,0].set_ylabel("True positive rate")
    fig.suptitle("Genome-wide ROC: original vs pad50 vs pad100 (shared variants)")
    fig.tight_layout();fig.savefig(out/"all_ranks_padding_roc.png",dpi=180);plt.close(fig)
    pd.DataFrame(rows).to_csv(out/"all_ranks_padding_roc_summary.csv",index=False)
    pd.DataFrame(boot).to_csv(out/"all_ranks_padding_roc_bootstrap.csv",index=False)
    print(pd.DataFrame(rows).to_string(index=False))

if __name__=="__main__":main()
