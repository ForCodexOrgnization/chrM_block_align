#!/usr/bin/env python3
"""Load new rank-specific stitched-MSA CI and compare it with old chrM CI.

Old input: 2026_10_10_mito_cons_anns_main.csv (POS, {rank}_rot_ref).
New input: block_results_v6/{rank}/conservation/conservation_ci.csv
( Position, Homo sapiens, CI_MSA, CI_callable, n_callable, ... ).

Old scores are 0-1, new CI values are percentages; new values are /100.
One row per human position is used for score/QC summaries. The variant
table retains all REF/ALT rows for ROC and consequence comparisons.
"""
import argparse
import json
import warnings
from pathlib import Path
import numpy as np
import pandas as pd

RANKS=("kingdom","phylum","subphylum","superclass","class",
       "superorder","order","suborder","infraorder","parvorder")
METRICS=("CI_MSA","CI_callable")
CORE=["POS","REF","ALT","symbol","consequence","clinvar_interp_new","mitomap_status_new"]

def read_old(path):
    df=pd.read_csv(path,low_memory=False)
    required={"POS",*(f"{r}_rot_ref" for r in RANKS)}
    missing=required-set(df.columns)
    if missing:raise ValueError(f"Old annotations missing columns: {sorted(missing)}")
    df["POS"]=pd.to_numeric(df["POS"],errors="coerce")
    if df.POS.isna().any() or ((df.POS<1)|(df.POS>16569)).any():
        raise ValueError("Old file has missing/out-of-range POS")
    df.POS=df.POS.astype(int)
    for rank in RANKS:
        col=f"{rank}_rot_ref"
        df[col]=pd.to_numeric(df[col],errors="coerce")
        valid=df[col].dropna()
        if len(valid) and not valid.between(0,1).all():
            raise ValueError(f"Old {col} must be on a 0-1 scale")
    return df

def read_new_rank(root,rank,metric):
    path=Path(root)/rank/"conservation"/"conservation_ci.csv"
    if not path.exists():return None,path
    cols=["Position","Homo sapiens","CI_MSA","CI_callable","n_callable","n_species","n_gap","n_ambiguous"]
    df=pd.read_csv(path)
    missing=set(cols)-set(df.columns)
    if missing:raise ValueError(f"{path}: missing {sorted(missing)}")
    if len(df)!=16569:raise ValueError(f"{path}: {len(df)} rows, expected 16569")
    df["Position"]=pd.to_numeric(df["Position"],errors="raise").astype(int)
    if df.Position.tolist()!=list(range(1,16570)):
        raise ValueError(f"{path}: coordinates must be unique and ordered 1..16569")
    for col in METRICS:
        df[col]=pd.to_numeric(df[col],errors="coerce")
        if not df[col].dropna().between(0,100).all():
            raise ValueError(f"{path}: {col} must be 0-100 percent")
    df=df.rename(columns={"Position":"POS","Homo sapiens":f"{rank}_human",
                          "n_callable":f"{rank}_n_callable","n_species":f"{rank}_n_species",
                          "n_gap":f"{rank}_n_gap","n_ambiguous":f"{rank}_n_ambiguous"})
    df[f"{rank}_new_ci"]=df[metric]/100.0
    df[f"{rank}_new_ci_msa"]=df["CI_MSA"]/100.0
    df[f"{rank}_new_ci_callable"]=df["CI_callable"]/100.0
    return df.drop(columns=["CI_MSA","CI_callable"]),path

def load_comparison(old_csv,results_dir,metric="CI_MSA",strict=False):
    old=read_old(old_csv)
    keep=["POS"]+[c for c in CORE if c in old and c!="POS"]
    # Avoid cartesian joins: new has one row per position, old retains variants.
    variant=old.copy()
    per_position=old.groupby("POS",sort=True,as_index=False).first()
    report=[]
    loaded=[]
    for rank in RANKS:
        fresh,path=read_new_rank(results_dir,rank,metric)
        if fresh is None:
            print(f"[PENDING] {rank}: {path}")
            report.append({"rank":rank,"status":"pending"})
            continue
        loaded.append(rank)
        variant=variant.merge(fresh,on="POS",how="left",validate="many_to_one")
        per_position=per_position.merge(fresh,on="POS",how="left",validate="one_to_one")
        old_col=f"{rank}_rot_ref";new_col=f"{rank}_new_ci"
        sub=per_position[[old_col,new_col]].dropna()
        if not len(sub):raise ValueError(f"No paired scores for {rank}")
        mismatches=0
        if "REF" in per_position:
            ref=per_position["REF"].astype(str).str.upper()
            human=per_position[f"{rank}_human"].astype(str).str.upper()
            mismatches=int(((ref!=human)&ref.isin(list("ACGT"))).sum())
            if mismatches:
                warnings.warn(f"{rank}: {mismatches} old REF bases differ from new human reference")
        zero_frac=float((sub[old_col]==0).mean())
        if zero_frac>.5:
            warnings.warn(f"{rank}: {zero_frac:.1%} of old CI values are zero; "
                          "confirm whether zero means unavailable before interpretation")
        report.append({"rank":rank,"status":"loaded","path":str(path),
                       "n_paired_positions":len(sub),
                       "old_zero_fraction":zero_frac,
                       "old_mean":float(sub[old_col].mean()),
                       "new_mean":float(sub[new_col].mean()),
                       "mean_delta_new_minus_old":float((sub[new_col]-sub[old_col]).mean()),
                       "pearson_r":float(sub.corr().iloc[0,1]) if sub[old_col].std()>0 and sub[new_col].std()>0 else None,
                       "reference_mismatch_positions":mismatches})
        print(f"[LOADED] {rank}: {len(fresh)} positions, {len(sub)} old/new pairs")
    if strict and len(loaded)!=len(RANKS):raise ValueError(f"Ranks still pending: {set(RANKS)-set(loaded)}")
    return variant,per_position,pd.DataFrame(report),loaded

def labels(df):
    clin=df.get("clinvar_interp_new",pd.Series("",index=df.index)).fillna("").astype(str).str.strip()
    mitomap=df.get("mitomap_status_new",pd.Series("",index=df.index)).fillna("").astype(str)
    yes=clin.isin(["Pathogenic","Pathogenic/Likely pathogenic","Likely pathogenic"])|mitomap.str.contains("Cfrm",regex=False)
    no=clin.eq("Benign")&~yes
    return pd.Series(np.where(yes,1,np.where(no,0,np.nan)),index=df.index)

def write_plots(variant,per_position,ranks,out,do_roc=True):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_auc_score,roc_curve
    # Per-position delta distributions, never inflated by 3 ALT variants at a site.
    if ranks:
        fig,ax=plt.subplots(figsize=(11,5))
        data=[(per_position[f"{r}_new_ci"]-per_position[f"{r}_rot_ref"]).dropna().values
              for r in ranks]
        ax.boxplot(data,labels=ranks,showfliers=False)
        ax.axhline(0,color="black",lw=.8)
        ax.set_ylabel("New CI - Old CI (0–1 scale)")
        ax.set_title("Per-position difference: new block MSA vs old chrM CI")
        ax.tick_params(axis="x",rotation=45)
        fig.tight_layout();fig.savefig(out/"score_delta_by_rank.png",dpi=180);plt.close(fig)
        # Old/new gene by rank summaries using distinct positions per symbol
        if "symbol" in per_position:
            tab=per_position.copy()
            mats={}
            for kind in ["old","new"]:
                columns=[f"{r}_rot_ref" if kind=="old" else f"{r}_new_ci" for r in ranks]
                m=tab.groupby("symbol",dropna=True)[columns].mean()
                m.columns=ranks
                mats[kind]=m
            for kind,m in mats.items():
                fig,ax=plt.subplots(figsize=(max(9,len(ranks)),max(6,len(m)*.25)))
                im=ax.imshow(m.values,aspect="auto",vmin=0,vmax=1,cmap="viridis")
                ax.set_xticks(range(len(ranks)),ranks,rotation=45,ha="right")
                ax.set_yticks(range(len(m)),m.index,fontsize=7)
                ax.set_title(f"Mean {kind.upper()} CI by gene and rank")
                fig.colorbar(im,ax=ax,label="CI (0–1)")
                fig.tight_layout();fig.savefig(out/f"gene_heatmap_{kind}.png",dpi=180);plt.close(fig)
    if do_roc and ranks:
        y=labels(variant)
        colors=plt.get_cmap("tab10")
        fig,axes=plt.subplots(1,2,figsize=(13,6),sharey=True)
        auc_rows=[]
        for i,kind in enumerate(["old","new"]):
            ax=axes[i]
            for j,rank in enumerate(ranks):
                col=f"{rank}_rot_ref" if kind=="old" else f"{rank}_new_ci"
                score=pd.to_numeric(variant[col],errors="coerce")
                ok=y.notna()&score.notna()
                if y[ok].nunique()<2:continue
                # Preserve original plotting-script ROC orientation (raw CI).
                auc=roc_auc_score(y[ok].astype(int),score[ok])
                fpr,tpr,_=roc_curve(y[ok].astype(int),score[ok])
                ax.plot(fpr,tpr,color=colors(j%10),label=f"{rank}: {auc:.3f}")
                auc_rows.append({"version":kind,"rank":rank,"auc_raw_ci":auc,"n_variants":int(ok.sum())})
            ax.plot([0,1],[0,1],"k--",alpha=.5)
            ax.set_title(kind.upper()+" CI ROC (raw score orientation)")
            ax.set_xlabel("False positive rate")
            ax.set_ylabel("True positive rate")
            ax.legend(fontsize=7)
        fig.tight_layout();fig.savefig(out/"paired_roc.png",dpi=180);plt.close(fig)
        pd.DataFrame(auc_rows).to_csv(out/"roc_auc_summary.csv",index=False)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--old-csv",required=True)
    p.add_argument("--results-dir",default="block_results_v6")
    p.add_argument("--new-metric",choices=METRICS,default="CI_MSA")
    p.add_argument("--out-dir",default="block_results_v6/ci_comparison")
    p.add_argument("--strict-all-ranks",action="store_true")
    p.add_argument("--no-plots",action="store_true")
    args=p.parse_args()
    out=Path(args.out_dir);out.mkdir(parents=True,exist_ok=True)
    variant,positions,summary,ranks=load_comparison(args.old_csv,args.results_dir,args.new_metric,args.strict_all_ranks)
    # Keep only relevant columns in the exported variant table.
    baseline=[c for c in CORE if c in variant]
    rank_cols=[]
    for rank in ranks:
        rank_cols += [f"{rank}_rot_ref",f"{rank}_new_ci",f"{rank}_new_ci_msa",
                      f"{rank}_new_ci_callable",f"{rank}_n_callable",f"{rank}_n_species",f"{rank}_n_gap"]
    variant[baseline+rank_cols].to_csv(out/"variant_old_new_ci.csv",index=False)
    positions[["POS"]+rank_cols].to_csv(out/"position_old_new_ci.csv",index=False)
    summary.to_csv(out/"comparison_summary.csv",index=False)
    if not args.no_plots:
        write_plots(variant,positions,ranks,out)
    (out/"comparison_metadata.json").write_text(json.dumps({
        "old_csv":args.old_csv,"new_ci_metric":args.new_metric,
        "new_ci_divided_by_100":True,
        "loaded_ranks":ranks,
        "pending_ranks":[r for r in RANKS if r not in ranks],
        "note":"Old zeros preserved, may represent missing data. Inspect old_zero_fraction."
    },indent=2)+"\n")
    print(f"Comparison saved to {out}. Loaded {len(ranks)}/10 ranks.")

if __name__=="__main__":main()
