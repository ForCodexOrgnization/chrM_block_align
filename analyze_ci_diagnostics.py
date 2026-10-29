#!/usr/bin/env python3
"""Diagnostics for old vs new chrM CI (0-100): boundaries, coverage, paired ROC.

Inputs are the output tables from compare_old_new_ci.py and canonical human
GenBank (data/sequence.gb). This does not rerun MAFFT or re-merge rank FASTAs.

Paired bootstrap samples human positions (clusters) with replacement. All
variants/labels at a sampled position receive the same multiplicity. AUC
differences always compare the same eligible variants for each rank.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

RANKS=("kingdom","phylum","subphylum","superclass","class",
       "superorder","order","suborder","infraorder","parvorder")
POSITIVE={"Pathogenic","Pathogenic/Likely pathogenic","Likely pathogenic"}

def boundaries_from_genbank(path, length=16569):
    from Bio import SeqIO
    record=SeqIO.read(path,"genbank")
    if len(record.seq)!=length:
        raise ValueError(f"Expected {length} bp human chrM GenBank")
    boundaries=set()
    for feature in record.features:
        if feature.type in {"CDS","rRNA","tRNA"} and feature.qualifiers.get("gene"):
            boundaries.add(int(feature.location.start)+1)
            boundaries.add(int(feature.location.end)%length+1)
    if not boundaries:raise ValueError("No annotated CDS/rRNA/tRNA boundaries found")
    return np.array(sorted(boundaries),dtype=int)

def circular_distance(pos, boundaries, length=16569):
    p=np.asarray(pos,dtype=int)[:,None]
    b=np.asarray(boundaries,dtype=int)[None,:]
    delta=np.abs(p-b)
    return np.minimum(delta,length-delta).min(axis=1)

def diagnose(positions,ranks,genbank,out,thresholds):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bounds=boundaries_from_genbank(genbank)
    dist=circular_distance(positions.POS.to_numpy(),bounds)
    boundary_bins=[-1,5,10,25,50,100,100000]
    boundary_labels=["0-5","6-10","11-25","26-50","51-100",">100"]
    cov_bins=[-0.001,*thresholds,1.000001]
    cov_labels=[f"{a:.0%}-{b:.0%}" for a,b in zip([0,*thresholds],[*thresholds,1])]
    long=[]
    for rank in ranks:
        old=f"{rank}_rot_ref";new=f"{rank}_new_ci"
        callable_col=f"{rank}_n_callable";total_col=f"{rank}_n_species"
        subset=positions[["POS",old,new,callable_col,total_col]].copy()
        subset["distance_to_boundary_bp"]=dist
        subset["abs_delta_ci"]=(subset[new]-subset[old]).abs()
        subset["delta_ci"]=subset[new]-subset[old]
        denominator=pd.to_numeric(subset[total_col],errors="coerce")
        subset["callable_fraction"]=pd.to_numeric(subset[callable_col],errors="coerce")/denominator.where(denominator>0)
        subset["rank"]=rank
        subset["boundary_bin"]=pd.cut(subset.distance_to_boundary_bp,boundary_bins,labels=boundary_labels)
        subset["coverage_bin"]=pd.cut(subset.callable_fraction,cov_bins,labels=cov_labels)
        long.append(subset[["POS","rank","distance_to_boundary_bp","boundary_bin","coverage_bin",
                             "callable_fraction","delta_ci","abs_delta_ci"]])
    long=pd.concat(long,ignore_index=True)
    long.to_csv(out/"position_diagnostics.csv",index=False)
    outputs={}
    for by,name in [("boundary_bin","boundary_summary"),("coverage_bin","coverage_summary")]:
        summaries=[]
        for (rank,binlabel),g in long.groupby(["rank",by],observed=False):
            vals=g["abs_delta_ci"].dropna()
            summaries.append({"rank":rank,by:str(binlabel),"n_positions":int(len(vals)),
                              "mean_abs_delta_ci":float(vals.mean()) if len(vals) else np.nan,
                              "median_abs_delta_ci":float(vals.median()) if len(vals) else np.nan,
                              "mean_delta_ci":float(g.delta_ci.mean()),
                              "mean_callable_fraction":float(g.callable_fraction.mean())})
        outputs[name]=pd.DataFrame(summaries)
        outputs[name].to_csv(out/f"{name}.csv",index=False)
    for name,category,category_order in [
        ("boundary","boundary_bin",boundary_labels),
        ("coverage","coverage_bin",cov_labels)]:
        fig,ax=plt.subplots(figsize=(12,6))
        for rank in ranks:
            sub=outputs[f"{name}_summary"]
            sub=sub[sub["rank"]==rank].set_index(category).reindex(category_order)
            ax.plot(range(len(category_order)),sub.mean_abs_delta_ci,marker="o",label=rank)
        ax.set_xticks(range(len(category_order)),category_order)
        ax.set_ylabel("Mean |New CI - Old CI| (percentage points)")
        ax.set_xlabel("Distance to nearest annotated gene boundary (bp)" if name=="boundary" else "Callable fraction")
        ax.set_title(f"Old/new CI differences by {name} (distinct human positions)")
        ax.legend(fontsize=8,ncol=2)
        fig.tight_layout();fig.savefig(out/f"{name}_abs_ci_difference.png",dpi=180);plt.close(fig)
    return outputs

def prepare_roc(variants,rank):
    old=f"{rank}_rot_ref";new=f"{rank}_new_ci"
    required={"POS",old,new,"clinvar_interp_new","mitomap_status_new"}
    missing=required-set(variants.columns)
    if missing:raise ValueError(f"Variant table missing {sorted(missing)}")
    clin=variants.clinvar_interp_new.fillna("").astype(str).str.strip()
    mito=variants.mitomap_status_new.fillna("").astype(str)
    positive=clin.isin(POSITIVE)|mito.str.contains("Cfrm",regex=False)
    negative=clin.eq("Benign")&~positive
    mask=(positive|negative)&variants[old].notna()&variants[new].notna()
    selected=variants.loc[mask,["POS",old,new]].copy()
    selected["label"]=positive[mask].astype(np.int8).to_numpy()
    # Exact duplicate variant rows should not artificially add observations.
    if {"REF","ALT"}.issubset(variants.columns):
        selected["REF"]=variants.loc[mask,"REF"].to_numpy()
        selected["ALT"]=variants.loc[mask,"ALT"].to_numpy()
        selected=selected.drop_duplicates(subset=["POS","REF","ALT","label"])
    else:
        selected=selected.drop_duplicates(subset=["POS",old,new,"label"])
    return selected

def bootstrap_auc(variants,ranks,n_boot,seed,out):
    from sklearn.metrics import roc_auc_score
    rng=np.random.default_rng(seed)
    auc_rows=[]
    all_boot=[]
    for rank in ranks:
        d=prepare_roc(variants,rank)
        old=f"{rank}_rot_ref";new=f"{rank}_new_ci"
        if d.empty or d.label.nunique()!=2:
            print(f"[SKIP] ROC {rank}: no valid paired positive/negative set")
            continue
        # Weighted AUC reproduces position-cluster sampling without constructing
        # a potentially enormous repeated dataframe.
        clusters, inverse=np.unique(d.POS.to_numpy(),return_inverse=True)
        y=d.label.to_numpy()
        score_old=d[old].to_numpy(dtype=float)
        score_new=d[new].to_numpy(dtype=float)
        point_old=roc_auc_score(y,score_old)
        point_new=roc_auc_score(y,score_new)
        differences=[]
        for b in range(n_boot):
            sampled=rng.integers(0,len(clusters),size=len(clusters))
            counts=np.bincount(sampled,minlength=len(clusters))
            weights=counts[inverse]
            if not np.any(weights[y==0]) or not np.any(weights[y==1]):
                continue
            old_auc=roc_auc_score(y,score_old,sample_weight=weights)
            new_auc=roc_auc_score(y,score_new,sample_weight=weights)
            differences.append(new_auc-old_auc)
            all_boot.append({"rank":rank,"replicate":b,"auc_old":old_auc,"auc_new":new_auc,
                             "delta_auc":new_auc-old_auc})
        if not differences:
            print(f"[WARN] No valid bootstrap replicates for {rank}")
            continue
        lo,hi=np.quantile(differences,[.025,.975])
        auc_rows.append({"rank":rank,"n_variants":len(d),"n_positions":len(clusters),
                         "auc_old":point_old,"auc_new":point_new,
                         "delta_auc":point_new-point_old,
                         "delta_auc_ci_2_5":lo,"delta_auc_ci_97_5":hi,
                         "bootstrap_replicates":len(differences),
                         "bootstrap_seed":seed})
        print(f"[ROC] {rank}: delta={point_new-point_old:+.4f} 95% CI [{lo:+.4f},{hi:+.4f}]")
    result=pd.DataFrame(auc_rows)
    result.to_csv(out/"paired_roc_bootstrap_summary.csv",index=False)
    pd.DataFrame(all_boot).to_csv(out/"paired_roc_bootstrap_replicates.csv",index=False)
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--comparison-dir",default="block_results_v6/ci_comparison")
    ap.add_argument("--human-genbank",default="data/sequence.gb")
    ap.add_argument("--out-dir",default=None)
    ap.add_argument("--bootstrap-reps",type=int,default=1000)
    ap.add_argument("--seed",type=int,default=20261029)
    ap.add_argument("--coverage-thresholds",nargs=2,type=float,default=[.5,.9],
                    metavar=("LOW","HIGH"))
    args=ap.parse_args()
    if not 0<args.coverage_thresholds[0]<args.coverage_thresholds[1]<1:
        ap.error("--coverage-thresholds must satisfy 0 < LOW < HIGH < 1")
    if args.bootstrap_reps<1:ap.error("--bootstrap-reps must be >= 1")
    base=Path(args.comparison_dir)
    out=Path(args.out_dir) if args.out_dir else base/"diagnostics"
    out.mkdir(parents=True,exist_ok=True)
    positions=pd.read_csv(base/"position_old_new_ci.csv")
    variants=pd.read_csv(base/"variant_old_new_ci.csv",low_memory=False)
    ranks=[r for r in RANKS if {f"{r}_rot_ref",f"{r}_new_ci",f"{r}_n_callable",f"{r}_n_species"}.issubset(positions.columns)]
    if not ranks:raise ValueError("No complete rank columns in comparison outputs")
    print("Ranks:",", ".join(ranks))
    diagnose(positions,ranks,args.human_genbank,out,args.coverage_thresholds)
    bootstrap_auc(variants,ranks,args.bootstrap_reps,args.seed,out)
    print("Outputs:",out)

if __name__=="__main__":main()
