#!/usr/bin/env python3
"""Rank-wise correlation of historical CI with pad100 Original/Mean/Min/Max.

Historical scores come from {rank}_rot_ref in the old annotation CSV.
This is an external consistency check, NOT a ground-truth accuracy measure.
Compare a strictly shared set of unique human positions per rank and group.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from compare_old_new_ci import RANKS

METHODS=("original","mean","min","max")
METRICS=("CI_MSA","CI_callable")


def historical_positions(path, rank, zero_missing):
    use=["POS",f"{rank}_rot_ref"]
    raw=pd.read_csv(path,usecols=use,low_memory=False)
    raw["POS"]=pd.to_numeric(raw.POS,errors="raise").astype(int)
    raw["old"]=pd.to_numeric(raw[use[1]],errors="coerce")
    if not raw.POS.between(1,16569).all():
        raise ValueError(f"{path}: positions out of bounds")
    # Duplicated alternate alleles share one position. Fail on conflicts.
    nuniq=raw.groupby("POS").old.nunique(dropna=True)
    if nuniq.gt(1).any():
        raise ValueError(f"Historical {rank} CI varies across ALT alleles at POS: {nuniq[nuniq>1].index[:8].tolist()}")
    old=raw.groupby("POS",as_index=False).old.first()
    if zero_missing:
        old.loc[old.old==0,"old"]=np.nan
    if not old.old.dropna().between(0,100).all():
        raise ValueError(f"Historical {rank} CI outside 0..100")
    return old


def load_rank(old_csv, root, rank, metric, zero_missing):
    old=historical_positions(old_csv,rank,zero_missing)
    d=old
    human_ref=None
    for method in METHODS:
        path=(root/rank/"conservation"/"conservation_ci.csv"
              if method=="original"
              else root/rank/"overlap_aware_ci"/method/"conservation_ci.csv")
        z=pd.read_csv(path)
        if (len(z)!=16569 or z.Position.tolist()!=list(range(1,16570)) or
                not {"Position","Homo sapiens",metric}.issubset(z)):
            raise ValueError(f"{path}: incorrect CI schema/position coverage")
        bases=z["Homo sapiens"].astype(str).str.upper().to_numpy()
        if human_ref is not None and not np.array_equal(bases,human_ref):
            raise ValueError(f"Human reference mismatch {path}")
        human_ref=bases
        t=pd.DataFrame({"POS":z.Position,method:pd.to_numeric(z[metric],errors="coerce")})
        if not t[method].dropna().between(0,100).all():
            raise ValueError(f"{path}: scores outside 0..100")
        d=d.merge(t,on="POS",how="inner",validate="one_to_one")
    meta=root/rank/"overlap_aware_ci"/"overlap_position_ci.csv"
    q=pd.read_csv(meta,usecols=["Position","is_overlap","genes"])
    if len(q)!=16569 or q.Position.tolist()!=list(range(1,16570)):
        raise ValueError(f"Bad metadata {meta}")
    d=d.merge(q.rename(columns={"Position":"POS"}),on="POS",validate="one_to_one")
    d["is_overlap"]=d.is_overlap.astype(int)
    d["score_changed"]=np.any(np.column_stack([
        ~np.isclose(d[m],d.original,rtol=0,atol=1e-8,equal_nan=True)
        for m in METHODS[1:]]),axis=1)
    # Scoring rules promise identical non-overlap positions.
    outside=d.is_overlap.eq(0)
    for m in METHODS[1:]:
        if not np.all(np.isclose(d.loc[outside,m],d.loc[outside,"original"],
                                 rtol=0,atol=1e-8,equal_nan=True)):
            raise ValueError(f"{rank}: {m} changes non-overlap positions")
    return d


def groups_for(d):
    result={"genome_wide":d,
            "gene_gene_overlap":d[d.is_overlap==1],
            "changed_overlap":d[(d.is_overlap==1)&d.score_changed]}
    for pair,subset in d.loc[d.is_overlap==1].groupby("genes"):
        result["pair:"+pair]=subset
    return result


def compare(subset,rank,metric,group,min_n):
    complete=subset.dropna(subset=["old",*METHODS]).copy()
    result=[]
    for method in METHODS:
        row=dict(rank=rank,metric=metric,group=group,method=method,
                 n_positions=int(len(complete)),
                 historical_zero_excluded=None,
                 pearson_r=np.nan,spearman_rho=np.nan,mae=np.nan,
                 rmse=np.nan,mean_signed_difference=np.nan,status="insufficient_positions")
        if len(complete)>=1:
            delta=complete[method].to_numpy(dtype=float)-complete.old.to_numpy(dtype=float)
            row.update(mae=float(np.abs(delta).mean()),
                       rmse=float(np.sqrt((delta**2).mean())),
                       mean_signed_difference=float(delta.mean()))
        if len(complete)<min_n:
            result.append(row);continue
        x=complete.old.to_numpy(dtype=float)
        y=complete[method].to_numpy(dtype=float)
        if np.std(x)==0 or np.std(y)==0:
            row["status"]="constant_score"
        else:
            row.update(pearson_r=float(pearsonr(x,y).statistic),
                       spearman_rho=float(spearmanr(x,y).statistic),
                       status="ok")
        result.append(row)
    return result,complete


def plot_overlap(d,rank,metric,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data=d[d.is_overlap.eq(1)].dropna(subset=["old",*METHODS])
    if data.empty:return
    fig,axes=plt.subplots(2,2,figsize=(10,9),sharex=True,sharey=True)
    for method,ax in zip(METHODS,axes.flat):
        ax.scatter(data.old,data[method],s=15,alpha=.7)
        ax.plot([0,100],[0,100],"k--",alpha=.4)
        ax.set(xlim=(0,100),ylim=(0,100),title=method,
               xlabel="Historical CI",ylabel="New CI")
        ax.grid(alpha=.15)
    fig.suptitle(f"{rank} {metric}: human gene-gene overlap only (n={len(data)})")
    fig.tight_layout()
    fig.savefig(out/f"{rank}_overlap_scatter.png",dpi=170)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--old-csv",default="data/2026_10_10_mito_cons_anns_main.csv")
    p.add_argument("--results-root",default="block_results_genbank_trnascan_mitos_fallback_v2/pad100")
    p.add_argument("--metrics",nargs="+",choices=METRICS,default=list(METRICS))
    p.add_argument("--ranks",nargs="+",choices=RANKS,default=list(RANKS))
    p.add_argument("--old-zero-is-missing",action="store_true",
                   help="Exclude historical zero values rather than treating as zero conservation")
    p.add_argument("--min-n",type=int,default=3)
    p.add_argument("--out-dir")
    a=p.parse_args()
    root=Path(a.results_root)
    out=Path(a.out_dir or root/"overlap_old_ci_correlation")
    out.mkdir(parents=True,exist_ok=True)
    rows=[]
    cohorts=[]
    for metric in a.metrics:
        metric_dir=out/metric
        metric_dir.mkdir(exist_ok=True)
        for rank in a.ranks:
            d=load_rank(a.old_csv,root,rank,metric,a.old_zero_is_missing)
            d["rank"]=rank;d["metric"]=metric
            cohorts.append(d)
            for group,sub in groups_for(d).items():
                values,complete=compare(sub,rank,metric,group,a.min_n)
                rows.extend(values)
            plot_overlap(d,rank,metric,metric_dir)
            print(metric,rank,"overlap n=",int(d.is_overlap.sum()),
                  "overlap old nonnull=",int(d[d.is_overlap==1].old.notna().sum()))
    summary=pd.DataFrame(rows)
    summary.to_csv(out/"correlation_summary.csv",index=False)
    pd.concat(cohorts,ignore_index=True).to_csv(out/"paired_position_scores.csv",index=False)
    (out/"metadata.json").write_text(json.dumps({
        "historical_reference":a.old_csv,
        "old_columns":"{rank}_rot_ref",
        "note":"Historical CI is a comparator, not truth; provenance and denominator need independent confirmation.",
        "one_row_per_human_position":True,
        "pairing":"shared complete cases across all four new methods per rank/group",
        "old_zero_is_missing":a.old_zero_is_missing,
        "min_n":a.min_n,
    },indent=2)+"\n")
    print("\nGene-gene overlap comparison:")
    display=summary[summary.group=="gene_gene_overlap"]
    print(display[["rank","metric","method","n_positions","pearson_r",
                   "spearman_rho","mae","status"]].to_string(index=False))
    print("Outputs:",out)


if __name__=="__main__":
    main()
