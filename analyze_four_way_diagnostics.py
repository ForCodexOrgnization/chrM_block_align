#!/usr/bin/env python3
"""Four-method, ten-rank diagnostics: AUC, gene type, boundaries, within gene,
core/flank CI shifts, and padded flank eligibility. Existing outputs read-only.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve
from Bio import SeqIO
from compare_all_ranks_padding_roc import load, METHODS, COMPARISONS, RANKS

def human_regions(gb):
    record=SeqIO.read(gb,"genbank");L=len(record.seq)
    if L!=16569:raise ValueError("Expected human chrM length 16569")
    gene=np.full(L,"intergenic",dtype=object)
    kind=np.full(L,"intergenic",dtype=object)
    priority=np.full(L,99,dtype=int)
    edges=[]
    # Human features only. A nucleotide with multiple features is assigned by
    # tRNA > rRNA > CDS priority, with reproducible feature-order tie breaks.
    levels={"tRNA":0,"rRNA":1,"CDS":2}
    for feature in record.features:
        if feature.type not in levels or not feature.qualifiers.get("gene"):continue
        g=feature.qualifiers["gene"][0]
        for part in feature.location.parts:
            start,end=int(part.start),int(part.end)
            edges.extend([(start%L)+1,(end%L)+1])
            for i in range(start,end):
                j=i%L
                if levels[feature.type]<priority[j]:
                    gene[j]=g;kind[j]=feature.type;priority[j]=levels[feature.type]
    b=np.asarray(sorted(set(edges)),dtype=int)
    pos=np.arange(1,L+1)[:,None];delta=np.abs(pos-b[None,:])
    distance=np.minimum(delta,L-delta).min(axis=1)
    band=np.select([distance<=5,distance<=25,distance<=100],
                   ["0-5","6-25","26-100"],default=">100")
    return pd.DataFrame({"POS":np.arange(1,L+1),"gene":gene,"gene_type":kind,
                         "boundary_bp":distance,"boundary_band":band})

def scores_auc(y,s,w):
    return float(roc_auc_score(y,s,sample_weight=w)) if w[y==0].sum()>0 and w[y==1].sum()>0 else np.nan

def within_gene(y,s,w,genes):
    # Only positive-negative pairs from the SAME gene enter this concordance.
    total=0.0;weighted=0.0
    for gene in np.unique(genes):
        if gene=="intergenic":continue
        idx=np.flatnonzero(genes==gene)
        plus=w[idx][y[idx]==1].sum();minus=w[idx][y[idx]==0].sum()
        if plus==0 or minus==0:continue
        pair_weight=plus*minus
        weighted+=pair_weight*scores_auc(y[idx],s[idx],w[idx])
        total+=pair_weight
    return weighted/total if total else np.nan

def diagnostic_roc(d,out,reps,seed,min_class):
    rng=np.random.default_rng(seed)
    records=[];boot=[]
    groupings={"genome_wide":np.ones(len(d),dtype=bool),
               "within_gene":d.gene.ne("intergenic").to_numpy()}
    for kind in ("CDS","tRNA","rRNA","intergenic"):
        groupings["type_"+kind]=d.gene_type.eq(kind).to_numpy()
    for band in ("0-5","6-25","26-100",">100"):
        groupings["boundary_"+band]=d.boundary_band.eq(band).to_numpy()
    for rank in RANKS:
        for group,mask in groupings.items():
            sub=d.loc[mask]
            y=sub.label.to_numpy(dtype=int)
            npos=int(y.sum());nneg=int(len(y)-npos)
            if min(npos,nneg)<min_class:
                records.append({"rank":rank,"group":group,"status":"insufficient_labels",
                                "n_variants":len(sub),"n_pathogenic":npos,"n_benign":nneg})
                continue
            scores={m:sub[f"{rank}_{m}"].to_numpy(dtype=float) for m in METHODS}
            genes=sub.gene.to_numpy()
            cluster,inv=np.unique(sub.POS.to_numpy(),return_inverse=True)
            def measure(w):
                if group=="within_gene":
                    return {m:within_gene(y,scores[m],w,genes) for m in METHODS}
                return {m:scores_auc(y,scores[m],w) for m in METHODS}
            point=measure(np.ones(len(sub)))
            draws={m:[] for m in METHODS}
            diff={f"{a}_minus_{b}":[] for a,b in COMPARISONS}
            for iteration in range(reps):
                chosen=rng.integers(len(cluster),size=len(cluster))
                w=np.bincount(chosen,minlength=len(cluster))[inv]
                v=measure(w)
                if not all(np.isfinite(x) for x in v.values()):continue
                for m in METHODS:draws[m].append(v[m])
                for a,b in COMPARISONS:diff[f"{a}_minus_{b}"].append(v[a]-v[b])
                boot.append({"rank":rank,"group":group,"replicate":iteration,
                             **{f"auc_{m}":v[m] for m in METHODS}})
            row={"rank":rank,"group":group,"status":"ok" if len(draws["original"]) else "no_valid_bootstraps",
                 "n_variants":len(sub),"n_positions":len(cluster),
                 "n_pathogenic":npos,"n_benign":nneg,
                 "valid_bootstraps":len(draws["original"])}
            for m in METHODS:
                vals=np.asarray(draws[m])
                row["auc_"+m]=point[m]
                row["auc_"+m+"_lo"]=np.quantile(vals,.025) if len(vals) else np.nan
                row["auc_"+m+"_hi"]=np.quantile(vals,.975) if len(vals) else np.nan
            for label,values in diff.items():
                a,b=label.split("_minus_")
                vals=np.asarray(values)
                row["delta_"+label]=point[a]-point[b]
                row["delta_"+label+"_lo"]=np.quantile(vals,.025) if len(vals) else np.nan
                row["delta_"+label+"_hi"]=np.quantile(vals,.975) if len(vals) else np.nan
            records.append(row)
    summary=pd.DataFrame(records)
    summary.to_csv(out/"stratified_four_way_auc.csv",index=False)
    pd.DataFrame(boot).to_csv(out/"stratified_four_way_bootstrap.csv",index=False)
    return summary

def draw_auc(summary,out):
    import matplotlib.pyplot as plt
    for group_prefix,filename,labels in (
        ("genome_wide","genome_auc_by_rank.png",["genome_wide"]),
        ("type_","gene_type_auc_by_rank.png",["type_CDS","type_tRNA","type_rRNA"]),
        ("boundary_","boundary_auc_by_rank.png",["boundary_0-5","boundary_6-25","boundary_26-100","boundary_>100"]),
        ("within_gene","within_gene_auc_by_rank.png",["within_gene"])):
        fig,axes=plt.subplots(len(labels),1,figsize=(13,3.3*len(labels)),sharex=True,squeeze=False)
        for ax,g in zip(axes.flat,labels):
            subset=summary[(summary.group==g)&(summary.status=="ok")].set_index("rank").reindex(RANKS)
            for method in METHODS:
                series=subset["auc_"+method]
                ax.plot(range(len(RANKS)),series,marker="o",label=method)
            ax.set_title(g);ax.set_ylabel("ROC AUC");ax.set_ylim(.45,1.05);ax.grid(alpha=.2)
            ax.legend(ncol=4,fontsize=8)
        axes[-1,0].set_xticks(range(len(RANKS)),RANKS,rotation=30,ha="right")
        fig.tight_layout();fig.savefig(out/filename,dpi=180);plt.close(fig)
    # Paired effect plot on shared labels with 95% CI.
    selected=summary[(summary.group=="genome_wide")&(summary.status=="ok")].set_index("rank").reindex(RANKS)
    fig,ax=plt.subplots(figsize=(13,5))
    comparisons=[("gene_only","original"),("pad50","gene_only"),("pad100","gene_only"),("pad100","pad50")]
    x=np.arange(len(RANKS))
    for j,(a,b) in enumerate(comparisons):
        tag=f"delta_{a}_minus_{b}"
        v=selected[tag].to_numpy(dtype=float)
        lo=selected[tag+"_lo"].to_numpy(dtype=float);hi=selected[tag+"_hi"].to_numpy(dtype=float)
        xx=x+(j-1.5)*.15
        ax.errorbar(xx,v,yerr=np.vstack([np.maximum(0,v-lo),np.maximum(0,hi-v)]),
                    fmt="o",capsize=2,label=f"{a} - {b}")
    ax.axhline(0,color="gray",ls="--")
    ax.set_xticks(x,RANKS,rotation=30,ha="right")
    ax.set_ylabel("Paired delta AUC (95% position bootstrap CI)")
    ax.legend(ncol=2);fig.tight_layout();fig.savefig(out/"genome_delta_auc_ci.png",dpi=180);plt.close(fig)

def position_differences(roots,out,regions):
    # Position-level changes independent of clinical labeling.
    rows=[]
    positions=regions[["POS","gene_type","boundary_band","gene","boundary_bp"]]
    for rank in RANKS:
        data=positions.copy()
        sources={
            "gene_only":Path(roots.gene_only_dir)/rank/"conservation"/"conservation_ci.csv",
            "pad50":Path(roots.padded_dir)/"pad50"/rank/"conservation"/"conservation_ci.csv",
            "pad100":Path(roots.padded_dir)/"pad100"/rank/"conservation"/"conservation_ci.csv"}
        for method,path in sources.items():
            src=pd.read_csv(path,usecols=["Position",roots.metric]).rename(columns={"Position":"POS",roots.metric:method})
            data=data.merge(src,on="POS",validate="one_to_one")
        # Core vs flank based on human annotated positions, not the owner status.
        data["region"]=np.where(data.gene_type=="intergenic","intergenic",
                         np.where(data.boundary_bp<=5,"gene_boundary_0_5","gene_core_or_near_core"))
        for a,b in (("pad50","gene_only"),("pad100","gene_only"),("pad100","pad50")):
            data["change"]=data[a]-data[b]
            for region,g in data.groupby("region"):
                rows.append({"rank":rank,"comparison":a+"_minus_"+b,"region":region,
                             "n_positions":len(g),"mean_delta_ci":g.change.mean(),
                             "mean_abs_delta_ci":g.change.abs().mean(),
                             "frac_abs_delta_gt10":(g.change.abs()>10).mean()})
    summary=pd.DataFrame(rows)
    summary.to_csv(out/"position_ci_change_by_region.csv",index=False)
    import matplotlib.pyplot as plt
    for comparison in summary.comparison.unique():
        subset=summary[summary.comparison==comparison]
        fig,ax=plt.subplots(figsize=(12,5))
        for region,g in subset.groupby("region"):
            vals=g.set_index("rank").reindex(RANKS)
            ax.plot(RANKS,vals.mean_abs_delta_ci,marker="o",label=region)
        ax.set_ylabel("Mean absolute CI difference (percentage points)")
        ax.set_title(comparison);ax.legend();ax.tick_params(axis="x",rotation=30)
        fig.tight_layout();fig.savefig(out/("ci_shift_"+comparison+".png"),dpi=180);plt.close(fig)

def eligibility(padded,out):
    # Entries are one row per gene/species, not per nucleotide.
    frames=[]
    for rank in RANKS:
        for pad in ("pad50","pad100"):
            path=Path(padded)/pad/rank/"flank_context_qc.tsv"
            if not path.exists():raise FileNotFoundError(path)
            d=pd.read_csv(path,sep="\t")
            if not {"gene","species","left_supported","right_supported"}.issubset(d.columns):
                raise ValueError(f"Bad flank QC columns: {path}")
            d["rank"]=rank;d["padding"]=pad
            frames.append(d)
    qc=pd.concat(frames,ignore_index=True)
    summary=qc.groupby(["rank","padding"],as_index=False).agg(
        n_gene_species=("species","size"),
        n_species=("species","nunique"),
        eligible_left=("left_supported","mean"),
        eligible_right=("right_supported","mean"))
    qc["both"]=qc.left_supported.astype(bool)&qc.right_supported.astype(bool)
    both=qc.groupby(["rank","padding"],as_index=False).both.mean().rename(columns={"both":"eligible_both"})
    summary=summary.merge(both,on=["rank","padding"])
    summary.to_csv(out/"flank_eligibility_by_rank.csv",index=False)
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(12,5))
    for pad,g in summary.groupby("padding"):
        ordered=g.set_index("rank").reindex(RANKS)
        ax.plot(RANKS,100*ordered.eligible_both,marker="o",label=pad)
    ax.set_ylabel("% annotated gene/species pairs with both neighboring genes supported")
    ax.tick_params(axis="x",rotation=30);ax.legend()
    fig.tight_layout();fig.savefig(out/"flank_eligibility_by_rank.png",dpi=180);plt.close(fig)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--old-csv",default="data/2026_10_10_mito_cons_anns_main.csv")
    p.add_argument("--gene-only-dir",default="block_results_v6")
    p.add_argument("--padded-dir",default="block_results_padded_v1")
    p.add_argument("--human-genbank",default="data/sequence.gb")
    p.add_argument("--metric",choices=("CI_MSA","CI_callable"),default="CI_MSA")
    p.add_argument("--bootstrap-reps",type=int,default=1000)
    p.add_argument("--seed",type=int,default=20261030)
    p.add_argument("--min-per-class",type=int,default=5)
    p.add_argument("--out-dir",default="block_results_padded_v1/four_way_diagnostics")
    a=p.parse_args()
    if a.bootstrap_reps<1:p.error("bootstrap-reps must be >=1")
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    d=load(a.old_csv,a.padded_dir,a.gene_only_dir,a.metric)
    regions=human_regions(a.human_genbank)
    d=d.merge(regions,on="POS",how="left",validate="many_to_one")
    d.to_csv(out/"shared_labeled_variants_with_regions.csv",index=False)
    summary=diagnostic_roc(d,out,a.bootstrap_reps,a.seed,a.min_per_class)
    draw_auc(summary,out)
    position_differences(a,out,regions)
    eligibility(a.padded_dir,out)
    print("Saved four-way diagnostics:",out)

if __name__=="__main__":main()
