#!/usr/bin/env python3
"""Three-way chrM ROC: old rot_ref, gene-only, padded, with clustered bootstrap.

All scores are 0-100. Uses identical eligible labeled variants per analysis.
Within-gene AUC counts only pathogenic/benign pairs in the SAME human gene,
weighted by the number of such pairs (not pooled across genes).
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.metrics import roc_auc_score, roc_curve

POSITIVE={"Pathogenic","Pathogenic/Likely pathogenic","Likely pathogenic"}
METHODS=("old","gene_only","padded")
PAIRS=(("gene_only","old"),("padded","old"),("padded","gene_only"))
WINDOWS=((0,5),(6,25),(26,100),(101,16569))

def annotation(path):
    rec=SeqIO.read(path,"genbank")
    L=len(rec.seq)
    boundaries=set();lookup={i:[] for i in range(1,L+1)}
    for f in rec.features:
        if f.type not in ("CDS","rRNA","tRNA"):continue
        name=f.qualifiers.get("gene",[])
        if not name:continue
        a,b=int(f.location.start),int(f.location.end)
        boundaries.update(((a%L)+1,(b%L)+1))
        for i in range(a,b):
            lookup[i%L+1].append((name[0],f.type))
    # Prefer the shortest specific annotation for overlapping canonical genes.
    gene=[];category=[]
    precedence={"tRNA":0,"rRNA":1,"CDS":2}
    for p in range(1,L+1):
        fs=sorted(lookup[p],key=lambda q:precedence[q[1]])
        gene.append(fs[0][0] if fs else "intergenic")
        category.append(fs[0][1] if fs else "intergenic")
    b=np.array(sorted(boundaries))
    pos=np.arange(1,L+1)[:,None]
    d=np.abs(pos-b[None,:])
    return pd.DataFrame({"POS":np.arange(1,L+1),"gene":gene,"category":category,
                         "boundary_distance":np.minimum(d,L-d).min(axis=1)})

def load(old_csv,gene_only,padded,rank,metric,genbank):
    source=pd.read_csv(old_csv,low_memory=False)
    old=f"{rank}_rot_ref"
    required={"POS","REF","ALT","clinvar_interp_new","mitomap_status_new",old}
    if not required.issubset(source):raise ValueError(f"Old CSV lacks: {required-set(source)}")
    source.POS=pd.to_numeric(source.POS,errors="raise").astype(int)
    source["old"]=pd.to_numeric(source[old],errors="coerce")
    for label,path in (("gene_only",gene_only),("padded",padded)):
        d=pd.read_csv(path,usecols=["Position",metric])
        if len(d)!=16569 or d.Position.nunique()!=16569 or not d.Position.between(1,16569).all():
            raise ValueError(f"{path} is not a complete 16569-position file")
        d=d.rename(columns={"Position":"POS",metric:label})
        source=source.merge(d,on="POS",how="left",validate="many_to_one")
    if not source[list(METHODS)].stack().dropna().between(0,100).all():
        raise ValueError("Expected CI values in range 0-100")
    clin=source.clinvar_interp_new.fillna("").astype(str).str.strip()
    mito=source.mitomap_status_new.fillna("").astype(str)
    positive=clin.isin(POSITIVE)|mito.str.contains("Cfrm",regex=False)
    negative=clin.eq("Benign")&~positive
    data=source.loc[(positive|negative)&source[list(METHODS)].notna().all(axis=1),
                    ["POS","REF","ALT",*METHODS]].copy()
    data["label"]=positive.loc[data.index].astype(int).to_numpy()
    data=data.drop_duplicates(subset=["POS","REF","ALT","label"])
    data=data.merge(annotation(genbank),on="POS",how="left",validate="many_to_one")
    return data

def auc(y,s,w=None):
    if w is None:w=np.ones(len(y))
    if w[y==0].sum()==0 or w[y==1].sum()==0:return np.nan
    return float(roc_auc_score(y,s,sample_weight=w))

def within_gene_auc(d,w=None):
    # Comparable pathogenic/benign pairs must be within the same gene.
    if w is None:w=np.ones(len(d))
    numerator=0.;pair_total=0.
    for _,idx in d.groupby("gene",sort=False).indices.items():
        y=d.label.to_numpy()[idx]
        a=d["_score"].to_numpy()[idx]
        wt=w[idx]
        p=wt[y==1].sum();n=wt[y==0].sum()
        if p==0 or n==0:continue
        pair=p*n
        value=auc(y,a,wt)
        numerator+=value*pair;pair_total+=pair
    return numerator/pair_total if pair_total else np.nan

def evaluate(d,group,replicates,rng,min_class,out):
    npos=int((d.label==1).sum());nneg=int((d.label==0).sum())
    counts={"group":group,"n_variants":len(d),"n_positions":d.POS.nunique(),
            "n_pathogenic":npos,"n_benign":nneg}
    if min(npos,nneg)<min_class:
        return [{**counts,"comparison":f"{b}-minus-{a}","status":"insufficient_labels"}
                for b,a in PAIRS]
    y=d.label.to_numpy()
    scores={m:d[m].to_numpy(dtype=float) for m in METHODS}
    genes=d.gene.to_numpy()
    positions,inv=np.unique(d.POS.to_numpy(),return_inverse=True)
    use_within=group=="within_gene"
    def calc(w):
        if use_within:
            local=d[["gene","label"]].copy()
            return {m:within_gene_auc(local.assign(_score=scores[m]),w) for m in METHODS}
        return {m:auc(y,scores[m],w) for m in METHODS}
    point=calc(np.ones(len(d)))
    draw={pair:[] for pair in PAIRS}
    for _ in range(replicates):
        sample=rng.integers(0,len(positions),size=len(positions))
        w=np.bincount(sample,minlength=len(positions))[inv]
        values=calc(w)
        for b,a in PAIRS:
            if np.isfinite(values[a]) and np.isfinite(values[b]):
                draw[(b,a)].append(values[b]-values[a])
    rows=[]
    for b,a in PAIRS:
        ds=np.asarray(draw[(b,a)])
        row={**counts,"comparison":f"{b}-minus-{a}","status":"ok" if len(ds)>0 else "no_valid_replicates",
             "auc_reference":point[a],"auc_candidate":point[b],
             "delta_auc":point[b]-point[a],"valid_replicates":len(ds),
             "ci_lower":float(np.quantile(ds,.025)) if len(ds) else np.nan,
             "ci_upper":float(np.quantile(ds,.975)) if len(ds) else np.nan}
        rows.append(row)
    return rows

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--old-csv",default="data/2026_10_10_mito_cons_anns_main.csv")
    p.add_argument("--gene-only-csv",default="block_results_v6/class/conservation/conservation_ci.csv")
    p.add_argument("--padded-csv",default="block_results_padded_v1/pad50/class/conservation/conservation_ci.csv")
    p.add_argument("--human-genbank",default="data/sequence.gb")
    p.add_argument("--rank",default="class")
    p.add_argument("--metric",choices=["CI_MSA","CI_callable"],default="CI_MSA")
    p.add_argument("--out-dir",default="block_results_padded_v1/pad50/class/roc_evaluation")
    p.add_argument("--bootstrap-reps",type=int,default=500)
    p.add_argument("--seed",type=int,default=20261029)
    p.add_argument("--min-per-class",type=int,default=5)
    a=p.parse_args()
    if a.bootstrap_reps<1:p.error("bootstrap-reps must be positive")
    out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
    d=load(a.old_csv,a.gene_only_csv,a.padded_csv,a.rank,a.metric,a.human_genbank)
    d.to_csv(out/"labeled_variant_scores.csv",index=False)
    groups={"genome_wide":d,"within_gene":d[d.gene!="intergenic"]}
    for category in ("CDS","tRNA","rRNA","intergenic"):
        groups["category_"+category]=d[d.category==category].copy()
    for lower,upper in WINDOWS:
        groups[f"boundary_{lower}_{upper}bp"]=d[d.boundary_distance.between(lower,upper)].copy()
    # Distinct variants for all comparisons, same positions within each group.
    rng=np.random.default_rng(a.seed)
    rows=[]
    for group,subset in groups.items():
        rows.extend(evaluate(subset,group,a.bootstrap_reps,rng,a.min_per_class,out))
    result=pd.DataFrame(rows)
    result.to_csv(out/"three_way_auc_summary.csv",index=False)
    print(result[["group","comparison","n_pathogenic","n_benign","delta_auc","ci_lower","ci_upper","status"]].to_string(index=False))
    # Standard genome-wide ROC plot uses the shared set of variants.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,6))
    for method in METHODS:
        if d.label.nunique()<2:continue
        fpr,tpr,_=roc_curve(d.label,d[method])
        ax.plot(fpr,tpr,label=f"{method}: {auc(d.label.to_numpy(),d[method].to_numpy()):.3f}")
    ax.plot([0,1],[0,1],"k--",alpha=.5)
    ax.set(xlabel="False positive rate",ylabel="True positive rate",title="chrM CI three-way ROC")
    ax.legend();fig.tight_layout();fig.savefig(out/"three_way_genome_roc.png",dpi=180);plt.close(fig)

if __name__=="__main__":main()
