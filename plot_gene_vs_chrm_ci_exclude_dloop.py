#!/usr/bin/env python3
"""Reproduce OLD vs NEW pad100 gene-vs-chrM CI scatterplots excluding D-loop.

20 panels: OLD/New pairs across ten ranks, one point per human reference base.
OLD X: <rank>_rot_ref from the historical annotation CSV.
NEW X: CI_MSA (or CI_callable) in the selected pad100 conservation_ci.csv.
Both Y: <rank>_ref_gene from the historical annotation CSV.

Before pairing, remove all reference positions in the human GenBank D-loop.
Missing scores are excluded pairwise using the SAME old/new intersection.
The default D-loop coordinates for NC_012920.1 are positions 16024..16569
and 1..576 (1-based inclusive), i.e., 1122 bp; the GenBank feature is
also inspected so a changed reference does not silently pass.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from plot_gene_vs_genbank_trnascan_pad100 import (
    RANKS, COLOR, old_scores, new_scores, categories_from_genbank, scatter
)

REFERENCE_LENGTH=16569
DLOOP_POSITIONS=set(range(16024,16570)) | set(range(1,577))

def get_dloop_positions(human_genbank):
    from Bio import SeqIO
    rec=SeqIO.read(human_genbank,"genbank")
    if len(rec.seq)!=REFERENCE_LENGTH:
        raise ValueError("Expected 16569-bp human mitochondrial reference")
    detected=set()
    for feature in rec.features:
        fields=[feature.type]
        for k in ("gene","note","standard_name","product","region_name"):
            fields.extend(map(str,feature.qualifiers.get(k,[])))
        label=" ".join(fields).lower().replace("_","").replace("-","")
        if "dloop" in label or "control region" in label:
            detected.update(int(i)+1 for i in feature.location)
    if detected and detected!=DLOOP_POSITIONS:
        raise ValueError(
            "Human GenBank D-loop differs from NC_012920.1 expected interval "
            f"(detected {len(detected)} bp, expected 1122). Review reference."
        )
    return DLOOP_POSITIONS

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--annotation-csv",
                   default="data/2026_10_10_mito_cons_anns_main.csv")
    p.add_argument("--human-genbank",default="data/sequence.gb")
    p.add_argument("--genbank-root",
                   default="block_results_genbank_trnascan_mitos_fallback_v2")
    p.add_argument("--metric",choices=("CI_MSA","CI_callable"),default="CI_MSA")
    p.add_argument("--out-dir",
                   default="gene_vs_chrm_ci_old_vs_new_pad100_no_dloop")
    a=p.parse_args()
    out=Path(a.out_dir)
    out.mkdir(parents=True,exist_ok=True)
    excluded=get_dloop_positions(a.human_genbank)
    old=old_scores(a.annotation_csv)
    categories=categories_from_genbank(a.human_genbank)
    retain=~categories.index.isin(excluded)
    fig,axes=plt.subplots(5,4,figsize=(21,17),sharex=True,sharey=True)
    rows=[]
    for i,rank in enumerate(RANKS):
        path=Path(a.genbank_root)/"pad100"/rank/"conservation"/"conservation_ci.csv"
        if not path.is_file():
            raise FileNotFoundError(path)
        new=new_scores(path,a.metric).reindex(categories.index)
        y=old[f"{rank}_ref_gene"].reindex(categories.index)
        x_old=old[f"{rank}_rot_ref"].reindex(categories.index)
        # Identical positions are used in both panels of each pair; remove
        # D-loop BEFORE computing correlations, not just from the rendered dots.
        paired=(retain & np.isfinite(y.to_numpy())
                & np.isfinite(x_old.to_numpy()) & np.isfinite(new.to_numpy()))
        colors=categories.to_numpy()[paired]
        yy=y.to_numpy()[paired]
        ax_old,ax_new=axes.flat[2*i:2*i+2]
        for method,ax,xx in (
            ("old",ax_old,x_old.to_numpy()[paired]),
            ("v2_pad100",ax_new,new.to_numpy()[paired]),
        ):
            n,pearson,spearman=scatter(ax,xx,yy,colors,
                                       f"{rank} ({'old' if method=='old' else 'new pad100'})")
            rows.append(dict(rank=rank,method=method,metric=a.metric,
                             n_reference_positions=n,n_dloop_excluded=len(excluded),
                             pearson_r=pearson,spearman_rho=spearman))
    for ax in axes.flat:
        ax.set_xlabel("chrM CI",fontsize=8)
        ax.set_ylabel("gene CI",fontsize=8)
        ax.tick_params(labelsize=7)
    legend=[Line2D([0],[0],marker="o",linestyle="",markersize=6,
                   color=COLOR[k],label=k) for k in COLOR]
    fig.legend(handles=legend,loc="upper center",ncol=4,frameon=False,
               bbox_to_anchor=(.5,.967),fontsize=9)
    fig.suptitle("Gene vs chrM CI scores: OLD vs NEW (pad100), D-loop excluded",
                 fontsize=13,fontweight="bold",y=.996)
    fig.tight_layout(rect=[0,0,1,.944])
    stem=f"gene_vs_chrm_ci_old_vs_new_pad100_no_dloop_{a.metric}"
    for ext in ("png","pdf"):
        fig.savefig(out/f"{stem}.{ext}",dpi=220)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(out/f"{stem}_correlations.csv",index=False)
    print(f"Excluded {len(excluded)} D-loop positions; saved plots and CSV to {out}")

if __name__=="__main__":
    main()
