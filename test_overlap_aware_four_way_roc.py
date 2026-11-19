#!/usr/bin/env python3
"""Regression checks for Original / Mean / Min / Max paired AUROC."""
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd

from compare_overlap_aware_four_way_roc import evaluate, score_table


def test_evaluation():
    df = pd.DataFrame({
        "POS":[100,101,102,103,104,105,106,107,108,109],
        "label":[1,0,1,0,1,0,1,0,1,0],
        "original":[90,10,80,20,75,25,85,15,70,30],
        "mean":[90,10,80,20,75,25,85,15,70,30],
        "min":[90,10,80,20,75,25,85,15,70,30],
        "max":[90,10,80,20,75,25,85,15,70,30],
    })
    result, draws, _ = evaluate(df,"class","genome_wide",30,1,5)
    assert result["status"] == "ok"
    assert result["n_pathogenic"] == 5 and result["n_benign"] == 5
    assert abs(result["auc_original"]-1.0)<1e-9
    assert abs(result["delta_mean_minus_original"])<1e-9
    assert len(draws)>0
    insufficient = evaluate(df.loc[df.label==1],"class","gene_gene_overlap",10,1,5)[0]
    assert insufficient["status"] == "insufficient_labels"


def test_score_table():
    with tempfile.TemporaryDirectory() as td:
        path=Path(td)/"score.csv"
        pd.DataFrame({
            "Position":range(1,16570),
            "Homo sapiens":["A"]*16569,
            "CI_MSA":[50.0]*16569,
        }).to_csv(path,index=False)
        d=score_table(path,"CI_MSA","original")
        assert len(d)==16569
        assert np.all(d.original==50)


if __name__=="__main__":
    test_evaluation()
    test_score_table()
    print("PASS overlap-aware paired ROC regression tests")
