#!/usr/bin/env python3
"""Synthetic checks for historical-vs-overlap CI comparison."""
import numpy as np
import pandas as pd
from compare_overlap_old_ci_correlation import compare, groups_for

def test_comparison():
    d=pd.DataFrame({
        "POS":[1,2,3,4,5],"old":[10.,20.,30.,40.,50.],
        "original":[11.,21.,31.,41.,51.],
        "mean":[10.,20.,30.,40.,50.],
        "min":[8.,18.,28.,38.,48.],
        "max":[12.,22.,32.,42.,52.],
        "is_overlap":[1,1,1,0,0],
        "score_changed":[True,True,True,False,False],
        "genes":["ATP6;ATP8"]*3+["",""],
    })
    groups=groups_for(d)
    assert len(groups["gene_gene_overlap"])==3
    assert len(groups["changed_overlap"])==3
    rows, paired=compare(groups["gene_gene_overlap"],"kingdom","CI_MSA","gene_gene_overlap",3)
    assert len(rows)==4 and len(paired)==3
    assert abs(rows[1]["pearson_r"]-1)<1e-8
    assert abs(rows[1]["mae"])<1e-8
    small,_=compare(d.iloc[:2],"kingdom","CI_MSA","test",3)
    assert all(r["status"]=="insufficient_positions" for r in small)

if __name__=="__main__":
    test_comparison()
    print("PASS old-CI correlation synthetic tests")
