#!/usr/bin/env python3
"""Synthetic regression checks for shared-species overlap alignment QC."""
from analyze_overlap_alignment_qc import score_pair

def test_score_pair():
    species={"homosapiens","s1","s2","s3"}
    a={"homosapiens":"A","s1":"A","s2":"C","s3":"-"}
    b={"homosapiens":"A","s1":"T","s2":"-","s3":"A"}
    x=score_pair(8527,"ATP6","ATP8",a,b,"A",species)
    assert x["n_species"]==4
    assert x["n_callable_a"]==3 and x["n_callable_b"]==3
    assert x["n_shared_callable"]==2
    assert x["n_only_a"]==1 and x["n_only_b"]==1
    assert x["n_discordant_shared"]==1
    assert abs(x["CI_MSA_a"]-50)<1e-8
    assert abs(x["CI_callable_a"]-200/3)<1e-8
    assert abs(x["CI_shared_a"]-100)<1e-8
    assert abs(x["CI_shared_b"]-50)<1e-8
    assert abs(x["shared_fraction"]-50)<1e-8
    assert abs(x["shared_discordance_pct"]-50)<1e-8

def test_no_shared():
    x=score_pair(4331,"TRNI","TRNQ",{"h":"A"},{"h":"-"},"A",{"h"})
    assert x["n_shared_callable"]==0
    assert x["shared_discordance_pct"] is None
    assert x["CI_shared_a"] is None
    assert x["CI_shared_b"] is None

if __name__=="__main__":
    test_score_pair()
    test_no_shared()
    print("PASS overlap alignment QC synthetic tests")
