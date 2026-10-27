#!/usr/bin/env python3
"""Dependency-free smoke tests for core normalization and coordinate semantics."""
from chrM_block_align import *
def test():
    assert len(EXPECTED)==37
    assert canonical_mitos("gene_nad4l_1")=="ND4L"
    assert canonical_mitos("rrnS")=="RNR1"
    assert canonical_trna("Ser","TGA")=="TRNS1"
    assert canonical_trna("Ser","GCT")=="TRNS2"
    assert canonical_trna("Leu","TAA")=="TRNL1"
    assert canonical_trna("Leu","TAG")=="TRNL2"
    assert rc("ACGT")=="ACGT"
    assert circular("ACGT",3,6)=="TAC"
    fs,amb=pick_features([Feature("ATP8",1,6,"+"),Feature("ATP8",1,6,"+")],20)
    assert "ATP8" in fs and not amb
    fs,amb=pick_features([Feature("ATP8",1,6,"+"),Feature("ATP8",12,17,"+")],20)
    assert "ATP8" in amb
    print("PASS: naming, circular extraction, duplicated locus tests")
if __name__=="__main__":test()
