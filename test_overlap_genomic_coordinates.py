#!/usr/bin/env python3
"""Regression tests for coordinate mapping without MAFFT."""
import tempfile
from pathlib import Path
from types import SimpleNamespace
from analyze_overlap_genomic_coordinates import (
    source_coordinates, project_columns, classify, physical_overlap
)
from chrM_block_align import Feature, rc

def test_coordinates():
    f=Feature("ATP6",2,6,"+")
    assert source_coordinates(f,8,0,True,True,"+")==[3,4,5,6]
    assert source_coordinates(Feature("ATP6",2,6,"-"),8,0,True,True,"+")==[6,5,4,3]
    assert source_coordinates(Feature("ATP6",2,6,"-"),8,0,True,True,"-")==[3,4,5,6]

def test_projection():
    with tempfile.TemporaryDirectory() as tmp:
        seq="AACGTGCA"
        path=Path(tmp)/"01_ATP6_pad100.aln.fa"
        # Human block wraps around circular boundaries, species has the same
        # reference sequence and a plus-strand feature.
        h=Feature("ATP6",0,4,"+")
        human="A"*300
        species="A"*300
        path.write_text(">Homo sapiens\n"+("A"*204)+"\n>sp1\n"+("A"*204)+"\n")
        calls,_=project_columns(
            path,"ATP6",h,{"sp1":{"ATP6":Feature("ATP6",0,4,"+")}},
            {"sp1":species},human,{1,2},lambda s,g:(True,True),100)
        assert calls[1]["sp1"]==("A",1)
        assert calls[2]["sp1"]==("A",2)

def test_categories():
    a=Feature("A",0,5,"+")
    b=Feature("B",3,8,"+")
    assert physical_overlap(a,b,10)
    assert classify(a,b,10,"A",4,"A",4)=="same_coordinate"
    assert classify(a,b,10,"A",4,"C",5)=="different_coordinate_with_gene_overlap"
    assert classify(a,Feature("B",6,9,"+"),10,"A",4,"C",8)=="no_physical_gene_overlap"
    assert classify(a,b,10,"-",None,"A",5)=="not_jointly_callable"
    assert classify(None,b,10,None,None,"A",5)=="missing_annotation"

if __name__=="__main__":
    test_coordinates()
    test_projection()
    test_categories()
    print("PASS coordinate overlap QC")
