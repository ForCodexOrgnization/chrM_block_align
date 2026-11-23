#!/usr/bin/env python3
"""Tests for multipart GenBank annotations and segmented sequence extraction."""
import csv
import tempfile
from pathlib import Path
import chrM_block_align as base
from run_genbank_trnascan_mitos_fallback import load_genbank_csv

def fixture(rows):
    with tempfile.TemporaryDirectory() as tmp:
        path=Path(tmp)/"split.csv"
        keys=["acc","Organism","Gene","G_locus_tag","start","end","strand","len","full_len","part_complete"]
        with path.open("w",newline="") as fp:
            writer=csv.DictWriter(fp,fieldnames=keys)
            writer.writeheader()
            writer.writerows(rows)
        alias={"nd5":frozenset(("ND5",)),"atp8":frozenset(("ATP8",))}
        result,_,invalid,_=load_genbank_csv(path,alias,[("ND5","+"),("ATP8","+")])
        return result,invalid

def row(gene,locus,start,end,part,total,strand="+",acc="NC_TEST"):
    return dict(acc=acc,Organism="Sample species",Gene=gene,
                G_locus_tag=locus,start=start,end=end,strand=strand,
                len=int(end)-int(start)+1,full_len=total,part_complete=part)

def test_two_part_origin():
    entries,invalid=fixture([
        row("ND5","locusA",17,20,"part 1/2",7),
        row("ND5","locusA",1,3,"part 2/2",7)
    ])
    assert not invalid,invalid
    x=entries["samplespecies"][0]
    assert x[0]=="ND5"
    assert x[6]==((16,20),(0,3)),x
    f=base.Feature("ND5",x[1]-1,x[2],x[3],x[6])
    seq="ABCDEFGHIJKLMNOPQRST"
    assert base.extract_feature(seq,f)=="QRSTABC"
    assert base.extract_feature(seq,f,2,2)=="OPQRSTABCDE"
    assert base.extract_feature(seq,base.Feature("ND5",16,23,"-",x[6]))=="QRSTABC"

def test_three_part_discontiguous():
    entries,invalid=fixture([
        row("ND5","locusB",4,6,"part 1/3",9),
        row("ND5","locusB",15,17,"part 2/3",9),
        row("ND5","locusB",1,3,"part 3/3",9),
    ])
    assert not invalid,invalid
    x=entries["samplespecies"][0]
    f=base.Feature("ND5",x[1]-1,x[2],x[3],x[6])
    seq="ABCDEFGHIJKLMNOPQRST"
    assert base.extract_feature(seq,f)=="DEFOPQABC"
    assert base.extract_feature(seq,f,1,1)=="CDEFOPQABCD"
    assert base.rc(base.extract_feature(seq,f))==base.rc("DEFOPQABC")

def test_incomplete_excluded_and_complete_unchanged():
    entries,invalid=fixture([
        row("ND5","bad",17,20,"part 1/2",7),
        row("ATP8","alone",5,8,"complete",4),
    ])
    assert any(reason=="incomplete_or_conflicting_split_locus" for _,reason in invalid)
    assert [x[0] for x in entries["samplespecies"]]==["ATP8"]
    assert entries["samplespecies"][0][6]==()

if __name__=="__main__":
    test_two_part_origin()
    test_three_part_discontiguous()
    test_incomplete_excluded_and_complete_unchanged()
    print("PASS split GenBank segmentation tests")
