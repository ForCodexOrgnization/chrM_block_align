#!/usr/bin/env python3
"""Basic circular flank-context checks for experimental padded pipeline."""
from chrM_block_align import Feature
from chrM_block_align_padded import neighbours,ordered

def main():
    names=["TRNT","TRNP","TRNF","RNR1"]
    fs={name:Feature(name,i*100,(i+1)*100,"+") for i,name in enumerate(names)}
    assert ordered(fs,400)==names
    assert neighbours(fs,"TRNP",fs,names,400)==(True,True)
    assert neighbours(fs,"TRNF",fs,names,400)==(True,True)
    # Local gene rearrangement invalidates both affected flanks.
    swap={"TRNT":fs["TRNT"],"TRNP":Feature("TRNP",200,300,"+"),
          "TRNF":Feature("TRNF",100,200,"+"),"RNR1":fs["RNR1"]}
    assert neighbours(swap,"TRNP",fs,names,400)==(False,True)
    assert neighbours(swap,"TRNF",fs,names,400)==(True,False)
    print("PASS padded flank context checks")

if __name__=="__main__":main()
