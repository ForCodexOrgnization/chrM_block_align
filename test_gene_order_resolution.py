#!/usr/bin/env python3
"""Small no-data tests for circular gene-order resolution."""
from resolve_gene_order_annotations import resolve,equivalent,candidates
def test():
    ref=[("A","+"),("B","+"),("C","-")]
    assert equivalent(ref,[("B","+"),("C","-"),("A","+")])
    rows=[({"A"},"+",1),({"B","D"},"+",2),({"C"},"-",3)]
    chosen,qc=resolve(rows,ref)
    assert chosen==["A","B","C"],(chosen,qc)
    assert qc[0]["status"]=="unique_exact_order"
    ambiguous=[({"A"},"+",1),({"B","D"},"+",2)]
    chosen,qc=resolve(ambiguous,ref)
    assert chosen[1] is None
    assert qc[0]["status"] in ("inferred_best_fit_not_accepted","ambiguous_best_fit")
    names={"trnaleu":frozenset(["TRNL1","TRNL2"])}
    assert candidates({"Gene":"tRNA-Leu"},names)=={"TRNL1","TRNL2"}
    print("PASS gene order resolution tests")
if __name__=="__main__":test()
