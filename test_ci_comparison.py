#!/usr/bin/env python3
"""Smoke test for comparison loader, including variant rows and pending ranks."""
import tempfile
from pathlib import Path
import pandas as pd
from compare_old_new_ci import RANKS,load_comparison

def run():
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);old=root/"old.csv";rdir=root/"results"
        o=pd.DataFrame({"POS":[1,1,2],"REF":["A","A","C"],"ALT":["G","T","A"],
                        **{f"{r}_rot_ref":[50,50,20] for r in RANKS}})
        o.to_csv(old,index=False)
        d=rdir/"class"/"conservation"
        d.mkdir(parents=True)
        pd.DataFrame({"Position":range(1,16570),
                      "Homo sapiens":["A","C"]+["G"]*16567,
                      "CI_MSA":[75,50]+[90]*16567,
                      "CI_callable":[100,60]+[95]*16567,
                      "n_callable":[2]*16569,"n_species":[3]*16569,
                      "n_gap":[1]*16569,"n_ambiguous":[0]*16569}).to_csv(d/"conservation_ci.csv",index=False)
        v,p,summary,ranks=load_comparison(old,rdir)
        assert ranks==["class"],ranks
        assert len(v)==3 and len(p)==2
        assert v["class_new_ci"].tolist()==[75,75,50]
        assert v["class_rot_ref"].tolist()==[50,50,20]
        assert p["class_new_ci"].tolist()==[75,50]
        assert summary.loc[summary["rank"]=="class","n_paired_positions"].iloc[0]==2
        print("PASS comparison parsing and variant/position joins")

if __name__=="__main__":run()
