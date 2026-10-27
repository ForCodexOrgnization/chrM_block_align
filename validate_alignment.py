#!/usr/bin/env python3
"""Validate a human-coordinate stitched MSA."""
import argparse,csv,json
from pathlib import Path
from chrM_block_align import fasta,norm
def validate(rank_dir,human_fasta):
    d=Path(rank_dir)
    expected=list(fasta(human_fasta))
    assert len(expected)==1 and len(expected[0][1])==16569
    msa=list(fasta(d/(d.name+"_stitched_aln.fasta")))
    assert msa,"Empty alignment"
    lengths={len(seq) for _,seq in msa}
    assert len(lengths)==1,"MSA records have different lengths"
    humans=[seq for name,seq in msa if norm(name)=="homosapiens"]
    assert len(humans)==1,"Expected one human"
    assert humans[0]==expected[0][1],"Human alignment sequence differs from canonical reference"
    with open(d/"human_coordinate_map.tsv") as f:positions=[int(r["human_position"]) for r in csv.DictReader(f,delimiter="\t")]
    assert positions==list(range(1,16570)),"Human mapping must be contiguous 1..16569"
    assert len(msa[0][1])==16569,"Unexpected alignment width"
    info=dict(pass_validation=True,rank=d.name,n_species=len(msa),human_bases=16569)
    (d/"validation.json").write_text(json.dumps(info,indent=2)+"\n")
    print("PASS",info)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--rank-dir",required=True)
    p.add_argument("--human-fasta",required=True)
    a=p.parse_args();validate(a.rank_dir,a.human_fasta)
