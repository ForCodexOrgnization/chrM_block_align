#!/usr/bin/env python3
"""Small end-to-end regression test for conservation scores."""
import csv
import tempfile
from pathlib import Path
from calculate_conservation import calculate

def test():
    with tempfile.TemporaryDirectory() as td:
        base=Path(td);msa=base/"tiny.fa";out=base/"scores";mapping=base/"map.tsv"
        msa.write_text(">Homo sapiens\nACGT\n>Species A\nAC-T\n>Species B\nATNT\n")
        mapping.write_text("human_position\n1\n2\n3\n4\n")
        calculate(msa,out,mapping,expected_length=4)
        with (out/"conservation_ci.csv").open() as f:
            rows=list(csv.DictReader(f))
        assert len(rows)==4
        assert [round(float(r["CI_MSA"]),4) for r in rows]==[100.,66.6667,33.3333,100.]
        assert [round(float(r["CI_callable"]),4) for r in rows]==[100.,66.6667,100.,100.]
        assert [int(r["n_gap"]) for r in rows]==[0,0,1,0]
        assert [int(r["n_ambiguous"]) for r in rows]==[0,0,1,0]
        with (out/"nucleotide_frequencies.csv").open() as f:
            freq=list(csv.DictReader(f))
        assert round(float(freq[2]["g"]),4)==50
        assert round(float(freq[2]["-"]),4)==50
        print("PASS conservation CI and nucleotide frequencies")

if __name__=="__main__":test()
