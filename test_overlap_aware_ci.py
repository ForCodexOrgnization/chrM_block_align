#!/usr/bin/env python3
"""Regression checks for overlap-aware CI (no MAFFT or Biopython required)."""
import csv
import tempfile
from pathlib import Path

from calculate_overlap_aware_ci import (
    aggregate, compute_overlap, context_ci, core_columns, PATTERN,
)
from chrM_block_align import Feature, rc

def test():
    species = {"homosapiens", "sp1", "sp2", "sp3"}
    human = "ACGT"
    baseline = [
        {"Position": str(i), "Homo sapiens": base, "CI_MSA": "25.00000000",
         "CI_callable": "25.00000000"}
        for i, base in enumerate(human, 1)
    ]
    at_pos2 = {
        "A": {2: {"homosapiens": "C", "sp1": "C", "sp2": "C", "sp3": "T"}},
        "B": {2: {"homosapiens": "C", "sp1": "C", "sp2": "T", "sp3": "-"}},
    }
    details, pos, scores = compute_overlap(
        at_pos2, {2: ["A", "B"]}, human, species, baseline)
    assert len(details) == 2
    assert len(pos) == 4
    assert [r["Position"] for r in scores["mean"]] == [1,2,3,4]
    # A = 3/4=75%, B = 2/4=50% including unavailable species.
    assert abs(float(scores["mean"][1]["CI_MSA"])-62.5) < 1e-6
    assert abs(float(scores["min"][1]["CI_MSA"])-50.0) < 1e-6
    assert abs(float(scores["max"][1]["CI_MSA"])-75.0) < 1e-6
    # Callable-only: A = 3/4=75%, B = 2/3=66.67%.
    assert abs(float(scores["mean"][1]["CI_callable"])-(75+200/3)/2) < 1e-6
    assert all(scores[m][0]["CI_MSA"] == "25.00000000" for m in scores)
    assert all(scores[m][2]["CI_callable"] == "25.00000000" for m in scores)

    # All contexts uncallable: callable-only falls back to prior score.
    no_calls = {"A": {2: {"homosapiens": "-", "sp1": "-"}},
                "B": {2: {"homosapiens": "-", "sp1": "-"}}}
    _, _, missing = compute_overlap(no_calls, {2:["A","B"]}, human,
                                    species, baseline)
    assert all(missing[m][1]["CI_callable"] == "25.00000000" for m in missing)

    # Verify human negative-strand MSA is flipped back to reference orientation.
    with tempfile.TemporaryDirectory() as td:
        d=Path(td)
        raw="ACGTACGT"
        feature=Feature("G",2,6,"-")
        path=d/"01_G_core.aln.fa"
        path.write_text(">Homo sapiens\\n"+rc(raw[2:6])+
                        "\\n>sp1\\n"+rc(raw[2:6])+"\\n")
        # Core sequence happens to be palindromic here.
        positions=core_columns(path,"G",feature,raw,{"homosapiens","sp1"})
        assert set(positions) == {3,4,5,6}
        assert positions[3]["sp1"] == raw[2]
        assert PATTERN.fullmatch("01_G_core.aln.fa")
        assert PATTERN.fullmatch("01_G_pad50.aln.fa")
    assert aggregate([None, 20., 40.], "mean") == 30.
    print("PASS overlap-aware Mean/Min/Max, baseline preservation, callable handling, MSA projection")

if __name__ == "__main__":
    test()
