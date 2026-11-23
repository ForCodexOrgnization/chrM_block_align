#!/usr/bin/env python3
"""Synthetic regression tests for control region / OriH interval QC."""
import csv
import tempfile
from pathlib import Path
from check_genbank_control_region_mitos2_orih import (
    intervals_misc, intersection_bp, parse_orih_gff, compare
)


def test_segments():
    row={"location":"join(950..1000,1..50)","start":"1","end":"1000"}
    assert intervals_misc(row)==[(950,1000),(1,50)]
    assert intersection_bp([(950,1000),(1,50)],[(980,1000)])==21
    assert intersection_bp([(1,100),(40,60),(90,110)],[(45,95)])==51
    assert intersection_bp([(1,10)],[(11,20)])==0


def test_small_end_to_end():
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp)
        d=root/"dloop.csv"
        with d.open("w",newline="") as fh:
            w=csv.DictWriter(fh,fieldnames=["class","Organism","start","end","acc"])
            w.writeheader()
            w.writerow(dict(class_="") if False else {"class":"D-loop","Organism":"Species one","start":"900","end":"1000","acc":"AC1"})
            w.writerow({"class":"D-loop","Organism":"Species one","start":"1","end":"50","acc":"AC1"})
        m=root/"misc.csv"
        with m.open("w",newline="") as fh:
            w=csv.DictWriter(fh,fieldnames=["misc_kind","Organism","Note","location","dloop_like","acc"])
            w.writeheader()
            w.writerow({"misc_kind":"other","Organism":"Species two","Note":"control region","location":"10..40","dloop_like":"True","acc":"AC2"})
            w.writerow({"misc_kind":"cdd_region","Organism":"Species three","Note":"control region","location":"10..50","acc":"AC3"})
        g=root/"mitos.gff"
        g.write_text("MT\tMITOS2\trep_origin\t20\t35\t.\t+\t.\tID=oriH;Name=heavy-strand origin\n")
        assert parse_orih_gff(g)==[(20,35,"MT")]
        s=root/"summary.tsv"
        s.write_text("species\tresult_gff\nSpecies one\tmitos.gff\nSpecies two\tmitos.gff\n")
        stats=compare(d,m,s,root/"out")
        assert stats["n_genbank_species"]==2
        assert stats["n_with_overlap"]==1
        assert stats["status_species_counts"]["no_overlap"]==1


if __name__=="__main__":
    test_segments()
    test_small_end_to_end()
    print("PASS GenBank control region x MITOS2 OriH QC tests")
