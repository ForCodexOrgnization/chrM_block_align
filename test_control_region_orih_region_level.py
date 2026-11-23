#!/usr/bin/env python3
"""Synthetic tests for separate CR1/CR2 region and species summaries."""
import csv
from pathlib import Path
from tempfile import TemporaryDirectory

from control_region_orih_region_level import analyze, classify, cr_label


def write(path, fields, rows, delim=","):
    with open(path, "w", newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter=delim)
        w.writeheader()
        w.writerows(rows)


def test_categories():
    assert classify("control region 1")=="full_control_region"
    assert classify("putative control region 2")=="putative_control_region"
    assert classify("control region CSB-1")=="subfeature"
    assert classify("control region coreTAS sequence")=="subfeature"
    assert classify("remnant control region")=="remnant"
    assert classify("AT-rich region")=="expanded_candidate"
    assert cr_label("control region 1")=="CR1"
    assert cr_label("control region II")=="CR2"


def test_region_level():
    with TemporaryDirectory() as tmp:
        root=Path(tmp)
        d=root/"new.csv"
        write(d,["acc","Organism","class","start","end","part_complete","Note","G_locus_tag"],[
            {"acc":"NC_1","Organism":"Species A","class":"D-loop","start":"10","end":"100","part_complete":"complete","Note":"control region 1"},
            {"acc":"NC_1","Organism":"Species A","class":"D-loop","start":"150","end":"200","part_complete":"complete","Note":"control region 2"},
            {"acc":"NC_2","Organism":"Species B","class":"D-loop","start":"900","end":"1000","part_complete":"part 1/2","Note":"control region","G_locus_tag":"locusB"},
            {"acc":"NC_2","Organism":"Species B","class":"D-loop","start":"1","end":"50","part_complete":"part 2/2","Note":"control region","G_locus_tag":"locusB"}
        ])
        misc=root/"misc.csv"
        write(misc,["Organism","acc","misc_kind","feature_index","Note","location","dloop_like"],[
            {"Organism":"Species A","acc":"NC_1","misc_kind":"other","feature_index":"1","Note":"control region CSB-1","location":"20..30"},
            {"Organism":"Species A","acc":"NC_1","misc_kind":"other","feature_index":"2","Note":"remnant control region","location":"210..230"}
        ])
        gff=root/"test.gff"
        gff.write_text("mt\tMITOS2\trep_origin\t160\t170\t.\t+\t.\tID=oriH;Name=heavy-strand origin\n")
        mitos=root/"mitos.tsv"
        write(mitos,["species","result_gff"],[
            {"species":"Species A","result_gff":"test.gff"},
            {"species":"Species B","result_gff":"test.gff"}
        ],"\t")
        summary=analyze(d,misc,mitos,root/"out")
        assert summary["n_species_with_both_cr1_cr2"]==1
        assert summary["n_invalid_split_dloop"]==0
        with open(root/"out"/"region_overlap.csv") as fh:
            regions=list(csv.DictReader(fh))
        cr1=next(r for r in regions if r["species"]=="Species A" and r["region_label"]=="CR1")
        cr2=next(r for r in regions if r["species"]=="Species A" and r["region_label"]=="CR2")
        assert cr1["overlaps_orih"]=="0"
        assert cr2["overlaps_orih"]=="1"
        assert sum(r["n_segments"]=="2" for r in regions)==1
        with open(root/"out"/"species_overlap.csv") as fh:
            species=list(csv.DictReader(fh))
        a=next(r for r in species if r["species"]=="Species A")
        assert a["n_full_regions_overlapping"]=="1"
        assert a["overlapping_labels"]=="CR2"
        assert a["duplicated_cr_candidate"]=="1"


if __name__=="__main__":
    test_categories()
    test_region_level()
    print("PASS CR1 CR2 region-level QC tests")
