#!/usr/bin/env python3
"""Region-level GenBank CR1/CR2 versus MITOS2 OriH QC.

Physical overlaps are provisional: matching species names do not verify that
GenBank and MITOS2 use the same assembly, orientation, or rotation.
"""
import argparse
import collections
import csv
import json
import re
from pathlib import Path

from check_genbank_control_region_mitos2_orih import (
    key, parse_orih_gff, intersection_bp, resolve_gff, intervals_misc
)

SUBFEATURE = re.compile(
    r"\b(?:csb[- _]?[123ivx]*|tas|coretas|termination.associated|"
    r"conserved.sequence.block|conserved.region|c.rich.sequence|"
    r"poly[cgt]|\bcore\b)\b", re.I)
REMNANT = re.compile(r"\b(?:pseudo[- ]?control|remnant|relict|degenerate)\b", re.I)
CONTROL = re.compile(r"\bcontrol[ -]?region\b|\bd[ -]?loop\b", re.I)
EXPANDED = re.compile(
    r"\ba[+&/]t[ -]?rich\b|\bat[ -]?rich\b|"
    r"\b(?:major|long)[ -]?(?:non[ -]?coding|noncoding)\b", re.I)
NUMBERED = re.compile(
    r"\bcontrol[ -]?region\s*(?:copy\s*)?(?:1|2|i|ii)\b|"
    r"\b(?:cr[- _]?[12]|cr[- _]?i{1,2})\b", re.I)
PART = re.compile(r"part\s+(\d+)\s*/\s*(\d+)", re.I)


def classify(text, dloop=False, dloop_like=False):
    """Whole-region annotations are distinct from internal motifs/remnants."""
    if REMNANT.search(text):
        return "remnant"
    if SUBFEATURE.search(text):
        return "subfeature"
    if dloop or CONTROL.search(text) or dloop_like:
        if re.search(r"\b(?:putative|probable|possible|predicted)\b", text, re.I):
            return "putative_control_region"
        return "full_control_region"
    if EXPANDED.search(text):
        return "expanded_candidate"
    return "unrelated"


def cr_label(text):
    # Exact numbered CR labels only; unnumbered duplicated copies stay unknown.
    if re.search(r"\bcontrol[ -]?region\s*(?:copy\s*)?(?:1|i)\b|\bcr[- _]?(?:1|i)\b(?!i)", text, re.I):
        return "CR1"
    if re.search(r"\bcontrol[ -]?region\s*(?:copy\s*)?(?:2|ii)\b|\bcr[- _]?(?:2|ii)\b", text, re.I):
        return "CR2"
    return "CR_unspecified"


def write_csv(path, rows, columns):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=columns)
        w.writeheader()
        w.writerows(rows)


def parse_dloop(path):
    complete = []
    split = collections.defaultdict(list)
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for idx, row in enumerate(csv.DictReader(fh), 2):
            if (row.get("class") or "").strip().casefold() != "d-loop":
                continue
            species = key(row.get("Organism"))
            acc = (row.get("acc") or "").strip()
            if not species:
                continue
            try:
                a, b = int(row["start"]), int(row["end"])
            except (KeyError, TypeError, ValueError):
                continue
            if not 1 <= a <= b:
                continue
            text = " ".join((row.get(x) or "") for x in
                            ("Gene", "alt_name", "alt_Product", "Note"))
            part = (row.get("part_complete") or "complete").strip()
            item = dict(species=species, display=row.get("Organism",""),
                        acc=acc, source="D-loop", source_id=f"line:{idx}",
                        classification=classify(text, dloop=True),
                        cr_label=cr_label(text), note=text, intervals=[(a,b)])
            match = PART.fullmatch(part)
            if match:
                number, total = map(int, match.groups())
                locus = (row.get("G_locus_tag") or "").strip()
                # Parts belong to a single source feature, not all D-loops
                # with the same species. Fail closed if grouping ambiguous.
                token = locus if locus and locus != "???" else (
                    row.get("Gene") or row.get("alt_name") or "").strip()
                split[(species, acc, token, total)].append((number, item))
            elif part.casefold() == "complete":
                complete.append(item)
    invalid=[]
    for group, parts in split.items():
        total=group[-1]
        parts.sort(key=lambda p:p[0])
        if len(parts)!=total or [p[0] for p in parts]!=list(range(1,total+1)):
            invalid.append(dict(group=str(group), problem="incomplete_or_duplicate_parts"))
            continue
        annotations={p["cr_label"] for _,p in parts}
        if len(annotations)!=1:
            invalid.append(dict(group=str(group), problem="conflicting_cr_labels"))
            continue
        first=parts[0][1].copy()
        first["source_id"]="split:"+str(group[2])
        first["intervals"]=[interval for _,p in parts for interval in p["intervals"]]
        complete.append(first)
    return complete, invalid


def parse_misc(path):
    result = []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for idx,row in enumerate(csv.DictReader(fh),2):
            if (row.get("misc_kind") or "").strip().casefold() != "other":
                continue
            text=" ".join((row.get(x) or "") for x in
                          ("Note","Product","Gene","other_qualifiers"))
            classification=classify(
                text, dloop_like=(row.get("dloop_like") or "").strip().casefold()=="true")
            if classification == "unrelated":
                continue
            spans=intervals_misc(row)
            species=key(row.get("Organism"))
            if not species or not spans:
                continue
            result.append(dict(
                species=species,display=row.get("Organism",""),
                acc=(row.get("acc") or "").strip(),
                source="misc_feature",
                source_id=(row.get("feature_index") or f"line:{idx}"),
                classification=classification,
                cr_label=cr_label(text),
                note=text,intervals=spans
            ))
    return result


def load_mitos(path):
    path=Path(path).resolve()
    result=collections.defaultdict(list)
    with open(path,newline="",encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh,delimiter="\t"):
            species=key(row.get("species"))
            if species:
                result[species].append(row)
    return path,result


def analyze(dloop,misc,mitos_summary,out_dir):
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    dloop_regions,invalid=parse_dloop(dloop)
    all_regions=dloop_regions+parse_misc(misc)
    mitos_path,mitos=load_mitos(mitos_summary)
    by_species=collections.defaultdict(list)
    for region in all_regions:
        by_species[region["species"]].append(region)
    full_types={"full_control_region","putative_control_region"}
    region_rows=[]
    species_rows=[]
    for species, regions in sorted(by_species.items()):
        full=[r for r in regions if r["classification"] in full_types]
        expanded=[r for r in regions if r["classification"]=="expanded_candidate"]
        noteonly=[r for r in regions if r["classification"] in {"subfeature","remnant"}]
        accs={r["acc"] for r in full if r["acc"]}
        gs=mitos.get(species,[])
        gff_values={r.get("result_gff","") for r in gs if r.get("result_gff")}
        if len(gff_values)>1:
            status="multiple_mitos_gffs"
            orih=[]
        elif not gs:
            status="absent_from_mitos_summary";orih=[]
        elif not gff_values:
            status="missing_gff_path";orih=[]
        else:
            path=resolve_gff(mitos_path,next(iter(gff_values)))
            if path is None:
                status="gff_not_found";orih=[]
            else:
                orih=parse_orih_gff(path)
                status="no_mitos2_orih" if not orih else "evaluated"
        # Multiple GenBank assemblies cannot safely be combined into one
        # positional comparison. Report separately with an explicit warning.
        if len(accs)>1 and status=="evaluated":
            status="multiple_genbank_accessions"
        overlapping=0
        overlapping_labels=set()
        for index,r in enumerate(regions,1):
            eligible=r["classification"] in full_types
            overlap=None
            if eligible and status=="evaluated":
                overlap=intersection_bp(r["intervals"],
                                        [(a,b) for a,b,_ in orih])
                if overlap:
                    overlapping+=1
                    overlapping_labels.add(r["cr_label"])
            region_rows.append(dict(
                species=r["display"],organism_key=species, accession=r["acc"],
                region_id=f"{r['source']}:{r['source_id']}:{index}",
                region_label=r["cr_label"],region_type=r["classification"],
                source=r["source"],note=r["note"],
                n_segments=len(r["intervals"]),
                intervals=";".join(f"{a}-{b}" for a,b in r["intervals"]),
                n_orih_features=len(orih),
                orih_intervals=";".join(f"{a}-{b}" for a,b,_ in orih),
                overlap_bp="" if overlap is None else overlap,
                overlaps_orih="" if overlap is None else int(overlap>0),
                status=status if eligible else "not_full_control_region"
            ))
        if status=="evaluated":
            status="overlap" if overlapping else "no_overlap"
        species_rows.append(dict(
            species=regions[0]["display"],organism_key=species,
            n_full_regions=len(full),n_expanded_candidates=len(expanded),
            n_subfeatures_or_remnants=len(noteonly),
            n_cr1=sum(r["cr_label"]=="CR1" for r in full),
            n_cr2=sum(r["cr_label"]=="CR2" for r in full),
            duplicated_cr_candidate=int(
                any(r["cr_label"]=="CR1" for r in full) and
                any(r["cr_label"]=="CR2" for r in full)),
            n_mitos_orih=len(orih),n_full_regions_overlapping=overlapping,
            overlapping_labels=";".join(sorted(overlapping_labels)),
            accession_count=len(accs),
            status=status if full else "no_full_control_region"
        ))
    region_columns=[
        "species","organism_key","accession","region_id","region_label",
        "region_type","source","note","n_segments","intervals",
        "n_orih_features","orih_intervals","overlap_bp","overlaps_orih","status"]
    species_columns=[
        "species","organism_key","n_full_regions","n_expanded_candidates",
        "n_subfeatures_or_remnants","n_cr1","n_cr2",
        "duplicated_cr_candidate","n_mitos_orih",
        "n_full_regions_overlapping","overlapping_labels","accession_count","status"]
    write_csv(out/"region_overlap.csv",region_rows,region_columns)
    write_csv(out/"species_overlap.csv",species_rows,species_columns)
    write_csv(out/"invalid_split_dloop.csv",invalid,["group","problem"])
    summary={
        "n_species_with_any_candidate_annotation":len(species_rows),
        "n_species_with_full_or_putative_control_region":sum(
            r["n_full_regions"]>0 for r in species_rows),
        "n_species_with_both_cr1_cr2":sum(
            r["duplicated_cr_candidate"] for r in species_rows),
        "n_full_regions":sum(r["n_full_regions"] for r in species_rows),
        "n_invalid_split_dloop":len(invalid),
        "region_type_counts":dict(collections.Counter(
            r["region_type"] for r in region_rows)),
        "species_status_counts":dict(collections.Counter(
            r["status"] for r in species_rows)),
        "coordinate_warning":"Overlap requires same underlying assembly, orientation and circular rotation. Species name alone does not establish this. Multiple GenBank accessions are excluded from positional evaluation.",
        "functional_warning":"OriH overlap does not demonstrate replication initiation activity."
    }
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dloop",default="data/kingdom_genbank_annotations_split_new.csv")
    p.add_argument("--misc",default="data/genbank_misc_ftr.csv")
    p.add_argument("--mitos-summary",default="data/mitos2_species_summary.tsv")
    p.add_argument("--out-dir",default="control_region_orih_qc/region_level")
    args=p.parse_args()
    analyze(args.dloop,args.misc,args.mitos_summary,args.out_dir)

if __name__=="__main__":
    main()
