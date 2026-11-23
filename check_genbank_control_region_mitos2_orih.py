#!/usr/bin/env python3
"""Compare GenBank control regions with MITOS2 OriH genomic intervals.

GenBank start/end and GFF positions are 1-based inclusive. Split intervals
are never replaced by their enclosing span. A positive overlap requires the
same species AND a valid coordinate match on the SAME genome assembly and
rotation; this script cannot independently verify assembly identity.
"""
import argparse
import collections
import csv
import json
from pathlib import Path
import re
from urllib.parse import unquote

CONTROL = re.compile(r"\bcontrol[ -]?region\b|\bd[ -]?loop\b", re.I)
AT_RICH = re.compile(r"\ba[+&/]t[ -]?rich(?: region)?\b|\bat[ -]?rich(?: region)?\b", re.I)
LONG_NONCODING = re.compile(r"\b(?:major|long)[ -]?(?:non[ -]?coding|noncoding) region\b", re.I)
ORIH = re.compile(
    r"(?<![a-z0-9])(?:ori[ _-]?h|o[ _-]?h|oh|"
    r"heavy[ _-]?strand[ _-]?(?:replication[ _-]?)?origin|"
    r"origin[ _-]?of[ _-]?(?:the[ _-]?)?heavy[ _-]?strand)(?![a-z0-9])", re.I
)
INTERVAL = re.compile(r"(?<!\d)(\d+)\s*\.\.\s*(\d+)(?!\d)")


def key(value):
    return re.sub(r"[^a-z0-9]", "", (value or "").casefold())


def intervals_misc(row):
    hits = [(int(a), int(b)) for a, b in INTERVAL.findall(row.get("location") or "")]
    if not hits:
        try:
            hits = [(int(row["start"]), int(row["end"]))]
        except (KeyError, TypeError, ValueError):
            return []
    return [(a, b) for a, b in hits if 1 <= a <= b]


def read_genbank(dloop_csv, misc_csv, expanded=False):
    regions = collections.defaultdict(set)
    sources = collections.defaultdict(set)
    accessions = collections.defaultdict(set)
    display = {}
    with open(dloop_csv, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            if (row.get("class") or "").strip().casefold() != "d-loop":
                continue
            species = key(row.get("Organism"))
            if not species:
                continue
            try:
                a, b = int(row["start"]), int(row["end"])
            except (ValueError, KeyError, TypeError):
                continue
            if a < 1 or b < a:
                continue
            regions[species].add((a, b))
            sources[species].add("D-loop")
            accessions[species].add(row.get("acc", ""))
            display[species] = row.get("Organism", "")
    with open(misc_csv, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            if (row.get("misc_kind") or "").strip() != "other":
                continue
            note = " ".join((row.get(k) or "") for k in
                            ("Note", "Product", "Gene", "other_qualifiers"))
            strict = (row.get("dloop_like") or "").casefold() == "true" or bool(CONTROL.search(note))
            tentative = bool(AT_RICH.search(note) or LONG_NONCODING.search(note))
            if not strict and not (expanded and tentative):
                continue
            species = key(row.get("Organism"))
            spans = intervals_misc(row)
            if not species or not spans:
                continue
            regions[species].update(spans)
            sources[species].add("misc_control_region" if strict else "misc_expanded_candidate")
            accessions[species].add(row.get("acc", ""))
            display[species] = row.get("Organism", "")
    return regions, sources, accessions, display


def parse_orih_gff(path):
    result = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            try:
                start, end = int(parts[3]), int(parts[4])
            except ValueError:
                continue
            annotations = []
            for part in parts[8].split(";"):
                if "=" in part:
                    annotations.append(unquote(part.split("=", 1)[1]))
            description = " ".join([parts[2], *annotations])
            if start >= 1 and end >= start and ORIH.search(description):
                result.add((start, end, parts[0]))
    return sorted(result)


def merged(intervals):
    spans = sorted(set(intervals))
    output = []
    for a, b in spans:
        if output and a <= output[-1][1] + 1:
            output[-1] = (output[-1][0], max(b, output[-1][1]))
        else:
            output.append((a, b))
    return output


def intersection_bp(a, b):
    left = merged(a)
    right = merged(b)
    i = j = total = 0
    while i < len(left) and j < len(right):
        x, y = left[i], right[j]
        total += max(0, min(x[1], y[1]) - max(x[0], y[0]) + 1)
        if x[1] < y[1]:
            i += 1
        else:
            j += 1
    return total


def resolve_gff(summary, value):
    if not value:
        return None
    path = Path(value)
    for candidate in (path, summary.parent / path):
        if candidate.is_file():
            return candidate
    return None


def compare(dloop, misc, mitos_summary, outdir, expanded=False):
    regions, sources, accessions, display = read_genbank(dloop, misc, expanded)
    summary = Path(mitos_summary).resolve()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    mitos = {}
    with open(summary, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            species = key(row.get("species"))
            if species:
                mitos.setdefault(species, row)
    rows = []
    for species in sorted(regions):
        row = mitos.get(species)
        path_string = (row or {}).get("result_gff", "")
        gff = resolve_gff(summary, path_string)
        orih = []
        if row is None:
            status = "absent_from_mitos_summary"
        elif not path_string:
            status = "missing_gff_path"
        elif gff is None:
            status = "gff_not_found"
        else:
            orih = parse_orih_gff(gff)
            status = "no_mitos2_orih" if not orih else "no_overlap"
        overlap_bp = intersection_bp(regions[species],
                                     [(a, b) for a, b, _ in orih]) if orih else 0
        if overlap_bp:
            status = "overlap"
        rows.append({
            "species": display[species], "organism_key": species,
            "genbank_accessions": ";".join(sorted(accessions[species])),
            "genbank_source": ";".join(sorted(sources[species])),
            "genbank_intervals": ";".join(f"{a}-{b}" for a,b in merged(regions[species])),
            "n_genbank_intervals": len(merged(regions[species])),
            "n_orih_features": len(orih),
            "orih_coordinates": ";".join(f"{a}-{b}" for a,b,_ in orih),
            "overlap_bp": overlap_bp, "status": status,
            "mitos2_gff": path_string,
        })
    with open(outdir/"species_overlap.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else
                                ["species","organism_key","status"])
        writer.writeheader()
        writer.writerows(rows)
    status_counts = dict(sorted(collections.Counter(row["status"] for row in rows).items()))
    source_counts = dict(sorted(collections.Counter(
        source for group in sources.values() for source in group).items()))
    metrics = {
        "n_genbank_species": len(regions),
        "n_mitos_summary_species": len(mitos),
        "source_species_counts": source_counts,
        "status_species_counts": status_counts,
        "n_with_mitos_orih": sum(row["n_orih_features"] > 0 for row in rows),
        "n_with_overlap": status_counts.get("overlap", 0),
        "expanded": expanded,
        "warning": "Coordinate overlap is provisional until assembly accession, orientation, and rotation match between the GenBank and MITOS2 inputs."
    }
    (outdir/"summary.json").write_text(json.dumps(metrics, indent=2)+"\n")
    print(json.dumps(metrics, indent=2))
    return metrics


def main():
    """Default CLI now runs region-level + species-level QC.

    Legacy helpers above remain importable for backward compatibility.
    """
    from control_region_orih_region_level import main as region_main
    region_main()


if __name__ == "__main__":
    main()
