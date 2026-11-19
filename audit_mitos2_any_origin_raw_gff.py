#!/usr/bin/env python3
"""Audit raw Kingdom MITOS2 GFF origin features and longest gene-free intervals.

Run from the chrM_block_align repository on the HPC system with MITOS2 GFFs
mounted. Outputs CSV files; no changes are made to upstream annotations.
"""
import argparse
import csv
import re
from collections import Counter
from pathlib import Path

def norm(value):
    return re.sub(r"[^a-z0-9]", "", str(value).casefold())

def read_tsv(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def attributes(raw):
    result = {}
    for entry in raw.split(";"):
        if "=" in entry:
            key, val = entry.split("=", 1)
            result[key.lower()] = val.strip()
    return result

def classify(feature_type, attrs):
    """Classify origin parent features; never count child transcript/exon rows."""
    ft = norm(feature_type)
    if ft in {"exon", "mrna", "transcript", "cds", "gene"}:
        return None
    origin_types = {"originofreplication", "reporigin",
                    "replicationorigin", "oriregion"}
    vals = {norm(attrs.get(key, "")) for key in
            ("name", "gene_id", "id", "product")}
    vals.discard("")
    oh = {"oh", "orih", "geneoh", "transcriptoh",
          "heavystrandorigin", "originofheavystrandreplication"}
    ol = {"ol", "oril", "geneol", "transcriptol",
          "lightstrandorigin", "originoflightstrandreplication"}
    generic = {"ori", "origin", "geneori", "transcriptori",
               "replicationorigin", "reporigin"}
    if ft not in origin_types and not (vals & (oh | ol | generic)):
        return None
    if vals & oh or any(re.fullmatch(r"oh[0-9]+", v) for v in vals):
        return "OriH"
    if vals & ol or any(re.fullmatch(r"ol[0-9]+", v) for v in vals):
        return "OriL"
    return "generic_or_unclassified_ori"

def segments(start, end, length):
    start, end, length = int(start), int(end), int(length)
    if end == 0:
        end = length  # longest-gap script uses 0 as genome-end sentinel
    if not (1 <= start <= length and 1 <= end <= length):
        raise ValueError("coordinates outside the mitochondrial genome")
    return ([(start, end)] if start <= end else
            [(start, length), (1, end)])

def interval_bp(start, end, length):
    return sum(b - a + 1 for a, b in segments(start, end, length))

def overlap_bp(a, b, c, d, length):
    return sum(max(0, min(y, v) - max(x, u) + 1)
               for x, y in segments(a, b, length)
               for u, v in segments(c, d, length))

def write_csv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", default="data/mitos2_species_summary.tsv")
    ap.add_argument("--longest", default=(
        "longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv"))
    ap.add_argument("--outdir", default="longest_all_gene_desert_qc/any_origin_audit")
    ap.add_argument("--gff-root", default=None,
                    help="Optional replacement root for code_*/species/result.gff")
    args = ap.parse_args()
    output = Path(args.outdir)

    gaps = {norm(row["species"]): row for row in read_tsv(args.longest)
            if row["rank"] == "kingdom"}
    if not gaps:
        raise ValueError("No Kingdom species in longest-interval table")
    gff_index = {}
    for row in read_tsv(args.summary):
        key = norm(row.get("species", ""))
        path = (row.get("result_gff") or "").strip()
        if key and path:
            if key in gff_index and gff_index[key] != path:
                raise ValueError("Conflicting GFF paths for " + key)
            gff_index[key] = path

    features = []
    species = []
    missing = []
    feature_types = Counter()
    for key, gap in gaps.items():
        name = gap["species"]
        raw_path = gff_index.get(key, "")
        path = Path(raw_path) if raw_path else None
        if path is not None and args.gff_root:
            match = re.search(r"(code_[^/]+/[^/]+/result\.gff)$", raw_path)
            if match:
                path = Path(args.gff_root) / match.group(1)
        status = "read"
        ori = []
        if path is None or not path.is_file():
            status = "missing_gff"
            missing.append({"species": name, "gff_path": raw_path})
        else:
            try:
                with path.open(encoding="utf-8", errors="replace") as handle:
                    for line in handle:
                        if not line.strip() or line.startswith("#"):
                            continue
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 9:
                            continue
                        feature_types[parts[2]] += 1
                        attrs = attributes(parts[8])
                        category = classify(parts[2], attrs)
                        if category is None:
                            continue
                        ori.append(dict(species=name, feature_type=parts[2],
                            category=category, start=int(parts[3]),
                            end=int(parts[4]), attributes=parts[8],
                            gff_path=str(path)))
            except (OSError, ValueError) as error:
                status = "error_" + type(error).__name__
                ori = []
                missing.append({"species": name, "gff_path": str(path)})

        # The longest-gap table uses sequence_length (not genome_bp).
        length = int(gap["sequence_length"])
        start = gap.get("gap_start_1based", "")
        end = gap.get("gap_end_1based", "")
        valid = gap.get("candidate_defined") == "1" and bool(start and end)
        for origin in ori:
            try:
                origin["origin_length_bp"] = interval_bp(
                    origin["start"], origin["end"], length)
                overlap = (overlap_bp(start, end, origin["start"],
                                      origin["end"], length)
                           if valid else None)
                origin["longest_overlap_bp"] = overlap if overlap is not None else ""
                origin["fully_contained"] = (
                    int(overlap == origin["origin_length_bp"])
                    if overlap is not None else "")
            except ValueError:
                origin["origin_length_bp"] = ""
                origin["longest_overlap_bp"] = ""
                origin["fully_contained"] = ""
            features.append(origin)
        categories = {origin["category"] for origin in ori}
        if status != "read":
            group = status
        elif not categories:
            group = "none"
        elif categories == {"OriH", "OriL"}:
            group = "OriH+OriL"
        elif len(categories) == 1:
            group = next(iter(categories)) + "_only"
        else:
            group = "+".join(sorted(categories))
        evaluable = [o for o in ori if isinstance(o["longest_overlap_bp"], int)]
        species.append(dict(species=name, status=status, origin_class=group,
            n_origin_parent_features=len(ori),
            n_OriH=sum(o["category"] == "OriH" for o in ori),
            n_OriL=sum(o["category"] == "OriL" for o in ori),
            n_generic=sum(o["category"] == "generic_or_unclassified_ori" for o in ori),
            longest_gap_bp=gap.get("longest_gap_bp", ""),
            interval_available=int(valid),
            any_origin_overlap=(int(any(o["longest_overlap_bp"] > 0
                                        for o in evaluable)) if evaluable else ""),
            any_origin_full=(int(any(o["fully_contained"] == 1
                                     for o in evaluable)) if evaluable else ""),
            all_origins_full=(int(all(o["fully_contained"] == 1
                                      for o in evaluable)) if evaluable else ""),
            gff_path=str(path) if path else ""))

    if len(species) != len(gaps):
        raise AssertionError("Species accounting mismatch")
    write_csv(output / "kingdom_origin_features_all_classes.csv", features,
              ["species", "feature_type", "category", "start", "end",
               "origin_length_bp", "longest_overlap_bp", "fully_contained",
               "attributes", "gff_path"])
    write_csv(output / "kingdom_species_any_origin_summary.csv", species,
              list(species[0]))
    no_oh = [row for row in species if row["n_OriH"] == 0]
    write_csv(output / "kingdom_species_no_OriH_raw_audit.csv", no_oh,
              list(species[0]))
    write_csv(output / "missing_gff.csv", missing,
              ["species", "gff_path"])
    print("Kingdom species:", len(species),
          "| raw GFFs read:", sum(r["status"] == "read" for r in species))
    print("Origin classification:", dict(Counter(r["origin_class"] for r in species)))
    print("Top raw GFF feature types:", feature_types.most_common(20))
    print("Origin feature labels:", dict(Counter(r["category"] for r in features)))
    print("Without OriH:", len(no_oh))
    print("Drosophila melanogaster:",
          [r for r in species if norm(r["species"]) == "drosophilamelanogaster"])
    print("Output:", output)

if __name__ == "__main__":
    main()
