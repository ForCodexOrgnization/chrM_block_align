#!/usr/bin/env python3
"""Derive Mean/Min/Max overlap-aware CI from existing gene-block MAFFT MSAs.

No MAFFT invocation, no edits to original stitched MSA/CI. On positions
annotated in >=2 human genes, score each gene's MSA independently and
aggregate those context scores by mean/min/max. Else reuse baseline CI
bit-for-bit. Supports core and padded gene alignment directories.

Important: CI_MSA includes gaps/unavailable species in the fixed baseline
species denominator; CI_callable excludes non-ACGT calls independently for
each gene. Missing callable context is excluded from aggregation. These
are gene-context conservation scores, not proof of dual-coding constraint.
"""
import argparse
import collections
import csv
import json
import re
from pathlib import Path

from chrM_block_align import EXPECTED, Feature, circular, rc, norm
from human_genbank import load_human_genbank

BASES = set("ACGT")
PATTERN = re.compile(r"^\\d+_(.+)_(core|pad(\\d+))\\.aln\\.fa$")

def read_fasta(path):
    from chrM_block_align import fasta
    records = dict(fasta(path))
    if not records:
        raise ValueError(f"Empty FASTA: {path}")
    if len({len(s) for s in records.values()}) != 1:
        raise ValueError(f"Unequal MSA lengths: {path}")
    return {norm(name): seq for name, seq in records.items()}

def read_baseline(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or [int(r["Position"]) for r in rows] != list(range(1, len(rows)+1)):
        raise ValueError("Baseline conservation rows must cover consecutive positions")
    return rows

def core_columns(path, gene, feature, human, species):
    """Return reference-position to normalized species-base call mappings."""
    aligned = read_fasta(path)
    if "homosapiens" not in aligned:
        raise ValueError(f"Human missing from {path}")
    match = PATTERN.fullmatch(path.name)
    if not match or match.group(1) != gene:
        raise ValueError(f"Bad gene-block filename: {path}")
    padding = int(match.group(3) or 0)
    # Both baseline and padded stitchers reverse-complement -strand MSAs
    # *after* MAFFT so alignments index human-positive genome coordinates.
    if feature.strand == "-":
        aligned = {name: rc(seq) for name, seq in aligned.items()}
    beg, end = feature.start-padding, feature.end+padding
    expected = circular(human, beg, end)
    if aligned["homosapiens"].replace("-", "") != expected:
        raise ValueError(f"Human gene/block reference mismatch: {gene}, {path}")
    # All call strings are projected onto human coordinates; human insertion
    # columns cannot correspond to a single human position and are excluded.
    coords = iter(j % len(human)+1 for j in range(beg, end))
    result = {}
    for j, nt in enumerate(aligned["homosapiens"]):
        if nt == "-":
            continue
        pos = next(coords)
        if (feature.start <= pos-1 < feature.end):
            result[pos] = {k: s[j] for k, s in aligned.items()}
    if next(coords, None) is not None:
        raise ValueError(f"Human coordinate count mismatch: {path}")
    if len(result) != feature.end-feature.start:
        raise ValueError(f"Incomplete gene core: {gene}")
    extra = set(aligned)-set(species)
    if extra:
        raise ValueError(f"Unexpected species in {path}: {list(extra)[:3]}")
    return result

def context_ci(calls, human_base, species):
    n = len(species)
    if not n:
        raise ValueError("No species in baseline")
    vals = [calls.get(k, "-").upper() for k in species]
    matches = sum(v == human_base for v in vals)
    callable_values = [v for v in vals if v in BASES]
    callable_count = len(callable_values)
    return {"CI_MSA": 100.0*matches/n,
            "CI_callable": (100.0*matches/callable_count if callable_count else None),
            "n_callable": callable_count, "n_match": matches,
            "n_gap": sum(v == "-" for v in vals)}

def aggregate(values, method):
    vals = [v for v in values if v is not None]
    if not vals:
        return None
    if method == "mean":
        return sum(vals)/len(vals)
    if method == "min":
        return min(vals)
    if method == "max":
        return max(vals)
    raise ValueError(method)

def compute_overlap(gene_calls, genes_at_pos, human, species, baseline):
    """Pure scoring core; suitable for synthetic tests."""
    if len(human) != len(baseline):
        raise ValueError("Human length and baseline disagree")
    details, positions, methods = [], [], {m: [] for m in ("mean", "min", "max")}
    for i, base in enumerate(human, 1):
        row = baseline[i-1]
        if row["Homo sapiens"] != base:
            raise ValueError(f"Baseline human base mismatch at {i}")
        genes = genes_at_pos.get(i, ())
        scores = []
        for gene in genes:
            calls = gene_calls[gene][i]
            score = context_ci(calls, base, species)
            scores.append((gene, score))
            details.append(dict(Position=i, Gene=gene, Human=base, **score))
        overlap = len(genes) >= 2
        summary = dict(Position=i, Human=base, n_human_genes=len(genes),
                       genes=";".join(genes), is_overlap=int(overlap))
        for metric in ("CI_MSA", "CI_callable"):
            vals = [s[metric] for _, s in scores]
            for method in methods:
                result = aggregate(vals, method) if overlap else None
                # No callable context: retain baseline and flag as fallback.
                if result is None:
                    original = row[metric]
                    result = float(original) if original not in ("", None) else None
                summary[f"{metric}_{method}"] = result
        summary["n_callable_gene_contexts"] = sum(
            s["CI_callable"] is not None for _,s in scores)
        if overlap:
            ci = [s["CI_callable"] for _,s in scores if s["CI_callable"] is not None]
            summary["CI_callable_range"] = max(ci)-min(ci) if len(ci)>1 else None
        else:
            summary["CI_callable_range"] = None
        positions.append(summary)
        for method in methods:
            r = {"Position": i, "Homo sapiens": base,
                 "CI_MSA": (f'{summary[f"CI_MSA_{method}"]:.8f}' if overlap else row["CI_MSA"]),
                 "CI_callable": (f'{summary[f"CI_callable_{method}"]:.8f}'
                                 if overlap and summary[f"CI_callable_{method}"] is not None
                                 else row["CI_callable"]),
                 "n_human_genes": len(genes), "genes": ";".join(genes),
                 "overlap_aware": int(overlap)}
            methods[method].append(r)
    return details, positions, methods

def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows for {path}")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def run(rank_dir, human_fasta, human_genbank, out_dir):
    rank_dir, out_dir = Path(rank_dir), Path(out_dir)
    baseline = read_baseline(rank_dir/"conservation"/"conservation_ci.csv")
    stitched = read_fasta(rank_dir/(rank_dir.name+"_stitched_aln.fasta"))
    if "homosapiens" not in stitched:
        raise ValueError("Stitched MSA lacks Homo sapiens")
    human = stitched["homosapiens"]
    from chrM_block_align import fasta
    reference = list(fasta(human_fasta))
    if len(reference) != 1 or reference[0][1] != human:
        raise ValueError("Human reference does not match stitched MSA")
    genes = load_human_genbank(human_genbank, human, Feature, EXPECTED)
    overlap = collections.defaultdict(list)
    for gene, f in genes.items():
        for p in range(f.start+1, f.end+1):
            overlap[p].append(gene)
    overlap = {p: sorted(names) for p,names in overlap.items() if len(names)>=2}
    if not overlap:
        raise ValueError("No human overlap positions found")
    block_dir = rank_dir/"gene_alignments"
    if not block_dir.is_dir():
        raise FileNotFoundError(f"No gene MSA directory: {block_dir}; rerun with --keep-intermediate")
    available = collections.defaultdict(list)
    for path in block_dir.glob("*.aln.fa"):
        match = PATTERN.fullmatch(path.name)
        if match and match.group(1) in genes:
            available[match.group(1)].append((0 if match.group(2)=="core" else 1,path))
    needed = {g for gs in overlap.values() for g in gs}
    missing = needed-set(available)
    if missing:
        raise FileNotFoundError(f"Missing gene-block MSAs for {sorted(missing)} in {block_dir}")
    gene_calls = {}
    for gene in needed:
        candidates = sorted(available[gene])
        if len(candidates) > 1:
            # A directory should correspond to one alignment method/padding.
            raise ValueError(f"Multiple gene-block MSAs for {gene}: {[str(x[1]) for x in candidates]}")
        gene_calls[gene] = core_columns(
            candidates[0][1],gene,genes[gene],human,set(stitched))
    details, positions, methods = compute_overlap(
        gene_calls, overlap, human, set(stitched), baseline)
    write_csv(out_dir/"overlap_gene_context_ci.csv",details)
    write_csv(out_dir/"overlap_position_ci.csv",positions)
    for method, rows in methods.items():
        write_csv(out_dir/method/"conservation_ci.csv",rows)
    report = {"rank_dir":str(rank_dir),"baseline_ci":str(rank_dir/"conservation"/"conservation_ci.csv"),
              "n_species":len(stitched),"n_overlap_positions":len(overlap),
              "overlap_genes":sorted(needed),"n_gene_context_rows":len(details),
              "methods":list(methods),
              "ci_definition":"gene-specific CI per existing MSA, aggregated per human overlap position",
              "non_overlap":"baseline CI unchanged","missing_callable":"ignore unavailable gene contexts; fallback baseline if none",
              "note":"For ROC, compare methods on the identical variants. These are gene-context scores, not independently replicated nucleotides."}
    out_dir.mkdir(parents=True,exist_ok=True)
    (out_dir/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rank-dir",required=True)
    p.add_argument("--human-fasta",default="data/Humans_mtdna.fasta")
    p.add_argument("--human-genbank",default="data/sequence.gb")
    p.add_argument("--out-dir",default=None)
    args=p.parse_args()
    out = args.out_dir or str(Path(args.rank_dir)/"overlap_aware_ci")
    run(args.rank_dir,args.human_fasta,args.human_genbank,out)

if __name__=="__main__":
    main()
