#!/usr/bin/env python3
"""QC gene-gene overlap conservation from cached gene-block MSAs (no MAFFT).

Reports gene-context calls, callability and gap fractions, and per-position
shared-species comparisons. Uses exactly the same human-coordinate projection
as calculate_overlap_aware_ci.py. Padding is alignment context, not a second
gene annotation.
"""
import argparse
import collections
import csv
import json
from pathlib import Path

from calculate_overlap_aware_ci import (
    BASES, PATTERN, core_columns, read_fasta, read_baseline, write_csv
)
from chrM_block_align import EXPECTED, Feature, fasta
from human_genbank import load_human_genbank


def _pct(n, d):
    return 100.0 * n / d if d else None


def score_pair(position, gene_a, gene_b, calls_a, calls_b, human_base, species):
    """Per-position diagnostics using a fixed species universe and shared ACGT subset."""
    a = {s: calls_a.get(s, "-").upper() for s in species}
    b = {s: calls_b.get(s, "-").upper() for s in species}
    sa = {s for s in species if a[s] in BASES}
    sb = {s for s in species if b[s] in BASES}
    both = sa & sb
    counts = collections.Counter((a[s] in BASES, b[s] in BASES) for s in species)
    a_matches = sum(a[s] == human_base for s in sa)
    b_matches = sum(b[s] == human_base for s in sb)
    a_shared = sum(a[s] == human_base for s in both)
    b_shared = sum(b[s] == human_base for s in both)
    disagree = sum(a[s] != b[s] for s in both)
    return dict(
        Position=position, genes=f"{gene_a};{gene_b}", Human=human_base,
        n_species=len(species), n_callable_a=len(sa), n_callable_b=len(sb),
        n_gap_a=sum(a[s] == "-" for s in species),
        n_gap_b=sum(b[s] == "-" for s in species),
        n_only_a=counts[(True, False)],
        n_only_b=counts[(False, True)],
        n_neither=counts[(False, False)],
        n_shared_callable=len(both),
        n_discordant_shared=disagree,
        shared_discordance_pct=_pct(disagree, len(both)),
        CI_MSA_a=_pct(a_matches, len(species)),
        CI_MSA_b=_pct(b_matches, len(species)),
        CI_callable_a=_pct(a_matches, len(sa)),
        CI_callable_b=_pct(b_matches, len(sb)),
        CI_shared_a=_pct(a_shared, len(both)),
        CI_shared_b=_pct(b_shared, len(both)),
        callable_fraction_a=_pct(len(sa), len(species)),
        callable_fraction_b=_pct(len(sb), len(species)),
        shared_fraction=_pct(len(both), len(species)),
    )


def read_core_maps(rank_dir, genes, overlap, human, species):
    needed = {g for names in overlap.values() for g in names}
    block_dir = rank_dir / "gene_alignments"
    if not block_dir.is_dir():
        raise FileNotFoundError(f"Cached gene-block alignments missing: {block_dir}")
    files = collections.defaultdict(list)
    for path in block_dir.glob("*.aln.fa"):
        match = PATTERN.fullmatch(path.name)
        if match and match.group(1) in needed:
            files[match.group(1)].append(path)
    result = {}
    for gene in sorted(needed):
        paths = files.get(gene, [])
        if len(paths) != 1:
            raise ValueError(f"Expected one alignment for {gene}; found: {paths}")
        result[gene] = core_columns(paths[0], gene, genes[gene], human, species)
    return result


def plot_pairs(rows, outdir, threshold):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    groups = collections.defaultdict(list)
    for row in rows:
        groups[row["genes"]].append(row)
    for genes, records in sorted(groups.items()):
        records.sort(key=lambda r: r["Position"])
        pos = [r["Position"] for r in records]
        a, b = genes.split(";")
        fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True, constrained_layout=True)
        for ax, metric, title in (
            (axes[0], "CI_MSA", "CI_MSA (fixed rank denominator)"),
            (axes[1], "CI_callable", "CI_callable (per-gene ACGT denominator)"),
            (axes[2], "callable_fraction", "Callable coverage (% of rank)"),
        ):
            ax.plot(pos, [r[f"{metric}_a"] for r in records], "o-", markersize=3, label=a)
            ax.plot(pos, [r[f"{metric}_b"] for r in records], "o-", markersize=3, label=b)
            ax.set_ylabel(title)
            ax.set_ylim(-5, 105)
            ax.grid(alpha=0.25)
            ax.legend()
        axes[2].plot(pos, [r["shared_fraction"] for r in records],
                     ":", linewidth=1.7, label="Both callable")
        axes[2].legend()
        axes[2].set_xlabel("Human chrM position (1-based)")
        fig.suptitle(f"{genes} | overlap-region alignment QC")
        fig.savefig(outdir / f"{a}_{b}_overlap_qc.png", dpi=160)
        plt.close(fig)


def run(rank_dir, outdir, human_fasta, human_genbank, watch_positions, min_shared):
    rank_dir, outdir = Path(rank_dir), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    baseline = read_baseline(rank_dir / "conservation" / "conservation_ci.csv")
    rank = rank_dir.name
    msa = read_fasta(rank_dir / f"{rank}_stitched_aln.fasta")
    human = msa["homosapiens"]
    ref = list(fasta(human_fasta))
    if len(ref) != 1 or ref[0][1] != human or len(baseline) != len(human):
        raise ValueError("Human reference / baseline CI / stitched MSA mismatch")
    species = set(msa)
    genes = load_human_genbank(human_genbank, human, Feature, EXPECTED)
    overlap = collections.defaultdict(list)
    for name, feature in genes.items():
        for p in range(feature.start + 1, feature.end + 1):
            overlap[p].append(name)
    overlap = {p: sorted(names) for p, names in overlap.items() if len(names) > 1}
    calls = read_core_maps(rank_dir, genes, overlap, human, species)
    rows = []
    for position, names in sorted(overlap.items()):
        for ia, a in enumerate(names):
            for b in names[ia+1:]:
                row = score_pair(position, a, b, calls[a][position],
                                 calls[b][position], human[position-1], species)
                row["watched"] = int(position in watch_positions)
                row["low_shared_coverage"] = int(
                    row["shared_fraction"] is None or row["shared_fraction"] < min_shared
                )
                row["abs_CI_MSA_delta"] = abs(row["CI_MSA_a"] - row["CI_MSA_b"])
                row["abs_CI_callable_delta"] = (
                    abs(row["CI_callable_a"] - row["CI_callable_b"])
                    if row["CI_callable_a"] is not None and row["CI_callable_b"] is not None
                    else None
                )
                rows.append(row)
    if not rows:
        raise ValueError("No gene-gene overlap positions")
    write_csv(outdir / "overlap_pairwise_qc.csv", rows)
    flagged = sorted(
        (r for r in rows if r["watched"] or r["low_shared_coverage"]),
        key=lambda r: (not r["watched"], -r["abs_CI_MSA_delta"], r["Position"])
    )
    if flagged:
        write_csv(outdir / "flagged_overlap_positions.csv", flagged)
    plot_pairs(rows, outdir, min_shared)
    summary = dict(rank=rank, n_species=len(species), n_overlap_positions=len(overlap),
                   n_pair_rows=len(rows), watched_positions=sorted(watch_positions),
                   n_low_shared_coverage=sum(r["low_shared_coverage"] for r in rows),
                   min_shared_pct=min_shared,
                   definitions="Shared means ACGT in both gene-context MSAs for the same species; no MAFFT rerun",
                   caveat="Shared calls are not guaranteed orthologous; inspect suspicious indel columns and annotation provenance")
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rank-dir", required=True)
    p.add_argument("--out-dir")
    p.add_argument("--human-fasta", default="data/Humans_mtdna.fasta")
    p.add_argument("--human-genbank", default="data/sequence.gb")
    p.add_argument("--watch-positions", type=int, nargs="*",
                   default=[4330, 4331, 8527, 8530, 8531, 8560, 8563, 10760])
    p.add_argument("--min-shared-pct", type=float, default=50)
    a = p.parse_args()
    if not 0 <= a.min_shared_pct <= 100:
        p.error("--min-shared-pct must be in [0,100]")
    run(a.rank_dir, a.out_dir or str(Path(a.rank_dir)/"overlap_alignment_qc"),
        a.human_fasta, a.human_genbank, set(a.watch_positions), a.min_shared_pct)


if __name__ == "__main__":
    main()
