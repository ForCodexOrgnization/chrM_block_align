#!/usr/bin/env python3
"""Four-way 100% OriH containment: longest, GenBank D-loop, GenBank CR, either.

- D-loop: ONLY class=D-loop from current split GenBank CSV.
- GenBank CR: all GenBank candidates through Tier 2 review, including D-loop.
- Either: logical OR of the THREE independent strategies (no region merging).
- OriH: union of MITOS2 parent OriH bases, same GFF rules as prior figures.
- All-species and OriH-predicted denominators match prior all-ranks figures.

Results are provisional until GenBank and MITOS2 assemblies/coordinate origins
are confirmed identical. Multiple GenBank accessions for a species make the
corresponding GenBank method unevaluable and are never considered successes.
"""
import argparse
import collections
import csv
import json
from pathlib import Path

from compare_pairwise_cr_to_longest_intervals import (
    RANKS, read_longest, read_mitos2_paths, read_orih_positions)
from compare_genbank_cr_to_longest_orih import (
    deduplicate, evaluate_regions, evaluate_longest, write_tsv)
from control_region_orih_region_level import parse_dloop, parse_misc

METHODS = ("Longest", "GenBank_Dloop", "GenBank_CR", "Either")


def score_methods(gap, ori, length, dloop_regions, cr_regions):
    """Each annotated CR copy is scored separately; no concatenation."""
    if ori is None or not ori:
        return {method: 0 for method in METHODS}, dict(
            longest_overlap="", dloop_overlap="", cr_overlap="")
    longest_hit, longest_bp = evaluate_longest(gap, ori, length)

    def annotation_score(regions):
        accessions = {r["acc"] for r in regions if r["acc"]}
        if len(accessions) > 1:
            return False, "", "multiple_genbank_accessions"
        yes, _, bp = evaluate_regions(regions, ori, length)
        return bool(yes), bp, "ok" if regions else "no_annotation"

    dloop_hit, dloop_bp, dloop_qc = annotation_score(dloop_regions)
    cr_hit, cr_bp, cr_qc = annotation_score(cr_regions)
    flags = {
        "Longest": int(longest_hit),
        "GenBank_Dloop": int(dloop_hit),
        "GenBank_CR": int(cr_hit),
        "Either": int(bool(longest_hit or dloop_hit or cr_hit)),
    }
    return flags, dict(longest_overlap=longest_bp, dloop_overlap=dloop_bp,
                       cr_overlap=cr_bp, dloop_qc=dloop_qc, cr_qc=cr_qc)


def compare(longest_tsv, dloop_csv, misc_csv, mitos_summary, out_dir, ranks=RANKS):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    gaps = read_longest(longest_tsv)
    gff_paths = read_mitos2_paths(mitos_summary)
    dloop_raw, invalid = parse_dloop(dloop_csv)
    # Same exact-interval de-duplication for both candidates, while retaining
    # D-loop-only provenance in an independent list.
    dloop = deduplicate(dloop_raw)
    cr = deduplicate(dloop_raw + parse_misc(misc_csv))
    ori_cache = {}
    species_output = []
    aggregate = []
    for rank in ranks:
        rows = []
        for (row_rank, key), gap in gaps.items():
            if row_rank != rank:
                continue
            length = int(gap.get("sequence_length") or 0)
            status = "missing_gff"
            ori = None
            gff_value = gff_paths.get(key)
            if not length:
                status = "missing_genome_length"
            elif gff_value:
                gff = Path(gff_value)
                if not gff.is_file():
                    status = "gff_not_found"
                else:
                    cache_key = (str(gff.resolve()), length)
                    if cache_key not in ori_cache:
                        try:
                            intervals, positions = read_orih_positions(gff, length)
                            ori_cache[cache_key] = (
                                "predicted" if intervals else "no_OriH", positions)
                        except (OSError, ValueError) as exc:
                            ori_cache[cache_key] = ("invalid_gff:" + str(exc), None)
                    status, ori = ori_cache[cache_key]
            dregs = dloop.get(key, [])
            cregs = cr.get(key, [])
            flags, bp = score_methods(gap, ori if status == "predicted" else None,
                                      length, dregs, cregs)
            rows.append(dict(
                rank=rank, species=gap.get("species", key), organism_key=key,
                genome_bp=length, OriH_status=status,
                OriH_total_bp=len(ori) if ori is not None else "",
                n_dloop_regions=len(dregs), n_genbank_cr_regions=len(cregs),
                n_dloop_accessions=len({r["acc"] for r in dregs if r["acc"]}),
                n_cr_accessions=len({r["acc"] for r in cregs if r["acc"]}),
                dloop_qc=bp.get("dloop_qc", "not_evaluable"),
                cr_qc=bp.get("cr_qc", "not_evaluable"),
                coordinate_status="provisional_unverified",
                **{f"{method}_full":flags[method] for method in METHODS},
                Longest_overlap_bp=bp["longest_overlap"],
                GenBank_Dloop_max_single_region_overlap_bp=bp["dloop_overlap"],
                GenBank_CR_max_single_region_overlap_bp=bp["cr_overlap"]
            ))
        if not rows:
            raise ValueError("No rank records in longest TSV: " + rank)
        species_output.extend(rows)
        n_pred = sum(r["OriH_status"] == "predicted" for r in rows)
        for denominator in ("All_species", "OriH_predicted_species"):
            subset = rows if denominator == "All_species" else [
                r for r in rows if r["OriH_status"] == "predicted"]
            n = len(subset)
            aggregate.append(dict(
                rank=rank, denominator=denominator,
                n_denominator=n, n_all_species=len(rows), n_OriH_predicted=n_pred,
                **{f"n_{method}":sum(row[f"{method}_full"] for row in subset)
                   for method in METHODS},
                **{f"pct_{method}":round(
                    100 * sum(row[f"{method}_full"] for row in subset) / n, 4
                ) if n else "" for method in METHODS}
            ))
    write_tsv(out / "species_100pct_containment.tsv", species_output,
              list(species_output[0]))
    write_tsv(out / "rank_100pct_containment.tsv", aggregate,
              list(aggregate[0]))
    summary = {
        "comparison": list(METHODS),
        "genbank_dloop": "Only GenBank class=D-loop records; reconstructed split parts",
        "genbank_cr": "GenBank classes through Tier 2 review, including D-loop",
        "either": "Logical OR, without combining two partial regions",
        "n_invalid_split_dloop": len(invalid),
        "warning": "GenBank vs MITOS2 assembly identity and circular coordinate orientation are not verified.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return species_output, aggregate


def plot(aggregate, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    out = Path(out_dir)
    labels = ("Longest", "GenBank D-loop", "GenBank CR", "Either")
    colors = ("#2778B5", "#F28B30", "#57A05B", "#D62929")
    for denominator, subtitle, filename in (
        ("All_species", "All species", "orih_dloop_cr_all_species"),
        ("OriH_predicted_species", "MITOS2 OriH-predicted species",
         "orih_dloop_cr_orih_predicted")
    ):
        data = {r["rank"]: r for r in aggregate if r["denominator"] == denominator}
        ranks = [rank for rank in RANKS if rank in data]
        x = np.arange(len(ranks))
        fig, ax = plt.subplots(figsize=(12.8, 5.5))
        # Three dedicated bands: title / subtitle / legend; no overlap.
        fig.subplots_adjust(left=.09, right=.985, bottom=.24, top=.72)
        width = .19
        for i, (method, label, color) in enumerate(zip(METHODS, labels, colors)):
            y = [float(data[rank][f"pct_{method}"]) for rank in ranks]
            ax.bar(x + (i - 1.5) * width, y, width, color=color, label=label)
        ax.set_ylim(0, 100)
        ax.set_yticks(range(0, 101, 20))
        ax.set_yticklabels([f"{y}%" for y in range(0, 101, 20)])
        ax.set_xticks(x)
        ax.set_xticklabels(ranks, rotation=35, ha="right")
        ax.set_ylabel("All species (%)" if denominator == "All_species"
                      else "MITOS2 OriH-predicted species (%)")
        ax.set_xlabel("Taxonomic rank")
        ax.grid(axis="y", linestyle="--", alpha=.25)
        ax.set_axisbelow(True)
        fig.text(.09, .95, "OriH containment across taxonomic ranks",
                 fontsize=15, fontweight="bold", ha="left", va="top")
        fig.text(.09, .89, "100% OriH coverage | " + subtitle,
                 fontsize=12, ha="left", va="top")
        fig.legend(*ax.get_legend_handles_labels(), loc="upper center",
                   bbox_to_anchor=(.5, .83), ncol=4, frameon=False, fontsize=10)
        for suffix in (".png", ".pdf"):
            fig.savefig(out / (filename + suffix),
                        dpi=220 if suffix == ".png" else None)
        plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--longest-tsv", default="longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv")
    p.add_argument("--dloop", default="data/kingdom_genbank_annotations_split_new.csv")
    p.add_argument("--misc", default="data/genbank_misc_ftr.csv")
    p.add_argument("--mitos-summary", default="data/mitos2_species_summary.tsv")
    p.add_argument("--out-dir", default="longest_all_gene_desert_qc/genbank_Dloop_CR_comparison")
    p.add_argument("--ranks", nargs="+", choices=RANKS, default=list(RANKS))
    args = p.parse_args()
    _, aggregate = compare(args.longest_tsv, args.dloop, args.misc,
                           args.mitos_summary, args.out_dir, args.ranks)
    plot(aggregate, args.out_dir)


if __name__ == "__main__":
    main()
