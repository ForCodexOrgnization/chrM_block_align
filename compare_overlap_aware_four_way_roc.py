#!/usr/bin/env python3
"""Compare pad100 Original vs Mean/Min/Max overlap-aware CI AUROC.

Reads the selected tRNAscan > GenBank > MITOS2 fallback pad100 outputs
without changing alignments or conservation tables. Applies existing
ClinVar/MITOMAP binary labels, identical complete-case variants per rank
and per group, and paired position-cluster bootstrap for AUC differences.
Raw conservation scores are used (no automatic score direction inversion).
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

from analyze_ci_diagnostics import POSITIVE
from compare_old_new_ci import RANKS

METHODS = ("original", "mean", "min", "max")
COMPARISONS = (("mean", "original"), ("min", "original"), ("max", "original"))
OUTPUT_NAME = "overlap_aware_roc"


def score_table(path, metric, method):
    if not path.is_file():
        raise FileNotFoundError(path)
    data = pd.read_csv(path)
    required = {"Position", "Homo sapiens", metric}
    if not required.issubset(data):
        raise ValueError(f"{path}: missing {sorted(required - set(data))}")
    if (len(data) != 16569 or
            not data["Position"].tolist() == list(range(1, 16570))):
        raise ValueError(f"{path}: expected ordered 16569 human positions")
    values = pd.to_numeric(data[metric], errors="coerce")
    if not values.dropna().between(0, 100).all():
        raise ValueError(f"{path}: {metric} must be 0..100")
    return pd.DataFrame({
        "POS": data["Position"].to_numpy(),
        method: values.to_numpy(),
        f"human_{method}": data["Homo sapiens"].astype(str).str.upper().to_numpy(),
    })


def load_variants(old_csv):
    raw = pd.read_csv(old_csv, low_memory=False)
    required = {"POS", "REF", "ALT", "clinvar_interp_new", "mitomap_status_new"}
    if not required.issubset(raw):
        raise ValueError(f"{old_csv}: missing {sorted(required-set(raw))}")
    clin = raw.clinvar_interp_new.fillna("").astype(str).str.strip()
    mito = raw.mitomap_status_new.fillna("").astype(str)
    positive = clin.isin(POSITIVE) | mito.str.contains("Cfrm", regex=False)
    negative = clin.eq("Benign") & ~positive
    labeled = raw.loc[positive | negative, ["POS", "REF", "ALT"]].copy()
    labeled["label"] = positive.loc[labeled.index].astype(int)
    labeled["POS"] = pd.to_numeric(labeled.POS, errors="raise").astype(int)
    labeled["REF"] = labeled.REF.astype(str).str.upper()
    labeled["ALT"] = labeled.ALT.astype(str).str.upper()
    labeled = labeled.drop_duplicates(["POS", "REF", "ALT", "label"])
    # Exclude individual alleles with conflicting labels, not whole positions.
    ambiguous = labeled.groupby(["POS", "REF", "ALT"])["label"].transform("nunique") > 1
    n_conflicts = int(ambiguous.sum())
    if n_conflicts:
        print(f"[WARN] Excluding {n_conflicts} conflicting variant-label rows")
    return labeled.loc[~ambiguous].copy()


def load_rank(root, rank, metric, variants):
    folder = root / rank
    paths = {
        "original": folder / "conservation" / "conservation_ci.csv",
        **{m: folder / "overlap_aware_ci" / m / "conservation_ci.csv"
           for m in METHODS if m != "original"},
    }
    tables = {method: score_table(path, metric, method)
              for method, path in paths.items()}
    # Trust annotation classification from the same CI computation.
    qcpath = folder / "overlap_aware_ci" / "overlap_position_ci.csv"
    qc = pd.read_csv(qcpath)
    if (len(qc) != 16569 or
            qc.Position.tolist() != list(range(1, 16570)) or
            not {"Position", "is_overlap"}.issubset(qc)):
        raise ValueError(f"Invalid overlap position metadata: {qcpath}")
    joined = variants.copy()
    for method in METHODS:
        joined = joined.merge(tables[method], on="POS", how="left",
                              validate="many_to_one")
    joined = joined.merge(
        qc[["Position", "is_overlap"]].rename(columns={"Position": "POS"}),
        on="POS", how="left", validate="many_to_one"
    )
    for method in METHODS:
        if not (tables[method][f"human_{method}"].values ==
                tables["original"]["human_original"].values).all():
            raise ValueError(f"{rank}: mismatched human reference for {method}")
    human = tables["original"]["human_original"]
    if not joined.REF.eq(human.iloc[joined.POS.to_numpy()-1].to_numpy()).all():
        print(f"[WARN] {rank}: some labeled REF alleles differ from human reference; excluding them")
        joined = joined.loc[
            joined.REF.eq(human.iloc[joined.POS.to_numpy()-1].to_numpy())
        ].copy()
    # Exact score preservation outside gene-gene overlap, as designed.
    not_overlap = qc.is_overlap.eq(0).to_numpy()
    original = tables["original"]["original"].to_numpy()
    for m in ("mean", "min", "max"):
        curr = tables[m][m].to_numpy()
        if not np.all(np.isclose(curr[not_overlap], original[not_overlap],
                                 rtol=0, atol=1e-8, equal_nan=True)):
            raise ValueError(f"{rank}: non-overlap scores changed for {m}")
    joined = joined.dropna(subset=list(METHODS) + ["is_overlap"]).reset_index(drop=True)
    joined["is_overlap"] = joined["is_overlap"].astype(int)
    joined["score_changed"] = np.any(np.column_stack([
        ~np.isclose(joined[m], joined["original"], rtol=0, atol=1e-8)
        for m in ("mean", "min", "max")
    ]), axis=1)
    joined["rank"] = rank
    return joined


def evaluate(d, rank, group, reps, seed, min_per_class):
    y = d.label.to_numpy(dtype=np.int8)
    n_pos = int(np.sum(y == 1))
    n_neg = int(np.sum(y == 0))
    base = dict(rank=rank, group=group, n_variants=len(d),
                n_positions=d.POS.nunique(),
                n_pathogenic=n_pos, n_benign=n_neg)
    if min(n_pos, n_neg) < min_per_class or len(d) == 0:
        return {**base, "status": "insufficient_labels"}, [], {}
    scores = {m:d[m].to_numpy(dtype=float) for m in METHODS}
    points = {m:float(roc_auc_score(y, scores[m])) for m in METHODS}
    clusters, inverse = np.unique(d.POS.to_numpy(), return_inverse=True)
    rng = np.random.default_rng(seed)
    draws = []
    for rep in range(reps):
        sample = rng.integers(len(clusters), size=len(clusters))
        weights = np.bincount(sample, minlength=len(clusters))[inverse]
        if not weights[y == 0].any() or not weights[y == 1].any():
            continue
        values = {m:float(roc_auc_score(y, scores[m], sample_weight=weights))
                  for m in METHODS}
        draws.append({"rank":rank, "group":group, "replicate":rep,
                      **{f"auc_{m}":values[m] for m in METHODS}})
    row = {**base, "status": "ok" if draws else "no_valid_bootstraps",
           **{f"auc_{m}":points[m] for m in METHODS},
           "valid_bootstraps":len(draws)}
    for m, baseline in COMPARISONS:
        name = f"{m}_minus_{baseline}"
        diff = np.array([r[f"auc_{m}"] - r[f"auc_{baseline}"] for r in draws])
        row[f"delta_{name}"] = points[m] - points[baseline]
        row[f"ci_lower_{name}"] = np.quantile(diff, 0.025) if len(diff) else np.nan
        row[f"ci_upper_{name}"] = np.quantile(diff, 0.975) if len(diff) else np.nan
    return row, draws, scores


def plot_rank(d, rank, metric, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    for ax, (name, subset) in zip(
        axs, (("genome_wide", d), ("gene_gene_overlap", d[d.is_overlap == 1]))
    ):
        if subset.label.nunique() < 2 or len(subset) < 2:
            ax.text(.5, .5, "Insufficient labels", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title(name)
            continue
        y = subset.label.to_numpy(dtype=int)
        for m in METHODS:
            fpr, tpr, _ = roc_curve(y, subset[m])
            auc = roc_auc_score(y, subset[m])
            ax.plot(fpr, tpr, label=f"{m}: {auc:.3f}")
        ax.plot([0, 1], [0, 1], "k--", alpha=.3)
        ax.set_title(f"{name} (P={int(sum(y))}, B={int(sum(1-y))})")
        ax.legend(fontsize=8)
        ax.set(xlabel="False positive rate", ylabel="True positive rate")
    fig.suptitle(f"{rank}: pad100 Original vs Mean / Min / Max ({metric}, raw CI)")
    fig.tight_layout()
    fig.savefig(out / f"{rank}_roc.png", dpi=170)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-root",
                   default="block_results_genbank_trnascan_mitos_fallback_v2/pad100")
    p.add_argument("--old-csv", default="data/2026_10_10_mito_cons_anns_main.csv")
    p.add_argument("--metric", choices=("CI_MSA", "CI_callable"), default="CI_MSA")
    p.add_argument("--ranks", nargs="+", choices=RANKS, default=list(RANKS))
    p.add_argument("--bootstrap-reps", type=int, default=1000)
    p.add_argument("--min-per-class", type=int, default=5)
    p.add_argument("--seed", type=int, default=20261119)
    p.add_argument("--out-dir")
    a = p.parse_args()
    if a.bootstrap_reps < 1 or a.min_per_class < 1:
        p.error("bootstrap-reps and min-per-class must be positive")
    root = Path(a.results_root)
    out = Path(a.out_dir or root / OUTPUT_NAME / a.metric)
    out.mkdir(parents=True, exist_ok=True)
    variants = load_variants(a.old_csv)
    summary, bootstrap, cohorts = [], [], []
    for rank in a.ranks:
        d = load_rank(root, rank, a.metric, variants)
        cohorts.append(d)
        groups = {
            "genome_wide":d,
            "gene_gene_overlap":d.loc[d.is_overlap == 1],
            "changed_overlap":d.loc[(d.is_overlap == 1) & d.score_changed],
        }
        for idx, (name, subset) in enumerate(groups.items()):
            row, draws, _ = evaluate(subset,rank,name,a.bootstrap_reps,
                                     a.seed + 11*a.ranks.index(rank) + idx,
                                     a.min_per_class)
            summary.append(row)
            bootstrap.extend(draws)
        plot_rank(d, rank, a.metric, out)
        print(f"{rank}: variants={len(d)}, overlap_variants={len(groups['gene_gene_overlap'])}, changed_overlap_variants={len(groups['changed_overlap'])}")
    pd.DataFrame(summary).to_csv(out/"auc_summary.csv",index=False)
    pd.DataFrame(bootstrap).to_csv(out/"auc_bootstrap.csv",index=False)
    pd.concat(cohorts,ignore_index=True).to_csv(out/"labeled_variant_scores.csv",index=False)
    (out/"metadata.json").write_text(json.dumps({
        "results_root":str(root),"label_source":a.old_csv,"metric":a.metric,
        "ranks":a.ranks,"methods":METHODS,"positive_labels":sorted(POSITIVE),
        "mitomap_positive":"contains Cfrm","negative_label":"ClinVar Benign and not positive",
        "score_orientation":"raw conservation CI; no inversion",
        "bootstrap":"paired human position-cluster resampling per rank and group",
        "min_per_class":a.min_per_class,"bootstrap_reps":a.bootstrap_reps,
        "note":"Overlap-only AUROC may be undefined; genome-wide difference is diluted when only a few sites change."
    }, indent=2)+"\n")
    print(pd.DataFrame(summary)[["rank","group","n_pathogenic","n_benign","status"]].to_string(index=False))
    print("Outputs:",out)


if __name__ == "__main__":
    main()
