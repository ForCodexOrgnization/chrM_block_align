# chrM block alignment (v6)

Human-anchored mitochondrial alignment across species with rearranged gene order.

## Design

- MITOS2 GFF `gene` (CDS) and `rRNA` annotations; tRNAscan CSV for tRNA.
- Independently align **complete** homologous genes, even if their genomic spans overlap.
- Intergenic segments between non-overlapping reference genes are handled separately, and only when flanking gene adjacency is supported.
- Each reference nucleotide is emitted **once** in the final stitched alignment. Overlap duplicates are reported.
- The human reference is an actual 16,569-bp sequence, not padded coordinates.
- Only aligned reference columns are emitted into stitched FASTA; insertion columns are retained in intermediate gene MSAs. The stitched output is a human-coordinate projection, **not** a physical genome alignment.

## Requirements

Python 3.9+, MAFFT 7.526. On McCleary:

```bash
module load MAFFT/7.526-GCC-13.3.0-with-extensions
mkdir -p logs
sbatch run_class.slurm
```

Expected input layout:

```text
data/
  unrot_fasta/class_chrm.faa
  mitos2_species_summary.tsv
  consolidated_trnascan-2.csv
  human_chrM_16569.fa
```

For direct execution:

```bash
python3 chrM_block_align.py --fasta-dir data/unrot_fasta --annotation-summary data/mitos2_species_summary.tsv --trnascan-csv data/consolidated_trnascan-2.csv --human-fasta data/human_chrM_16569.fa --output-dir block_results_v6 --ranks class --threads 4 --keep-intermediate
python3 validate_alignment.py --rank-dir block_results_v6/class --human-fasta data/human_chrM_16569.fa
```

## Important caveats

- A missing or ambiguous ortholog is represented as a gap (unavailable) rather than treated as a substitution.
- Intergenic blocks with rearranged neighbors are not assumed homologous. Candidate flanks are added in MAFFT batches (default 100); unusually long non-human flanks (>2500 bp by default) are excluded individually and counted in `intergenic_qc.tsv`. Existing gene alignments are reused when the cached sample IDs and ungapped sequences match exactly.
- Separate gene alignments preserve the complete overlap, but the final chrM projection chooses one aligned result per human nucleotide. It records overlap conflicts and does not automatically validate orthology in ambiguous duplications.
- The Homo sapiens sequence in every rank FASTA, the --human-fasta sequence, MITOS2 GFF, and tRNAscan coordinates must all refer to the same exact 16,569-bp human reference. No liftover from 16,565 bp is performed. The script checks FASTA identity and GFF region length; tRNAscan coordinates must be independently verified.
- Review generated QC and validate on real McCleary data before computing CI. This is research code, not a clinically validated pipeline.

## Re-running after v6 intergenic update

Run the same `sbatch run_class.slurm` command with `--keep-intermediate` enabled (already present in the sample SLURM script). Complete cached `*.aln.fa` files are validated and reused; intergenic blocks previously saved with human-only results are regenerated when their expected sequence membership differs. Check `intergenic_qc.tsv` and `block_summary.tsv`. The intended human reference file for the current dataset is `data/Humans_mtdna.fasta`.

The former `--max-flank-sequences` threshold is no longer used to skip intergenic MSAs; `--flank-batch-size` (default 100) controls progressive addition.

## Canonical human annotation (required)

Install Biopython (`python3 -m pip install biopython`), copy the matched `NC_012920.1` GenBank file to `data/sequence.gb`, and use `data/Humans_mtdna.fasta` as the 16569-bp reference. Set `--human-genbank data/sequence.gb`. All **37 human features** (13 CDS, 2 rRNAs, 22 tRNAs) are parsed directly from this GenBank; the old human MITOS2 GFF and human tRNAscan records are ignored. Other species still use MITOS2 GFF and tRNAscan. The GenBank nucleotide sequence must match the human FASTA exactly. If changing gene boundaries, any affected cached alignments are validated and recalculated based on input sequences.

## Conservation scores from stitched MSA

The SLURM workflow now runs `calculate_conservation.py` after `validate_alignment.py`. To run **without rerunning MAFFT**, execute:

```bash
python3 calculate_conservation.py \
  --msa block_results_v6/class/class_stitched_aln.fasta \
  --coordinate-map block_results_v6/class/human_coordinate_map.tsv \
  --out-dir block_results_v6/class/conservation
python3 test_conservation.py
```

Outputs:
- `conservation/conservation_ci.csv`: one row for each human position with reference nucleotide, `CI_MSA` (percent matches across **all** species including Homo sapiens, matching the supplied legacy algorithm), `CI_callable` (percent matches among A/C/G/T only), and denominators/missingness counts.
- `conservation/nucleotide_frequencies.csv`: a/c/g/t/- percentages over all valid A/C/G/T/gap characters, ignoring N, together with counts.
- `conservation/conservation_summary.json`: definition and coverage summary.

**Interpretation:** Gaps in the stitched MSA also represent missing annotations or unaligned blocks, not necessarily biological deletions. `CI_MSA` counts them in its denominator, as the legacy code did. For biologically interpretable scores, inspect `n_callable` and `CI_callable`. This tool does not generate the enormous per-species MSA CSV by default, because the stitched FASTA already contains those sequences.

## Compare old chrM CI with new stitched MSA CI

The original annotations file (`2026_10_10_mito_cons_anns_main.csv`) has `POS` and ten old `{rank}_rot_ref` scores (0–100). It contains multiple REF/ALT variants at many positions. The new pipeline writes `block_results_v6/{rank}/conservation/conservation_ci.csv`, one row per human position, with `CI_MSA` and `CI_callable` percentages (0–100). Both old and new scores are preserved on the 0–100 scale in tables and plots.

Place the old CSV under `data/2026_10_10_mito_cons_anns_main.csv` and run:

```bash
python3 compare_old_new_ci.py \
  --old-csv data/2026_10_10_mito_cons_anns_main.csv \
  --results-dir block_results_v6 \
  --new-metric CI_MSA \
  --out-dir block_results_v6/ci_comparison
```

This script skips unfinished ranks (unless `--strict-all-ranks` is set), retains all variant rows for downstream ClinVar ROC, and makes a distinct-position table for comparisons. Output: `variant_old_new_ci.csv`, `position_old_new_ci.csv`, `comparison_summary.csv`, `paired_roc.png`, `score_delta_by_rank.png`, `gene_heatmap_old.png`, `gene_heatmap_new.png`, and `roc_auc_summary.csv`. Uses pandas, numpy, matplotlib, scikit-learn (install if missing). `--no-plots` only requires pandas/numpy.

**Methodological caution:** old `{rank}_rot_ref` values may include zeros in the supplied annotation CSV. Confirm whether zero encodes genuinely zero conservation or unavailable/not computed scores; the script flags this instead of changing values. ROC uses the original plotting script's raw-score orientation, so an AUC below 0.5 may mean conservation and pathogenicity are inversely correlated, not that the predictor has no signal.


### CI diagnostic analyses: boundaries, coverage, and paired ROC

Once `compare_old_new_ci.py` has created `block_results_v6/ci_comparison/position_old_new_ci.csv` and `variant_old_new_ci.csv`, run:

```bash
mkdir -p logs
sbatch run_ci_diagnostics.slurm
# Or interactively with a suitable Python environment:
python3 analyze_ci_diagnostics.py --comparison-dir block_results_v6/ci_comparison \
  --human-genbank data/sequence.gb --bootstrap-reps 1000
```

The diagnostic script skips missing ranks (e.g., kingdom while running), reads canonical human CDS/rRNA/tRNA boundaries from `data/sequence.gb` including circular wrap-around, and computes per-position absolute old/new CI differences by distance-to-boundary bins and callability (`n_callable/n_species`, default <50%, 50%-90%, >=90%). Outputs in `ci_comparison/diagnostics/`: `boundary_summary.csv`, `boundary_abs_ci_difference.png`, `coverage_summary.csv`, `coverage_abs_ci_difference.png`, `position_diagnostics.csv`, `paired_roc_bootstrap_summary.csv`, and `paired_roc_bootstrap_replicates.csv`.

Paired ROC AUC bootstrap resamples **human positions as clusters**, taking all associated REF/ALT variants together, and computes **new minus old AUC** on identical labeled variants per rank; default 1000 replicates, fixed seed. The 95% interval is a bootstrap percentile interval and is not a formal significance test. CI scores remain in the original **0–100** scale (AUC itself is 0–1). All required Python dependencies: `biopython pandas numpy matplotlib scikit-learn`. The SLURM job uses the configured Hail environment; change `PYTHON` there if these libraries are not installed.

### Preliminary intergenic-length survey (before padding)

Run the small **class-only** annotation survey without MAFFT:

```bash
git pull origin main
mkdir -p logs
sbatch run_intergenic_survey.slurm
```

Or directly: `python3 survey_intergenic_lengths.py --ranks class`.
Outputs: `block_results_v6/intergenic_length_survey/class_intergenic_summary.tsv`, `class_intergenic_observations.tsv`, `class_skipped_species.tsv`, and `all_ranks_intergenic_summary.tsv`. All intervals use local species circular chrM coordinates and human GenBank gene-pair order. `gap_le_50_bp_fraction` means the entire gap fits within **one** 50-bp flank; `gap_le_100_bp_fraction` fits one 100-bp flank (or two 50-bp flanks); `gap_le_200_bp_fraction` fits two 100-bp flanks. Statistics exclude negative overlaps and absent/rearranged/nonadjacent pairs, which are reported separately. The survey does not establish homology of flanks and should not alone determine padding length.

### Experimental padded gene-block alignment (baseline untouched)

**Safety:** This is a new, experimental implementation. It does not modify `chrM_block_align.py`, `run_all_ranks.slurm`, or `block_results_v6/`. The ongoing kingdom baseline job can continue. Test `class` first.

```bash
git pull origin main
mkdir -p logs
python3 test_padded_alignment.py
sbatch run_padded_class.slurm                # +/-50 bp; default
# For +/-100 bp, submit separately after the first test:
PAD=100 sbatch --export=ALL,PAD=100 run_padded_class.slurm
```

`chrM_block_align_padded.py` extracts full annotated genes plus +/-50 or +/-100 bp flanks in each species, including circular wrap-around. A nonhuman flank is only extracted if the expected neighboring gene and strand are supported. Padded MSAs use the canonical GenBank human reference to project onto 16,569 positions; gene core candidates take priority over padding, then closest anchor wins with deterministic gene order ties. The independent TRNP-TRNF intergenic MSA is retained and used as a fallback; ordinary IGR alignment is performed only when padded blocks leave reference columns uncovered. Note that this is a first experimental implementation, not yet biologically validated: strict adjacency gating may discard useful flanks, and the position-based owner rule should be reviewed against alternative region-level stitching before making biological claims.

Output: `block_results_padded_v1/pad50/class/` or `pad100/class/`, including stitched FASTA, coordinate map, `stitch_decisions.tsv`, `overlap_conflicts.tsv`, `flank_context_qc.tsv`, `intergenic_qc.tsv`, conservation scores and validation report. Original CI scoring remains exact human-base matches divided by **all species**, including gaps/IUPAC ambiguity as nonmatches.
