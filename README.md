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

### Three-way CI ROC: old versus gene-only versus +/-50 bp padded

Run `mkdir -p logs && sbatch run_three_way_roc.slurm` after the class padded conservation CSV is ready. The new script `compare_three_way_roc.py` uses **identical clinically labeled REF/ALT variants** for all three 0-100 CI measures: `class_rot_ref` (old), `block_results_v6/class/conservation/conservation_ci.csv` (gene-only), `block_results_padded_v1/pad50/class/conservation/conservation_ci.csv` (padded). ClinVar P/LP or MITOMAP `Cfrm` are positive; ClinVar Benign without `Cfrm` are negative. Duplicate variant records are removed. All AUC differences use paired **position-clustered bootstrap** (default 1000), 95% percentile intervals. Groups: genome-wide, protein coding (CDS), tRNA, rRNA, intergenic, human GenBank gene-boundary windows **0-5, 6-25, 26-100, >100 bp**, and **within-gene** AUC comparing pathogenic/benign variants *only within the same human gene* (weighted by the available within-gene positive-negative pairs). Groups with <5 examples of either label are marked `insufficient_labels` rather than assigned a misleading AUC. Outputs: `three_way_auc_summary.csv`, `labeled_variant_scores.csv`, `three_way_genome_roc.png`.

**Interpretation:** gene-category composition and differences in label availability can bias global comparisons. Boundary-specific estimates may have few labels; report the counts and confidence intervals. Within-gene AUC does not remove all ascertainment bias. AUC tests pathogenicity ranking, not alignment orthology correctness. Existing baseline alignment scripts and output directories remain untouched.

### All-rank genome-wide ROC after kingdom baseline completion

The standalone `compare_all_ranks_genome_roc.py` script compares **original `{rank}_rot_ref`** with **block v6 `CI_MSA`** for all ten taxonomic ranks, including kingdom. It requires `block_results_v6/{rank}/conservation/conservation_ci.csv` for each rank. All ranks are scored on the **same complete-case set of uniquely labeled POS/REF/ALT variants**, with paired 1000-replicate position-cluster bootstrap CIs. It checks CI in 0–100 and retains the previous ClinVar/MITOMAP label definitions.

```bash
git pull origin main
mkdir -p logs
sbatch run_all_ranks_genome_roc.slurm
```

Files are created in `block_results_v6/all_ranks_genome_roc/`: `all_ranks_roc_summary.csv`, `all_ranks_roc_bootstrap.csv`, `all_ranks_roc.png`, `shared_labeled_variants.csv`. This only reads the completed alignment/conservation outputs; it does not rerun alignment or edit original files. Note: positive/negative clinical label composition can differ across gene categories. This genome-wide test is an overall benchmark and does not replace within-gene or boundary-specific ROC.

### Run both +/-50 bp and +/-100 bp padding across all ten ranks

One submission creates a **20-task SLURM array** (10 ranks x 2 padding choices), with at most two tasks running concurrently. Each task runs the alignment, structural validation, and conservation calculation separately. The same script automatically maps array indices 0-9 to pad50 and 10-19 to pad100. No PAD environment variable is necessary.

```bash
git pull origin main
mkdir -p logs
sbatch run_padded_all_ranks.slurm
squeue -u "$USER"
```

Outputs are isolated in `block_results_padded_v1/pad50/<rank>/` and `block_results_padded_v1/pad100/<rank>/`. Existing validated output is skipped. Partially completed output causes that task to stop without overwriting it, so review any error log before resubmission. Old `block_results_v6/` results are unchanged. Each concurrent task requests 4 CPUs, 96 GiB of RAM, and up to 3 days; confirm these are allowed by your HPC partition. Higher ranks may take substantially longer, and the padded strategy is still experimental.


### Four-method all-rank stratified diagnostics

The **read-only** `analyze_four_way_diagnostics.py` compares original CI, v6 gene-only, pad50 and pad100, using shared ClinVar/MITOMAP-labeled variants at every rank. It produces genome-wide, gene-type (CDS/tRNA/rRNA/intergenic), nearest-human-gene-boundary (0–5, 6–25, 26–100, >100 bp), and **within-gene** pathogenicity concordance AUC (only positive/negative variant pairs from the same gene; weighted by eligible pair counts). All 6 pairwise method effects use a **paired position-cluster bootstrap**, with 95% percentile CIs. Sparse groups (<5 positive or <5 negative) are marked insufficient rather than assigned AUC values.

It also plots position-level mean absolute CI changes (gene boundary, other gene positions, intergenic) and padded-flank eligibility from each rank's `flank_context_qc.tsv`. The eligibility denominator is *annotated gene/species pairs*, not all species and not aligned bases. Position categories use human GenBank annotations; overlapping genes are assigned to tRNA before rRNA before CDS for grouping. These diagnostics are exploratory, not independent verification of orthology.

```bash
git pull origin main
mkdir -p logs
sbatch run_four_way_diagnostics.slurm
```

Outputs: `block_results_padded_v1/four_way_diagnostics/` with `stratified_four_way_auc.csv`, `stratified_four_way_bootstrap.csv`, `position_ci_change_by_region.csv`, `flank_eligibility_by_rank.csv`, and seven+ PNG figures. The initial four-way genome-wide ROC from `compare_all_ranks_padding_roc.py` remains separate.

### Experimental GenBank + MITOS2 block pipeline (separate from current methods)

This new **independent** pipeline uses GenBank gene-coordinate annotations where available, MITOS2 gene/rRNA coordinates as fallback, and does **not use tRNAscan**. The reference human always uses the canonical 16,569-bp human GenBank (`data/sequence.gb`). All gene names are mapped through the uploaded synonym table `genefiltfinal5_2026_09_24.csv` (canonical names are columns, alternate names are cells). **Two distinct CSV files are required**:

1. `data/genefiltfinal5_2026_09_24.csv`: the supplied 92-row gene synonym table (copy/upload to McCleary; this file is not automatically committed).
2. `data/kingdom_genes_v2.csv`: export the *actual* GenBank feature/coordinate worksheet shown in the screenshot as CSV. It must contain `Organism`, `Gene`, `start`, `end`, `strand` and ideally `acc`. The synonym CSV by itself has **no positions** and cannot replace this file. Coordinates are assumed 1-based inclusive and in the same circular sequence orientation/origin as each `data/unrot_fasta/{rank}_chrm.faa`; if not, perform sequence-origin mapping first. Organism values must match rank FASTA organism names after punctuation normalization. Multiple accessions for one organism require manual resolution before use, otherwise ambiguous loci are excluded.

```bash
git pull origin main
mkdir -p data logs
# copy BOTH CSVs to data/ first
sbatch run_genbank_mitos_all_ranks.slurm
```

This submits ten rank tasks, each sequentially running gene-only, pad50, pad100. Outputs: `block_results_genbank_mitos_v1/{gene_only,pad50,pad100}/{rank}` with stitched alignment, coordinate map, method-specific QC; standard validation and conservation summaries are produced when all three methods complete. Annotation provenance: `block_results_genbank_mitos_v1/annotation_qc/{rank}/annotation_source_audit.tsv`, `ambiguous_gene_annotations.tsv`, `unmapped_gene_names.tsv`, `invalid_rows.tsv`. Run these as **experimental** results; compare their gene coverage and coordinate consistency against the original MITOS2+tRNAscan pipeline before scientific interpretation.

**Safety:** The original `chrM_block_align.py`, `chrM_block_align_padded.py`, all prior output folders, and their SLURM launchers are unchanged. This script calls the original methods with new in-memory annotation features, without modifying them. Existing output directories are not overwritten. No species/genome rotations are inferred by the script; incorrect GenBank-to-FASTA coordinate origin would invalidate alignment even if structural checks pass.

### Gene-order-assisted annotation resolution (new GenBank + MITOS2 pipeline only)

`run_genbank_mitos_blocks.py` now imports `resolve_gene_order_annotations.py` before extracting features. It recognizes synonyms from **Gene, G_Symbol, G_locus_tag, G_Sequence, G_tag, alt_name, alt_synonym, alt_Product**, keeping candidate gene identities when the fields conflict or tRNA Leu/Ser subtype is ambiguous. Candidate assignment tests the **complete circular gene order with strand** against the canonical human GenBank annotation (rotations and whole-genome reverse orientation are permitted for **annotation comparison only**). A gene name is resolved automatically only when a **unique exact-order solution** exists. Inexact best-fit assignments, ties, and >4096 candidate combinations remain unresolved and are recorded; they are **not** forced to human order. Genuine rearrangements are therefore not automatically rewritten.

Per-rank QC: `block_results_genbank_mitos_v1/annotation_qc/<rank>/gene_order_resolutions.tsv`, alongside provenance, ambiguous locus, invalid row, and unmapped-name tables. Genes unresolved by GenBank name matching may still use unambiguous MITOS2 protein/rRNA fallback through the normal annotation source precedence. No source coordinates or FASTAs are rotated by gene-order matching, and the GenBank annotation CSV must still refer to the same actual genomic coordinate origin as the FASTA. The main `chrM_block_align.py` and `chrM_block_align_padded.py` remain unchanged.

Smoke test: `python test_gene_order_resolution.py`. This is *conservative, annotation-order-informed resolution*, not independent confirmation of homology. Check `gene_order_resolutions.tsv` and coordinate matching before full-scale deployment. Existing GenBank+MITOS2 experimental output is **not** overwritten on reruns; use a fresh output directory for a new version.
