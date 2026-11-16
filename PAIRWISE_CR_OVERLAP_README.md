# Pairwise human-CR versus longest gene-free interval (all species)

Run with `sbatch run_compare_pairwise_cr_to_longest_intervals.slurm`.

The script **retains every species** in each requested taxonomic rank and
calculates three overlap checks, including two checks against the species' longest unfiltered
gene-free interval, all in unrotated mitochondrial coordinates:

1. **CR overlap**: Human NC_012920.1 CR positions 16024..16569 and 1..576
   (1,122 human bases) are selected by counting human *non-gap* nucleotides
   in each full-length pairwise alignment. The species base in each such
   column is mapped to the original circular genome after identifying a
   unique exact rotation of the ungapped aligned species sequence. Species
   insertions in human-gap columns are not counted. An overlap requires >=1
   mapped species base within the longest interval. This is **not evidence
   of homology by itself** for distant taxa. The overlap fraction
   `mapped_CR_coverage_by_gap` is the fraction of mapped species CR
   nucleotides within the longest interval (not the human CR mapping fraction).
2. **OriH overlap**: taken independently from the MITOS2 OriH columns in
   `longest_gene_free_interval_by_species.tsv`. A missing prediction or
   GFF is marked `unavailable`, **not** `no overlap`.

3. **Direct pairwise CR versus MITOS2 OriH overlap**: read the original
   per-species MITOS2 GFF from `data/mitos2_species_summary.tsv`.
   Only parent `origin_of_replication` / `rep_origin` features labeled
   `OH` / `OriH` are used, matching
   `analyze_longest_all_gene_desert.py`; child OH-b exons and OriL
   do not contribute. All OriH feature bases are unioned, respecting
   circular wrap, and intersected directly with the individual mapped
   CR species positions. This comparison does **not** require an
   existing longest gene-free interval. `CR_OriH_overlap_bp` and
   `CR_OriH_overlap_any` report exact base-pair intersection and a
   1/0 flag. Unavailable means missing OriH annotation, invalid CR
   alignment, or no mapped CR bases — not a biological negative.

4. **OriH-guided candidate selection** (both orders):
   - Longest-first: select the longest interval if it overlaps OriH;
     otherwise choose the pairwise-mapped CR if *it* directly overlaps OriH.
   - CR-first: select mapped CR if it directly overlaps OriH; otherwise
     choose the longest interval if it overlaps OriH.

   Since both use the same logical OR of two criteria, **the total number
   selected is the same**. However, the primary/fallback assignments and
   the genomic regions selected can differ. Both-positive and one-positive
   species are counted separately. All species stay in the denominator,
   and the CR coverage threshold is QC only. A mapped CR is a set of
   aligned species bases; it does not itself define a single contiguous
   interval with inserted unaligned nucleotides.

The --min-cr-mapped-fraction cutoff (default 0.50) **does not filter species
or overlap calculations**. Validated alignments below cutoff are labeled
`CR_mapping_group=below_threshold`; those above are
`at_or_above_threshold`. Invalid/ambiguous alignments are
`no_valid_mapping`. The single exact-rotation ambiguity is retained in
species output with CR overlap unavailable. Species without a definable
longest interval also remain in the output; overlap is unavailable. A
validated pairwise alignment with zero mapped human CR bases is retained,
but the overlap is unavailable, not a negative result.

An interval with `gap_end_1based=0` ends at the last genome base; the
comparison handles this legacy circular sentinel.

## Outputs

All files are produced in
`longest_all_gene_desert_qc/pairwise_CR_validation/`:

- `all_ranks_species_pairwise_CR_vs_longest.tsv`: **every** rank/species,
  CR mapping percentage and quality tier, CR overlap columns,
  independent OriH prediction/overlap columns, direct CR/OriH overlap,
  selected source for each priority order, and joint category.
  `CR_overlap_any` and `OriH_overlaps_longest`: 1=yes, 0=no,
  empty=unavailable.
- `all_ranks_summary_pairwise_CR_vs_longest.tsv`: denominators, counts
  of CR and OriH overlaps **separately**, below-threshold CR hits,
  and four-way intersection where both can be evaluated.
- `all_ranks_CR_overlap_by_coverage.tsv`: CR overlap by rank/coverage tier.
- `all_ranks_OriH_overlap.tsv`: OriH status and longest-interval overlap
  by rank.
- `all_ranks_CR_OriH_joint_by_coverage.tsv`: joint CR/OriH outcomes
  stratified by mapping coverage, including unavailable values.
- `all_ranks_support_categories.tsv`: joint categories with all species.
- `all_ranks_OriH_guided_fallback_strategies.tsv`: per-rank first-class,
  fallback, total selection and overlap counts for both priority orders.
- `pairwise_mapping_QC.tsv`: rotation validation details and counts.
- `duplicate_rank_species_names.tsv`: any ambiguous labels.

**Denominators:** `n_total_species` includes all rank species;
`n_CR_overlap_evaluable` requires a valid mapping with >=1 mapped CR base
and a longest interval; `n_OriH_overlap_evaluable` requires a predicted
OriH, available GFF, and longest interval. Do not interpret unavailable as
no overlap. The historical `n_mapped_CR_usable` column continues to
represent CR coverage >=50% for comparison with earlier results.

The new GFF input defaults to `data/mitos2_species_summary.tsv`.
Use `--mitos2-summary /full/path/summary.tsv` if stored elsewhere.
The SLURM wrapper supports the `MITOS2_SUMMARY` environment variable.

Tests: `python -m unittest test_compare_pairwise_cr_to_longest_intervals.py`.
