# Longest gene-free interval vs GenBank control region — 100% OriH containment

Replaces **pair-wise mapped CR** with **current GenBank D-loop + misc_feature control-region annotations** in the comparative plot. It does not remove or overwrite the original pair-wise analysis.

Input data:
- `longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv` (previously computed longest interval, 10 ranks)
- `data/kingdom_genbank_annotations_split_new.csv` (D-loop, including multipart segments)
- `data/genbank_misc_ftr.csv` (control-region misc_features; excludes subfeatures/remnants/AT-rich-only candidates from strict CR)
- `data/mitos2_species_summary.tsv` (MITOS2 GFF parent OH/OriH annotations)

## Scientific definition

The denominator is the distinct union of parent OriH bases on an mtDNA sequence. A method succeeds only if **a single physical candidate region** contains **every** OriH base. GenBank CR1 and CR2 remain separate; combining CR1+CR2 to reach 100% is NOT permitted. The longest gap and GenBank CR also cannot be combined to claim individual-method containment.

- **Longest**: one longest gene-free interval contains all OriH bases.
- **GenBank CR**: at least one individually annotated full or putative control region contains all OriH bases.
- **Both**: each of those methods succeeds independently.
- **Either**: at least one succeeds independently.

Two denominators per rank: (1) MITOS2 OriH-predicted species and (2) all species in the longest-gap rank TSV. A species with no OriH prediction counts as a failure in all-species percentage, *not* as biological evidence of noncontainment; its status is explicit in the species table.

**Caution:** Source accession, genome orientation and circular rotation must match across longest-gap TSV, GenBank coordinates, and MITOS2 GFF. The current comparison is a provisional physical-coordinate comparison unless provenance is independently verified. Multiple GenBank accessions for one species are flagged and not scored as a GenBank success. The longest-gene-desert input may have been built from older GenBank annotations; the comparison uses that historical longest result for direct continuity, without silently recomputing it.

## McCleary

```bash
cd /vast/palmer/scratch/lake_nicole/lt692/chrM_block_align
git pull origin main
mkdir -p logs
python3 test_compare_genbank_cr_to_longest_orih.py
sbatch run_compare_genbank_cr_to_longest_orih.slurm
```

Or directly:

```bash
python3 compare_genbank_cr_to_longest_orih.py \
  --longest-tsv longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv \
  --dloop data/kingdom_genbank_annotations_split_new.csv \
  --misc data/genbank_misc_ftr.csv \
  --mitos-summary data/mitos2_species_summary.tsv \
  --out-dir longest_all_gene_desert_qc/genbank_CR_comparison
```

## Outputs

- `species_100pct_containment.tsv`: one row per rank/species, OriH status/bp, each success flag, matching GenBank CR labels and intervals, coordinate caution, and raw overlaps.
- `rank_100pct_containment.tsv`: rank/denominator/sample sizes, counts and success percentages for Longest, GenBank CR, Both, Either.
- `orih_containment_genbank_predicted.png` and `.pdf`: bar chart among OriH-predicted species.
- `orih_containment_genbank_all_species.png` and `.pdf`: bar chart among all species.
- `summary.json`: method and quality caveats.

The code reuses the existing **parent OriH GFF filters**, legacy circular gene-desert interval conventions, and ten-rank species cohorts from `compare_pairwise_cr_to_longest_intervals.py`.

These plots measure recovery of MITOS2 OriH predictions; they do **not** establish which region is an active replication origin.
