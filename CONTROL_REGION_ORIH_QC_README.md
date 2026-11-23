# GenBank D-loop / control region vs MITOS2 OriH

The analysis merges:
- `data/kingdom_genbank_annotations_split_new.csv`: features with `class=D-loop`.
- `data/genbank_misc_ftr.csv`: `misc_kind=other` and notes/flags identifying control-region/D-loop features. CDD structural annotations are excluded.
- `data/mitos2_species_summary.tsv`: species-to-`result_gff` paths for MITOS2 predicted OriH features.

## Run on McCleary

```bash
cd /vast/palmer/scratch/lake_nicole/lt692/chrM_block_align
git pull origin main
test -s data/kingdom_genbank_annotations_split_new.csv
test -s data/genbank_misc_ftr.csv
test -s data/mitos2_species_summary.tsv
python3 test_control_region_orih_qc.py
mkdir -p logs
sbatch run_control_region_orih_qc.slurm
```

For interactive use (strict set):

```bash
python3 check_genbank_control_region_mitos2_orih.py \
  --dloop data/kingdom_genbank_annotations_split_new.csv \
  --misc data/genbank_misc_ftr.csv \
  --mitos-summary data/mitos2_species_summary.tsv \
  --out-dir control_region_orih_qc/strict
```

Add `--expanded` to also include tentative AT-rich/long noncoding annotations (saved by SLURM in `control_region_orih_qc/expanded`).

The outputs in each directory are:
- `summary.json`: distinct GenBank species, source counts, MITOS2 presence, OriH and overlap status counts.
- `species_overlap.csv`: one row per normalized species with GenBank coordinates (including disjoint multipart regions), MITOS2 OriH coordinates, intersected bp and status.

Statuses distinguish `overlap`, `no_overlap`, `no_mitos2_orih`, `gff_not_found`, `missing_gff_path`, and `absent_from_mitos_summary`.

**Important:** The code detects MITOS2 OriH/OH by GFF feature type and attributes, and requires coordinate systems to be directly comparable. Species-name matching is not sufficient evidence that GenBank and MITOS2 used identical accessions, orientation or circular rotation. Therefore computed genomic intersection is *provisional* until assembly provenance is verified. The script does not infer OriH when MITOS2 provides none.

Input CSV files are large and remain on McCleary; this commit does not copy them into GitHub.

## Updated CR1 / CR2 region-level analysis

The default CLI now uses `control_region_orih_region_level.py` and produces:

- `control_region_orih_qc/region_level/region_overlap.csv`: each separately annotated control region, CR1/CR2 label when present, source, note, physical segments, overlap with OriH, and classification.
- `control_region_orih_qc/region_level/species_overlap.csv`: one row per organism, CR1 and CR2 counts, duplicated-CR candidate flag, count of full regions overlapping OriH, overlapping region labels, and QC status.
- `control_region_orih_qc/region_level/invalid_split_dloop.csv`: incomplete/malformed split D-loop groups.
- `control_region_orih_qc/region_level/summary.json`: aggregate counts and cautions.

Full and putative control regions are used for the primary OriH overlap. Short functional subfeatures (CSB/TAS/core), remnants and expanded AT-rich candidates are **retained for audit but excluded from primary full-region overlap**. CR1/CR2 are never joined into a single continuous interval across distant genomic positions. Multi-segment annotations are retained as multiple source intervals belonging to a single region. When multiple GenBank accessions exist for a species, positional comparisons are marked ambiguous rather than silently merged.

Run tests with `python3 test_control_region_orih_region_level.py`. The legacy parser helpers remain importable, but `check_genbank_control_region_mitos2_orih.py` now dispatches to the region-level CLI. Previous strict/expanded output directories are left unchanged.

**Limitations:** Repeated or coextensive annotations from D-loop and misc_feature can describe the same underlying region; region counts are *annotation counts*, not automatically de-duplicated biological copies. Identify unique CR copies by accession plus interval and supporting locus annotations before interpreting duplication as biological evidence. Matching GenBank/MITOS2 assemblies and coordinate origins remains a prerequisite for biological overlap interpretation.
