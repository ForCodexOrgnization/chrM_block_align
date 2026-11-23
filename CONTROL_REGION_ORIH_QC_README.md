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
