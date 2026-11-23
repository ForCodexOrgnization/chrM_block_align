#!/bin/bash
#SBATCH --job-name=genbank_orih_100
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/genbank_orih_100_%j.out
#SBATCH --error=logs/genbank_orih_100_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON="${PYTHON:-/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python}"
"${PYTHON}" test_compare_genbank_cr_to_longest_orih.py
"${PYTHON}" compare_genbank_cr_to_longest_orih.py \
  --longest-tsv longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv \
  --dloop data/kingdom_genbank_annotations_split_new.csv \
  --misc data/genbank_misc_ftr.csv \
  --mitos-summary data/mitos2_species_summary.tsv \
  --ranks kingdom phylum subphylum superclass class superorder order suborder infraorder parvorder \
  --out-dir longest_all_gene_desert_qc/genbank_CR_comparison

"${PYTHON}" diagnose_genbank_cr_orih_containment.py \
  --species-tsv longest_all_gene_desert_qc/genbank_CR_comparison/species_100pct_containment.tsv \
  --out-dir longest_all_gene_desert_qc/genbank_CR_comparison
