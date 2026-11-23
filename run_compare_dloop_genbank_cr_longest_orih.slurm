#!/bin/bash
#SBATCH --job-name=dloop_cr_orih
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/dloop_cr_orih_%j.out
#SBATCH --error=logs/dloop_cr_orih_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON="${PYTHON:-/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python}"
"${PYTHON}" test_compare_dloop_genbank_cr_longest_orih.py
"${PYTHON}" compare_dloop_genbank_cr_longest_orih.py \
  --longest-tsv longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv \
  --dloop data/kingdom_genbank_annotations_split_new.csv \
  --misc data/genbank_misc_ftr.csv \
  --mitos-summary data/mitos2_species_summary.tsv \
  --ranks kingdom phylum subphylum superclass class superorder order suborder infraorder parvorder \
  --out-dir longest_all_gene_desert_qc/genbank_Dloop_CR_comparison
