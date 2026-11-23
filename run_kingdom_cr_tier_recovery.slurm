#!/bin/bash
#SBATCH --job-name=kingdom_cr_tiers
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=logs/kingdom_cr_tiers_%j.out
#SBATCH --error=logs/kingdom_cr_tiers_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON="${PYTHON:-/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python}"
"${PYTHON}" test_genbank_control_region_vocabulary.py
"${PYTHON}" test_kingdom_cr_tier_recovery.py
"${PYTHON}" audit_genbank_control_region_vocabulary.py \
  --genbank data/kingdom_genbank_annotations_split_new.csv \
  --misc data/genbank_misc_ftr.csv \
  --kingdom-tsv longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv \
  --out-dir control_region_vocabulary_audit
"${PYTHON}" analyze_kingdom_cr_tier_recovery.py \
  --genbank data/kingdom_genbank_annotations_split_new.csv \
  --misc data/genbank_misc_ftr.csv \
  --longest-tsv longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv \
  --mitos-summary data/mitos2_species_summary.tsv \
  --out-dir control_region_vocabulary_audit/kingdom_recovery
