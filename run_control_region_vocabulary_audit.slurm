#!/bin/bash
#SBATCH --job-name=cr_vocab_audit
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/cr_vocab_audit_%j.out
#SBATCH --error=logs/cr_vocab_audit_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON="${PYTHON:-python3}"
"${PYTHON}" test_genbank_control_region_vocabulary.py
"${PYTHON}" audit_genbank_control_region_vocabulary.py \
  --genbank data/kingdom_genbank_annotations_split_new.csv \
  --misc data/genbank_misc_ftr.csv \
  --kingdom-tsv longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv \
  --out-dir control_region_vocabulary_audit
