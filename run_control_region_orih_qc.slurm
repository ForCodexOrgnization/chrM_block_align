#!/bin/bash
#SBATCH --job-name=orih_control_qc
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/orih_control_qc_%j.out
#SBATCH --error=logs/orih_control_qc_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON="${PYTHON:-python3}"
"${PYTHON}" test_control_region_orih_qc.py
"${PYTHON}" check_genbank_control_region_mitos2_orih.py \
  --dloop data/kingdom_genbank_annotations_split_new.csv \
  --misc data/genbank_misc_ftr.csv \
  --mitos-summary data/mitos2_species_summary.tsv \
  --out-dir control_region_orih_qc/strict
"${PYTHON}" check_genbank_control_region_mitos2_orih.py \
  --dloop data/kingdom_genbank_annotations_split_new.csv \
  --misc data/genbank_misc_ftr.csv \
  --mitos-summary data/mitos2_species_summary.tsv \
  --expanded \
  --out-dir control_region_orih_qc/expanded
