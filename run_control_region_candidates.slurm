#!/bin/bash
#SBATCH --job-name=chrM_CR_candidate
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=logs/control_region_candidates_%j.out
#SBATCH --error=logs/control_region_candidates_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON="${PYTHON:-/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python}"

"${PYTHON}" analyze_control_region_candidates.py \
  --genbank-annotations data/kingdom_genes_v2.csv \
  --gene-aliases data/genefiltfinal5_2026_09_24.csv \
  --trnascan-csv data/consolidated_trnascan-2.csv \
  --mitos2-summary data/mitos2_species_summary.tsv \
  --min-gap-bp 300 \
  --output-dir control_region_candidate_qc
