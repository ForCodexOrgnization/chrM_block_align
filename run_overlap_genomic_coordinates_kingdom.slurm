#!/bin/bash
#SBATCH --job-name=overlap_coords
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=logs/overlap_coords_%j.out
#SBATCH --error=logs/overlap_coords_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON="${PYTHON:-/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python}"
"${PYTHON}" analyze_overlap_genomic_coordinates.py \
  --rank-dir block_results_genbank_trnascan_mitos_fallback_v2/pad100/kingdom \
  --genes ATP6 ATP8 \
  --positions 8527 8530 8531 8563
