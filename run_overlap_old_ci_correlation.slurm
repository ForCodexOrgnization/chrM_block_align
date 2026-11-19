#!/bin/bash
#SBATCH --job-name=overlap_corr
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/overlap_corr_%j.out
#SBATCH --error=logs/overlap_corr_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
PYTHON="${PYTHON:-/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python}"
"${PYTHON}" compare_overlap_old_ci_correlation.py \
  --results-root block_results_genbank_trnascan_mitos_fallback_v2/pad100 \
  --old-csv data/2026_10_10_mito_cons_anns_main.csv
