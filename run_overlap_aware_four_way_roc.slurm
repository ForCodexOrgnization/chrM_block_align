#!/bin/bash
#SBATCH --job-name=overlap_4roc
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/overlap_4roc_%j.out
#SBATCH --error=logs/overlap_4roc_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
PYTHON="${PYTHON:-/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python}"
ROOT="block_results_genbank_trnascan_mitos_fallback_v2/pad100"
for METRIC in CI_MSA CI_callable; do
  "${PYTHON}" compare_overlap_aware_four_way_roc.py \
    --results-root "${ROOT}" \
    --metric "${METRIC}" \
    --bootstrap-reps 1000 \
    --min-per-class 5
done
