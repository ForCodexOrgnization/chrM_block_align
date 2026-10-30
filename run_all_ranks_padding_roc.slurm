#!/bin/bash
#SBATCH --job-name=chrM_padROC
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/all_ranks_padding_roc_%j.out
#SBATCH --error=logs/all_ranks_padding_roc_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python
"${PYTHON}" -c "import pandas,numpy,sklearn,matplotlib;print('ROC dependencies OK')"
"${PYTHON}" compare_all_ranks_padding_roc.py \
 --old-csv data/2026_10_10_mito_cons_anns_main.csv \
 --padded-dir block_results_padded_v1 \
 --gene-only-dir block_results_v6 \
 --metric CI_MSA \
 --bootstrap-reps 1000 \
 --out-dir block_results_padded_v1/all_ranks_padding_roc
