#!/bin/bash
#SBATCH --job-name=chrM_4way_diag
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/four_way_diag_%j.out
#SBATCH --error=logs/four_way_diag_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python
"${PYTHON}" -c "import numpy,pandas,sklearn,matplotlib,Bio;print('dependencies OK')"
"${PYTHON}" analyze_four_way_diagnostics.py \
 --old-csv data/2026_10_10_mito_cons_anns_main.csv \
 --gene-only-dir block_results_v6 \
 --padded-dir block_results_padded_v1 \
 --human-genbank data/sequence.gb \
 --metric CI_MSA \
 --bootstrap-reps 1000 \
 --out-dir block_results_padded_v1/four_way_diagnostics
