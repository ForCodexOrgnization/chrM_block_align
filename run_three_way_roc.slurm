#!/bin/bash
#SBATCH --job-name=chrM_3way_ROC
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/three_way_roc_%j.out
#SBATCH --error=logs/three_way_roc_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python
"${PYTHON}" -c "import pandas,numpy,matplotlib,sklearn,Bio; print('Dependencies OK')"
"${PYTHON}" compare_three_way_roc.py \
  --rank class \
  --old-csv data/2026_10_10_mito_cons_anns_main.csv \
  --gene-only-csv block_results_v6/class/conservation/conservation_ci.csv \
  --padded-csv block_results_padded_v1/pad50/class/conservation/conservation_ci.csv \
  --human-genbank data/sequence.gb \
  --metric CI_MSA \
  --bootstrap-reps 1000 \
  --out-dir block_results_padded_v1/pad50/class/roc_evaluation
