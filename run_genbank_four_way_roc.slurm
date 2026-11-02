#!/bin/bash
#SBATCH --job-name=chrM_GB_ROC
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/genbank_four_way_roc_%j.out
#SBATCH --error=logs/genbank_four_way_roc_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python
"${PYTHON}" compare_genbank_four_way_roc.py \
 --old-csv data/2026_10_10_mito_cons_anns_main.csv \
 --results-dir block_results_genbank_mitos_v1 \
 --out-dir block_results_genbank_mitos_v1/roc_comparison \
 --metric CI_MSA \
 --bootstrap-reps 1000
