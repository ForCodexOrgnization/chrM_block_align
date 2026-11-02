#!/bin/bash
#SBATCH --job-name=gene_CI_corr
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/gene_ci_corr_%j.out
#SBATCH --error=logs/gene_ci_corr_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python
"${PYTHON}" plot_gene_ci_correlation.py \
 --annotation-csv data/2026_10_10_mito_cons_anns_main.csv \
 --new-root block_results_genbank_mitos_v1 \
 --methods gene_only pad50 pad100 \
 --out-dir block_results_genbank_mitos_v1/gene_ci_correlations
