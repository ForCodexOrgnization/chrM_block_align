#!/bin/bash
#SBATCH --job-name=all_gene_ci_corr
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/all_gene_ci_corr_%j.out
#SBATCH --error=logs/all_gene_ci_corr_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python
"${PYTHON}" compare_gene_ci_mitos2_genbank.py \
  --annotation-csv data/2026_10_10_mito_cons_anns_main.csv \
  --mitos2-gene-only-root block_results_v6 \
  --mitos2-padded-root block_results_padded_v1 \
  --genbank-root block_results_genbank_mitos_v1 \
  --metric CI_MSA \
  --out-dir gene_ci_comparison_all_pipelines
