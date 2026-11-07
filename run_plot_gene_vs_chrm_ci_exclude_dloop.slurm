#!/bin/bash
#SBATCH --job-name=ci_no_dloop
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ci_no_dloop_%j.out
#SBATCH --error=logs/ci_no_dloop_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

PYTHON="/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python"

"${PYTHON}" plot_gene_vs_chrm_ci_exclude_dloop.py \
  --genbank-root block_results_genbank_trnascan_mitos_fallback_v2 \
  --metric CI_MSA \
  --out-dir gene_vs_chrm_ci_no_dloop
