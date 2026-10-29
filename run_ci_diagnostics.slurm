#!/bin/bash
#SBATCH --job-name=chrM_CI_diagnostics
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/ci_diagnostics_%j.out
#SBATCH --error=logs/ci_diagnostics_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg

PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python
"${PYTHON}" -c "import Bio,pandas,numpy,matplotlib,sklearn; print('Dependencies OK')"
"${PYTHON}" analyze_ci_diagnostics.py \
  --comparison-dir block_results_v6/ci_comparison \
  --human-genbank data/sequence.gb \
  --bootstrap-reps 1000 \
  --out-dir block_results_v6/ci_comparison/diagnostics
