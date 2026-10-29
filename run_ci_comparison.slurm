#!/bin/bash
#SBATCH --job-name=chrM_CI_compare
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ci_compare_%j.out
#SBATCH --error=logs/ci_compare_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

# Limit numerical-library threading and avoid interactive graphical backends.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPLBACKEND=Agg

# Existing Python environment on McCleary.
PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python

echo "Starting CI comparison: $(date)"
echo "Python: ${PYTHON}"
"${PYTHON}" -c "import pandas, numpy, matplotlib, sklearn; print('Dependencies OK')"

"${PYTHON}" compare_old_new_ci.py \
    --old-csv data/2026_10_10_mito_cons_anns_main.csv \
    --results-dir block_results_v6 \
    --new-metric CI_MSA \
    --out-dir block_results_v6/ci_comparison

echo "CI comparison completed: $(date)"
