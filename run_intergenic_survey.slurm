#!/bin/bash
#SBATCH --job-name=chrM_IGR_survey
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/igr_survey_%j.out
#SBATCH --error=logs/igr_survey_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python
"${PYTHON}" -c "import Bio, pandas, numpy; print('Dependencies OK')"
"${PYTHON}" survey_intergenic_lengths.py \
  --ranks class \
  --out-dir block_results_v6/intergenic_length_survey
