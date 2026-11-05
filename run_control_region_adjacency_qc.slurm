#!/bin/bash
#SBATCH --job-name=chrM_ctrl_adj
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=logs/control_region_adjacency_%j.out
#SBATCH --error=logs/control_region_adjacency_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON="${PYTHON:-/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python}"

# Uses resolved canonical annotations exactly like the existing pipelines.
# No MAFFT or cached MSA is needed; does not modify alignment outputs.
"${PYTHON}" analyze_control_region_adjacency.py --pipeline genbank \
  --output-dir control_region_adjacency_qc_genbank_trnascan
"${PYTHON}" analyze_control_region_adjacency.py --pipeline mitos2 \
  --output-dir control_region_adjacency_qc_mitos2
