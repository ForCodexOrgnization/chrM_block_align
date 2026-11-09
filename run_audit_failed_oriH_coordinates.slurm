#!/bin/bash
#SBATCH --job-name=audit_OriH_coords
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/audit_OriH_coords_%j.out
#SBATCH --error=logs/audit_OriH_coords_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON="/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python"
"${PYTHON}" audit_failed_oriH_gff_coordinates.py \
  --failures strict_orih_gene_desert_qc/failure_diagnostics/OriH_predicted_but_not_selected.tsv \
  --ranks kingdom \
  --out-dir strict_orih_gene_desert_qc/coordinate_audit
