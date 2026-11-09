#!/bin/bash
#SBATCH --job-name=gb_coord_audit
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/gb_coord_audit_%j.out
#SBATCH --error=logs/gb_coord_audit_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON="/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python"

"${PYTHON}" audit_genbank_coordinate_parsing.py \
  --species-tsv strict_orih_gene_desert_qc/all_ranks_species_strict_OriH.tsv \
  --ranks kingdom \
  --only-zero-gap \
  --out-dir strict_orih_gene_desert_qc/genbank_coordinate_audit_zero_gap
