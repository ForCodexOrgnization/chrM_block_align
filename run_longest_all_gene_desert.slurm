#!/bin/bash
#SBATCH --job-name=longest_gene_desert
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=logs/longest_gene_desert_%j.out
#SBATCH --error=logs/longest_gene_desert_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON="/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python"

"${PYTHON}" analyze_longest_all_gene_desert.py \
  --fasta-dir data/unrot_fasta \
  --genbank-annotations data/kingdom_genes_v2.csv \
  --trnascan-csv data/consolidated_trnascan-2.csv \
  --mitos2-summary data/mitos2_species_summary.tsv \
  --out-dir longest_all_gene_desert_qc

"${PYTHON}" plot_longest_all_gene_desert.py \
  --input longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv \
  --out-dir longest_all_gene_desert_qc/figures
