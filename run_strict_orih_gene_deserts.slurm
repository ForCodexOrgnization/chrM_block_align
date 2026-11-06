#!/bin/bash
#SBATCH --job-name=strict_OriH_desert
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/strict_OriH_desert_%j.out
#SBATCH --error=logs/strict_OriH_desert_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON="/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python"

"${PYTHON}" analyze_strict_orih_gene_deserts.py \
  --fasta-dir data/unrot_fasta \
  --genbank-annotations data/kingdom_genes_v2.csv \
  --gene-aliases data/genefiltfinal5_2026_09_24.csv \
  --trnascan-csv data/consolidated_trnascan-2.csv \
  --mitos2-summary data/mitos2_species_summary.tsv \
  --out-dir strict_orih_gene_desert_qc

"${PYTHON}" plot_strict_orih_gene_deserts.py \
  --input strict_orih_gene_desert_qc/all_ranks_species_strict_OriH.tsv \
  --out-dir strict_orih_gene_desert_qc/figures
