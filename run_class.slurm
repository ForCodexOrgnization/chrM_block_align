#!/bin/bash
#SBATCH --job-name=chrM_block_class
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/chrM_%j.out
#SBATCH --error=logs/chrM_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module load MAFFT/7.526-GCC-13.3.0-with-extensions
python3 chrM_block_align.py \
  --fasta-dir data/unrot_fasta \
  --annotation-summary data/mitos2_species_summary.tsv \
  --trnascan-csv data/consolidated_trnascan-2.csv \
  --human-fasta data/Humans_mtdna.fasta \
  --output-dir block_results_v6 \
  --ranks class \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --keep-intermediate
python3 validate_alignment.py --rank-dir block_results_v6/class --human-fasta data/Humans_mtdna.fasta
