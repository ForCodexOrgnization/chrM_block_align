#!/bin/bash
#SBATCH --job-name=chrM_blocks
#SBATCH --array=0-8%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/chrm_%A_%a.out
#SBATCH --error=logs/chrm_%A_%a.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module load MAFFT/7.526-GCC-13.3.0-with-extensions

# Class already completed. Run the nine remaining ranks.
RANKS=(
    kingdom
    phylum
    subphylum
    superclass
    superorder
    order
    suborder
    infraorder
    parvorder
)
RANK="${RANKS[${SLURM_ARRAY_TASK_ID}]}"

echo "Starting rank: ${RANK}"
python3 chrM_block_align.py \
  --fasta-dir data/unrot_fasta \
  --annotation-summary data/mitos2_species_summary.tsv \
  --trnascan-csv data/consolidated_trnascan-2.csv \
  --human-fasta data/Humans_mtdna.fasta \
  --human-genbank data/sequence.gb \
  --output-dir block_results_v6 \
  --ranks "${RANK}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --flank-batch-size 100 \
  --keep-intermediate

python3 validate_alignment.py \
  --rank-dir "block_results_v6/${RANK}" \
  --human-fasta data/Humans_mtdna.fasta

python3 calculate_conservation.py \
  --msa "block_results_v6/${RANK}/${RANK}_stitched_aln.fasta" \
  --coordinate-map "block_results_v6/${RANK}/human_coordinate_map.tsv" \
  --out-dir "block_results_v6/${RANK}/conservation"

echo "Completed rank: ${RANK}"
