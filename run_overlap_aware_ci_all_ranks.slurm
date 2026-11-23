#!/bin/bash
#SBATCH --job-name=overlap_CI
#SBATCH --partition=day
#SBATCH --array=0-9%3
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/overlap_ci_%A_%a.out
#SBATCH --error=logs/overlap_ci_%A_%a.err
set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
PYTHON="${PYTHON:-python3}"
SOURCE_ROOT="${SOURCE_ROOT:-block_results_genbank_trnascan_mitos_fallback_split_v3/pad100}"
RANKS=(kingdom phylum subphylum superclass class superorder order suborder infraorder parvorder)
RANK="${RANKS[${SLURM_ARRAY_TASK_ID}]}"
echo "Scoring existing tRNAscan > GenBank > MITOS2 fallback pad100 alignment: ${SOURCE_ROOT}/${RANK}"

"${PYTHON}" calculate_overlap_aware_ci.py \
  --rank-dir "${SOURCE_ROOT}/${RANK}" \
  --human-fasta data/Humans_mtdna.fasta \
  --human-genbank data/sequence.gb
