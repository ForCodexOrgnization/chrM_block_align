#!/bin/bash
#SBATCH --job-name=overlap_qc
#SBATCH --partition=day
#SBATCH --array=0-9%3
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/overlap_qc_%A_%a.out
#SBATCH --error=logs/overlap_qc_%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON="${PYTHON:-/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python}"
SOURCE_ROOT="${SOURCE_ROOT:-block_results_genbank_trnascan_mitos_fallback_v2/pad100}"
RANKS=(kingdom phylum subphylum superclass class superorder order suborder infraorder parvorder)
RANK="${RANKS[${SLURM_ARRAY_TASK_ID}]}"
"${PYTHON}" analyze_overlap_alignment_qc.py \
  --rank-dir "${SOURCE_ROOT}/${RANK}" \
  --human-fasta data/Humans_mtdna.fasta \
  --human-genbank data/sequence.gb \
  --min-shared-pct 50
