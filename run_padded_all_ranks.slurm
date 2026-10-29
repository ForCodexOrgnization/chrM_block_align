#!/bin/bash
#SBATCH --job-name=chrM_pad_all
#SBATCH --partition=week
#SBATCH --array=0-19%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/chrm_pad_all_%A_%a.out
#SBATCH --error=logs/chrm_pad_all_%A_%a.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module load MAFFT/7.526-GCC-13.3.0-with-extensions
PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python
RANKS=(kingdom phylum subphylum superclass class superorder order suborder infraorder parvorder)
RANK_INDEX=$((SLURM_ARRAY_TASK_ID % 10))
PAD_INDEX=$((SLURM_ARRAY_TASK_ID / 10))
RANK="${RANKS[$RANK_INDEX]}"
PADS=(50 100)
PAD="${PADS[$PAD_INDEX]}"
OUT="block_results_padded_v1/pad${PAD}/${RANK}"
if [[ -f "${OUT}/validation.json" && -f "${OUT}/conservation/conservation_ci.csv" ]]; then
  echo "SKIP existing validated output: ${OUT}"
  exit 0
fi
if [[ -d "${OUT}/gene_alignments" || -f "${OUT}/${RANK}_stitched_aln.fasta" ]]; then
  echo "STOP: partial or existing results at ${OUT}. Investigate before resubmission; no overwrite." >&2
  exit 3
fi
echo "RUN rank=${RANK} pad=${PAD} output=${OUT}"
"${PYTHON}" chrM_block_align_padded.py --ranks "${RANK}" --padding "${PAD}" --output-dir "block_results_padded_v1/pad${PAD}" --threads "${SLURM_CPUS_PER_TASK}"
"${PYTHON}" validate_alignment.py --rank-dir "${OUT}" --human-fasta data/Humans_mtdna.fasta
"${PYTHON}" calculate_conservation.py --msa "${OUT}/${RANK}_stitched_aln.fasta" --coordinate-map "${OUT}/human_coordinate_map.tsv" --out-dir "${OUT}/conservation"
echo "DONE ${RANK}"
