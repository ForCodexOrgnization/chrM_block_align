#!/bin/bash
#SBATCH --job-name=gb_CI_backfill
#SBATCH --partition=day
#SBATCH --array=0-9%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=logs/genbank_ci_backfill_%A_%a.out
#SBATCH --error=logs/genbank_ci_backfill_%A_%a.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python
RANKS=(kingdom phylum subphylum superclass class superorder order suborder infraorder parvorder)
RANK="${RANKS[${SLURM_ARRAY_TASK_ID}]}"
ROOT="${GENBANK_RESULTS_ROOT:-block_results_genbank_mitos_v1}"
for METHOD in gene_only pad50 pad100; do
  DIR="${ROOT}/${METHOD}/${RANK}"
  CI="${DIR}/conservation/conservation_ci.csv"
  MSA="${DIR}/${RANK}_stitched_aln.fasta"
  MAP="${DIR}/human_coordinate_map.tsv"
  if [[ -s "${CI}" ]]; then
    echo "Already exists, skipping: ${CI}"
    continue
  fi
  if [[ ! -s "${MSA}" || ! -s "${MAP}" ]]; then
    echo "MISSING INPUT for ${METHOD}/${RANK}: ${MSA} or ${MAP}" >&2
    exit 2
  fi
  "${PYTHON}" validate_alignment.py --rank-dir "${DIR}" --human-fasta data/Humans_mtdna.fasta
  "${PYTHON}" calculate_conservation.py --msa "${MSA}" --coordinate-map "${MAP}" --out-dir "${DIR}/conservation"
  [[ -s "${CI}" ]] || { echo "CI generation failed: ${CI}" >&2; exit 3; }
done
