#!/bin/bash
#SBATCH --job-name=pairwise_cr
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=logs/pairwise_cr_%j.out
#SBATCH --error=logs/pairwise_cr_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

PYTHON="/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python"
PAIRWISE_DIR="/home/ak3238/project/Cons/chrm/pwa_rotated/pwa_aligned"
ROTINDEX="${ROTINDEX:-data/kingdom_rotindex.csv}"
MITOS2_SUMMARY="${MITOS2_SUMMARY:-data/mitos2_species_summary.tsv}"

if [ ! -f "${MITOS2_SUMMARY}" ]; then
  echo "ERROR: Missing MITOS2 GFF summary: ${MITOS2_SUMMARY}" >&2
  exit 1
fi

if [ ! -f "${ROTINDEX}" ]; then
  echo "ERROR: Missing rotation index: ${ROTINDEX}" >&2
  echo "Copy kingdom_rotindex.csv to data/ or submit with --export=ALL,ROTINDEX=/full/path/kingdom_rotindex.csv" >&2
  exit 1
fi

"${PYTHON}" compare_pairwise_cr_to_longest_intervals.py \
  --pairwise-dir "${PAIRWISE_DIR}" \
  --rotindex "${ROTINDEX}" \
  --mitos2-summary "${MITOS2_SUMMARY}" \
  --ranks kingdom phylum subphylum superclass class superorder order suborder infraorder parvorder \
  --min-cr-mapped-fraction 0.5 \
  --out-dir longest_all_gene_desert_qc/pairwise_CR_validation
