#!/bin/bash
#SBATCH --job-name=GB_TRI_FALL200
#SBATCH --partition=day
#SBATCH --array=0-9%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=logs/gb_trifall_v2_pad100_%A_%a.out
#SBATCH --error=logs/gb_trifall_v2_pad100_%A_%a.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
module load MAFFT/7.526-GCC-13.3.0-with-extensions

PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python

RANKS=(
  kingdom phylum subphylum superclass class
  superorder order suborder infraorder parvorder
)
RANK="${RANKS[${SLURM_ARRAY_TASK_ID}]}"

OUT="block_results_genbank_trnascan_mitos_fallback_v2"
RDIR="${OUT}/pad100/${RANK}"

# Keep existing completed results intact.
if [[ -f "${RDIR}/${RANK}_stitched_aln.fasta" ]]; then
  echo "Existing stitched alignment, refusing to overwrite: ${RDIR}" >&2
  exit 3
fi

"${PYTHON}" run_genbank_trnascan_mitos_fallback.py \
  --genbank-annotations data/kingdom_genes_v2.csv \
  --gene-aliases data/genefiltfinal5_2026_09_24.csv \
  --trnascan-csv data/consolidated_trnascan-2.csv \
  --mitos2-summary data/mitos2_species_summary.tsv \
  --ranks "${RANK}" \
  --approaches pad100 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --output-dir "${OUT}" \
  --resume-incomplete

"${PYTHON}" validate_alignment.py \
  --rank-dir "${RDIR}" \
  --human-fasta data/Humans_mtdna.fasta

"${PYTHON}" calculate_conservation.py \
  --msa "${RDIR}/${RANK}_stitched_aln.fasta" \
  --coordinate-map "${RDIR}/human_coordinate_map.tsv" \
  --out-dir "${RDIR}/conservation"
