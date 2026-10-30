#!/bin/bash
#SBATCH --job-name=chrM_GB_MITOS
#SBATCH --partition=day
#SBATCH --array=0-9%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=logs/genbank_mitos_%A_%a.out
#SBATCH --error=logs/genbank_mitos_%A_%a.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module load MAFFT/7.526-GCC-13.3.0-with-extensions
PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python
RANKS=(kingdom phylum subphylum superclass class superorder order suborder infraorder parvorder)
RANK="${RANKS[${SLURM_ARRAY_TASK_ID}]}"
GENBANK_CSV="${GENBANK_CSV:-data/kingdom_genes_v2.csv}"
ALIASES_CSV="${ALIASES_CSV:-data/genefiltfinal5_2026_09_24.csv}"
OUT="block_results_genbank_mitos_v1"
for approach in gene_only pad50 pad100; do
  if [[ -e "${OUT}/${approach}/${RANK}" ]]; then
    echo "Existing output, refusing to overwrite: ${OUT}/${approach}/${RANK}" >&2
    exit 3
  fi
done
"${PYTHON}" run_genbank_mitos_blocks.py \
 --genbank-annotations "${GENBANK_CSV}" \
 --gene-aliases "${ALIASES_CSV}" \
 --ranks "${RANK}" \
 --threads "${SLURM_CPUS_PER_TASK}" \
 --output-dir "${OUT}"
for approach in gene_only pad50 pad100; do
  RDIR="${OUT}/${approach}/${RANK}"
  "${PYTHON}" validate_alignment.py --rank-dir "${RDIR}" --human-fasta data/Humans_mtdna.fasta
  "${PYTHON}" calculate_conservation.py --msa "${RDIR}/${RANK}_stitched_aln.fasta" \
    --coordinate-map "${RDIR}/human_coordinate_map.tsv" --out-dir "${RDIR}/conservation"
done
