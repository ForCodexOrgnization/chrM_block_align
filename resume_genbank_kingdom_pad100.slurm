#!/bin/bash
#SBATCH --job-name=gb_kingdom_pad100
#SBATCH --partition=day
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=logs/gb_kingdom_pad100_%j.out
#SBATCH --error=logs/gb_kingdom_pad100_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module load MAFFT/7.526-GCC-13.3.0-with-extensions
PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python
ROOT=block_results_genbank_mitos_v1
DIR="${ROOT}/pad100/kingdom"
if [[ ! -s "${DIR}/kingdom_stitched_aln.fasta" ]]; then
    "${PYTHON}" run_genbank_mitos_blocks.py \
      --genbank-annotations data/kingdom_genes_v2.csv \
      --gene-aliases data/genefiltfinal5_2026_09_24.csv \
      --ranks kingdom --approaches pad100 --threads "${SLURM_CPUS_PER_TASK}" \
      --output-dir "${ROOT}" --resume-incomplete
fi
for METHOD in gene_only pad50 pad100; do
    RDIR="${ROOT}/${METHOD}/kingdom"
    MSA="${RDIR}/kingdom_stitched_aln.fasta"
    MAP="${RDIR}/human_coordinate_map.tsv"
    CI="${RDIR}/conservation/conservation_ci.csv"
    if [[ ! -s "${MSA}" || ! -s "${MAP}" ]]; then
        echo "Missing completed stitched alignment or coordinate map: ${RDIR}" >&2
        exit 2
    fi
    if [[ ! -s "${CI}" ]]; then
        "${PYTHON}" validate_alignment.py --rank-dir "${RDIR}" --human-fasta data/Humans_mtdna.fasta
        "${PYTHON}" calculate_conservation.py --msa "${MSA}" --coordinate-map "${MAP}" --out-dir "${RDIR}/conservation"
    else
        echo "Keeping existing CI: ${CI}"
    fi
done
echo "Kingdom all three methods completed and CI files present."
