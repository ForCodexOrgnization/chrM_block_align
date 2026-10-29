#!/bin/bash
#SBATCH --job-name=chrM_pad_class
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/chrm_pad_%j.out
#SBATCH --error=logs/chrm_pad_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
module load MAFFT/7.526-GCC-13.3.0-with-extensions
PYTHON=/vast/palmer/home.mccleary/njl27/shared/mt_calling_inputs/required_software/hail-conda/bin/python
PAD="${PAD:-50}"
OUT="block_results_padded_v1/pad${PAD}"
"${PYTHON}" chrM_block_align_padded.py --ranks class --padding "${PAD}" --output-dir "${OUT}" --threads "${SLURM_CPUS_PER_TASK}"
"${PYTHON}" validate_alignment.py --rank-dir "${OUT}/class" --human-fasta data/Humans_mtdna.fasta
"${PYTHON}" calculate_conservation.py --msa "${OUT}/class/class_stitched_aln.fasta" --coordinate-map "${OUT}/class/human_coordinate_map.tsv" --out-dir "${OUT}/class/conservation"
