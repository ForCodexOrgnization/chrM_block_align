# Species-coordinate overlap QC

Read-only coordinate tracing for cached tRNAscan > GenBank > MITOS2 fallback / pad100 MAFFT blocks.
Reconstructs source annotation priority and maps species-aligned nucleotides to their original 1-based FASTA coordinates. No MAFFT rerun and no changes to CI files.

## Run on McCleary

    git pull origin main
    python3 test_overlap_genomic_coordinates.py
    mkdir -p logs
    sbatch run_overlap_genomic_coordinates_kingdom.slurm

Or customize:

    python3 analyze_overlap_genomic_coordinates.py --rank-dir block_results_genbank_trnascan_mitos_fallback_v2/pad100/kingdom --genes ATP6 ATP8 --positions 8527 8530 8531 8563

## Outputs

- coordinate_overlap_qc/species_coordinate_qc.csv: one row per human position and species with both gene-context calls and coordinates, selected gene boundaries, physical overlap and QC status.
- coordinate_overlap_qc/coordinate_status_summary.csv: category counts by human position.
- coordinate_overlap_qc/summary.json: provenance and caveats.

Statuses: same_coordinate; different_coordinate_with_gene_overlap; no_physical_gene_overlap; missing_annotation; not_jointly_callable.
Only jointly ACGT-callable positions are assigned a genomic coordinate comparison; gaps have no coordinate. Original source FASTA coordinates are circular and 1-based. Matching physical bases do not prove cross-species orthology.
Input genes must truly overlap at each specified human position; padding does not count as an annotated overlap.
Results remain exploratory until reviewed on McCleary.
