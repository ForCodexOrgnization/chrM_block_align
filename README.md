# chrM block alignment (v6)

Human-anchored mitochondrial alignment across species with rearranged gene order.

## Design

- MITOS2 GFF `gene` (CDS) and `rRNA` annotations; tRNAscan CSV for tRNA.
- Independently align **complete** homologous genes, even if their genomic spans overlap.
- Intergenic segments between non-overlapping reference genes are handled separately, and only when flanking gene adjacency is supported.
- Each reference nucleotide is emitted **once** in the final stitched alignment. Overlap duplicates are reported.
- The human reference is an actual 16,569-bp sequence, not padded coordinates.
- Only aligned reference columns are emitted into stitched FASTA; insertion columns are retained in intermediate gene MSAs. The stitched output is a human-coordinate projection, **not** a physical genome alignment.

## Requirements

Python 3.9+, MAFFT 7.526. On McCleary:

```bash
module load MAFFT/7.526-GCC-13.3.0-with-extensions
mkdir -p logs
sbatch run_class.slurm
```

Expected input layout:

```text
data/
  unrot_fasta/class_chrm.faa
  mitos2_species_summary.tsv
  consolidated_trnascan-2.csv
  human_chrM_16569.fa
```

For direct execution:

```bash
python3 chrM_block_align.py --fasta-dir data/unrot_fasta --annotation-summary data/mitos2_species_summary.tsv --trnascan-csv data/consolidated_trnascan-2.csv --human-fasta data/human_chrM_16569.fa --output-dir block_results_v6 --ranks class --threads 4 --keep-intermediate
python3 validate_alignment.py --rank-dir block_results_v6/class --human-fasta data/human_chrM_16569.fa
```

## Important caveats

- A missing or ambiguous ortholog is represented as a gap (unavailable) rather than treated as a substitution.
- Intergenic blocks with rearranged neighbors are not assumed homologous; memory guard limits large flank MSAs.
- Separate gene alignments preserve the complete overlap, but the final chrM projection chooses one aligned result per human nucleotide. It records overlap conflicts and does not automatically validate orthology in ambiguous duplications.
- The Homo sapiens sequence in every rank FASTA, the --human-fasta sequence, MITOS2 GFF, and tRNAscan coordinates must all refer to the same exact 16,569-bp human reference. No liftover from 16,565 bp is performed. The script checks FASTA identity and GFF region length; tRNAscan coordinates must be independently verified.
- Review generated QC and validate on real McCleary data before computing CI. This is research code, not a clinically validated pipeline.
