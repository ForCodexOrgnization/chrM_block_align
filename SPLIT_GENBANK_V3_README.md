# Split-aware GenBank annotations (v3)

## New input

`data/kingdom_genbank_annotations_split_new.csv` (provided on McCleary). This 323 MB source file is NOT committed to GitHub. It has a final `part_complete` column: `complete`, `part 1/2`, `part 2/2`, and sometimes `part 1/3` through `part 3/3`.

From the supplied CSV: 576,332 rows, 385 multipart loci, of which 382 have two pieces and three have three pieces. Some multipart features are D-loop/control region and are not used for CDS/tRNA/rRNA blocks.

## How it works

1. GenBank records sharing accession, normalized species, gene/source locus label and ordered part IDs are joined into ONE feature before gene-name resolution and fallback selection.
2. Each fragment's physical 0-based half-open genome interval is retained, and the gene sequence is formed by concatenating fragment sequences in part order. Only then is a minus-strand gene reverse-complemented.
3. The sequence between fragments is **not** part of the aligned gene. This matters for three-part ND5 loci with a discontinuous middle region.
4. Padded blocks add up to 100 bp at the two external sequence endpoints if neighboring genes support those flanks; no padding is inserted between parts.
5. Incomplete, mixed-strand, mismatched-length or contradictory sets of parts are excluded rather than silently treated as separate genes.
6. Existing tRNAscan > GenBank > MITOS2 fallback precedence is preserved; the human reference always uses the canonical human GenBank.
7. `split_loci_qc.tsv` contains reconstructed feature coordinates and segment counts; per-rank annotation QC and alignment validation continue to run.

## McCleary commands

    cd /vast/palmer/scratch/lake_nicole/lt692/chrM_block_align
    git pull origin main
    python3 test_split_genbank_annotations.py
    test -s data/kingdom_genbank_annotations_split_new.csv
    mkdir -p logs
    sbatch run_genbank_trnascan_mitos_fallback_pad100.slurm

The revised job writes to `block_results_genbank_trnascan_mitos_fallback_split_v3/pad100/{rank}` and never overwrites the original `block_results_genbank_trnascan_mitos_fallback_v2` outputs.

After the v3 MSA and CI jobs finish, downstream overlap-aware CI, alignment QC, AUROC and correlation runner defaults also target the new v3 output root. If you are intentionally analyzing older v2 outputs, override the source path in the CLI.

**Caution:** The source FASTA for each species must be the same unrotated mitochondrial sequence to which the CSV's `start`/`end` coordinates refer. Validate per-rank logs, split_loci_qc.tsv, and alignments before interpreting revised CI. Newly joined multipart genes can change cached alignments; do not reuse v2 scores as if they came from this annotation.
