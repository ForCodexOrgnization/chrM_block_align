# Overlap-aware CI: Mean / Min / Max (shared MAFFT alignments)

This is a **post-processing experiment**, not a change to the baseline stitching
algorithm. All three scoring choices read the **same existing gene-block MAFFT
alignments**, so MAFFT is not rerun.

## Motivation

Human mitochondrial genes overlap at several positions (for example,
ATP8/ATP6 and ND4L/ND4). The baseline `chrM_block_align.py` stitches each human
position from the first processed gene block, retaining alternative gene-block
calls only in a QC count. The padded method likewise chooses a single owner.
For highly divergent species, gene rearrangements can place the two homologs
at different genomic loci. Their two gene-specific alignments therefore need
not map to the same species nucleotide. This experiment does **not** force
species-level call agreement.

Instead, compute one **gene-context CI per human overlap position and gene**
and combine the resulting percentages. The output is one CI per human
reference position under each rule.

- **Mean:** arithmetic mean of available gene-context CIs.
- **Min:** minimum of available gene-context CIs.
- **Max:** maximum of available gene-context CIs.
- **Original:** untouched baseline conservation output for comparison.

Only **positions covered by at least two human GenBank genes** change.
All other positions keep the existing baseline CI string verbatim.

Both `CI_MSA` (fixed rank-wide species denominator, gaps and ambiguity count as
non-matches) and `CI_callable` (ACGT denominator, separately per gene) are
computed. For `CI_callable`, unavailable gene contexts are omitted from
Mean/Min/Max. If *none* are callable, the old baseline score is retained and
the QC table reports zero callable gene contexts. The same species may
contribute to both gene contexts, so the two CIs are **not** independent
measurements or separately counted species.

## Run on pre-existing class gene-only output

```bash
cd /vast/palmer/scratch/lake_nicole/lt692/chrM_block_align
git pull origin main
python3 test_overlap_aware_ci.py
python3 calculate_overlap_aware_ci.py \
  --rank-dir block_results_v6/class \
  --human-fasta data/Humans_mtdna.fasta \
  --human-genbank data/sequence.gb
```

The script requires:

1. `<rank-dir>/<rank>_stitched_aln.fasta` (unchanged)
2. `<rank-dir>/conservation/conservation_ci.csv` (original)
3. `<rank-dir>/gene_alignments/*_GENE_core.aln.fa` or
   `*_GENE_pad50.aln.fa` / `*_GENE_pad100.aln.fa` for every human
   overlap gene (cached MAFFT outputs)
4. `data/sequence.gb`, Biopython, and the matching human FASTA.

**Important:** If a baseline job ran without `--keep-intermediate`, the
`gene_alignments` directory may have been removed. The original stitched
FASTA does **not** retain discarded overlapping-gene calls; in that case
the intermediate alignments must be restored or regenerated. Do **not**
attempt to infer gene-context CI from the stitched FASTA alone.

For padded results, choose their rank directory, e.g.:

```bash
python3 calculate_overlap_aware_ci.py \
  --rank-dir block_results_padded_v1/pad50/class
```

For **all ten ranks** with existing complete alignments and baseline CI:

```bash
mkdir -p logs
sbatch run_overlap_aware_ci_all_ranks.slurm
# or reuse padded alignments:
SOURCE_ROOT=block_results_padded_v1/pad50 sbatch --export=ALL,SOURCE_ROOT run_overlap_aware_ci_all_ranks.slurm
```

The SLURM runner fails clearly for ranks missing an original CI or gene MSA;
it never rebuilds MAFFT alignments.

## Output

Created in `<rank-dir>/overlap_aware_ci/`:

| File | Description |
|---|---|
| `overlap_gene_context_ci.csv` | One row per human overlap position per gene: gene-specific CIs, callable count, match count |
| `overlap_position_ci.csv` | Every human position, overlap status, overlapping genes, Mean/Min/Max CI, callable gene-context count, callable CI range |
| `mean/conservation_ci.csv` | One row per human position, standard `Position`, `Homo sapiens`, `CI_MSA`, `CI_callable` columns |
| `min/conservation_ci.csv` | Same schema, minimum aggregator |
| `max/conservation_ci.csv` | Same schema, maximum aggregator |
| `summary.json` | Inputs, number of overlap positions, method and fallback definitions |

The three method-specific CSV files carry gene-context scores, **not**
stitched-FASTA base counts. They intentionally do not copy the baseline
`n_callable`, `n_match`, etc., since those counts would no longer describe
the aggregated CIs. Use the gene-context QC file to inspect denominators.

## Evaluation design

Compare **Original / Mean / Min / Max** on exactly the same labeled ClinVar
variants per taxonomic rank. Non-overlap positions are identical by design.
Report overlap-only variant counts; overlap-specific AUROC may be undefined or
unstable if there are too few pathogenic/benign labels. A genome-wide AUROC
difference may be small because only a limited set of human positions changes.

**Limitations:** Mean/Min/Max are exploratory summaries of two gene-specific
nucleotide conservation estimates, not direct proof of selection in two
reading frames. Alignment quality may differ between genes. Do not select
a final method based only on training/evaluation on the same ClinVar labels.
