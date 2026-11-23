# Kingdom-wide GenBank control-region vocabulary audit

Use the latest inputs:
- `data/kingdom_genbank_annotations_split_new.csv`
- `data/genbank_misc_ftr.csv`
- `longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv` (Kingdom species cohort)

The script reads *all D-loop features* and `misc_kind=other` features. It scans `Gene`, `G_Symbol`, `Note`, relevant product/name fields and `other_qualifiers`, retaining the original Note for manual review. CDD annotations are omitted. This is a nomenclature audit, **not** a proof of functional control regions.

Classification hierarchy:

1. `explicit_cr`, `putative_cr`: D-loop, control region, numbered CR (when explicitly written as control region), known typos such as `contral region`.
2. `numbered_cr_abbreviation`, `source_dloop_like`: review-tier candidate based on short CR1/CR2 abbreviation or upstream D-loop-like flag.
3. `at_rich_alternative`: A+T-rich, A-T-rich, AT-rich, AT-loop, AT-region.
4. `long_noncoding_alternative`: major/long/large noncoding or intergenic region.
5. `subfeature`: CSB/TAS/CoRe/SSB and related control-region motifs; not a whole control region.
6. `remnant`: pseudo or remnant CR; not primary full-region evidence.
7. `replication_origin_only`: origin annotation alone is not a control-region boundary.
8. `unrelated`: excluded.

Run:

```bash
git pull origin main
mkdir -p logs
python3 test_genbank_control_region_vocabulary.py
sbatch run_control_region_vocabulary_audit.slurm
```

Or:

```bash
python3 audit_genbank_control_region_vocabulary.py \
  --genbank data/kingdom_genbank_annotations_split_new.csv \
  --misc data/genbank_misc_ftr.csv \
  --kingdom-tsv longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv \
  --out-dir control_region_vocabulary_audit
```

Outputs:
- `annotation_vocabulary.tsv`: raw Note, annotation category, counts, unique species/accessions and example text.
- `candidate_annotations.tsv`: individual candidate annotations with original location and evidence tier.
- `species_annotation_tiers.tsv`: cohort-wide presence flags and best **annotation-only** tier.
- `summary.json`: counts, including misspellings.

`--all-genbank-species` audits all species in source files regardless of Kingdom FASTA cohort.

**Important:** An abbreviated `CR1` can mean other things outside mitochondrial context and is accordingly review-tier rather than automatically definitive. The scan is not intended to replace the existing region-level CR deduplication or OriH coordinate-checking pipeline. Review new vocabulary terms before adding them to the final candidate selection pipeline. After review, separately compare independent gene-free candidates and annotated control-region candidates against MITOS2 OriH on validated assembly-matched coordinates. Do not use MITOS2 OriH both to construct and independently validate the same boundaries.
