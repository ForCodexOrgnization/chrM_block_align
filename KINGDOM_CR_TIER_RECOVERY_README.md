# Kingdom-wide hierarchical control region discovery and OriH evaluation

This analysis applies the reviewed GenBank vocabulary rules and then evaluates how much each new source expands CR candidate coverage over the **Kingdom cohort only**.

**Source inputs**
- `data/kingdom_genbank_annotations_split_new.csv` (D-loop annotations; split features reconstructed)
- `data/genbank_misc_ftr.csv` (`misc_kind=other`, including Note spellings)
- `longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv` (existing longest gene-free interval; not recomputed)
- `data/mitos2_species_summary.tsv` (parent MITOS2 OriH GFF)

**Ordered candidate classes**
1. `tier1_explicit`: D-loop, explicit control region including numbered CR1/CR2, typo `contral/conrol`, and putative CR (retained with raw note)
2. `tier1_review`: abbreviated CR1/CR2 or D-loop-like flags without explicit full-region name
3. `tier2_at_rich`: A+T-rich/AT-rich/AT-loop/AT-region alternatives
4. `tier2_review`: NCR, noncoding/intergenic candidates requiring manual evaluation
5. `tier3_longest`: existing longest unfiltered gene-free interval, retained even if no explicit CR annotation

Remnant and short subfeatures (CSB, TAS, SSB, CoRe) are excluded as standalone control regions; an explicit source D-loop feature is retained. Multiple CR copies at distinct genomic positions remain distinct. Identical regions on the same accession are de-duplicated. Candidate assignment is independent of MITOS2 OriH. A 100% OriH containment hit requires **one** individual candidate to cover the *entire union of OriH bases*; CR1/CR2 are not concatenated.

The cumulative candidate assignment counts include species without MITOS2 predictions. The cumulative OriH figures use only predicted-OriH species as their denominator, and report availability separately. The analysis is **provisional** unless GenBank reference and MITOS2 sequence/circular coordinates are independently verified. The current optional FASTA manifest can be used in `--require-verified-coordinates` mode; without it, there will be no verified OriH success (which is appropriate).

## McCleary

```bash
cd /vast/palmer/scratch/lake_nicole/lt692/chrM_block_align
git pull origin main
mkdir -p logs
python3 test_genbank_control_region_vocabulary.py
python3 test_kingdom_cr_tier_recovery.py
sbatch run_kingdom_cr_tier_recovery.slurm
```

## Outputs

`control_region_vocabulary_audit/kingdom_recovery/`:

- `kingdom_candidates.tsv`: every independently sourced candidate (CR1, CR2, alternative name, longest gap), original notes, interval and provisional OriH containment
- `kingdom_species_tier_recovery.tsv`: per Kingdom species source availability, first available annotation tier, fallback assignment, OriH status and individual/cumulative success flags
- `kingdom_cumulative_recovery.tsv`: per tier assignment counts, **newly added species**, cumulative 100% OriH containment counts, and new OriH-supported species
- `kingdom_cumulative_recovery.png` / `.pdf`: two percentage trajectories (all-species candidate assignment and OriH-predicted-species containment)
- `summary.json`: numeric summary and limitations

The vocabulary audit itself is written to `control_region_vocabulary_audit/annotation_vocabulary.tsv` and companion files.

Important limitations: A region matching MITOS2 OriH does not establish an experimentally active replication origin. A provisional coordinate match does not validate identical assemblies. **Do not treat 100% OriH recovery as independent annotation truth without FASTA and accession checks**. The gene-free interval is reused from earlier results, possibly built using older gene annotations. `tier2_review` terms (e.g. NCR, generic intergenic region) need additional gene-neighborhood QC before they are called biological control regions.
