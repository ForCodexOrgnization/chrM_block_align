# Longest vs GenBank D-loop vs GenBank CR — exact OriH containment

This is a **separate** four-bar analysis and does not overwrite the previous Longest / GenBank CR / Both / Either plots.

Each candidate succeeds when **a single candidate interval** contains **100% of all MITOS2 parent OriH bases**. Intervals from CR1 and CR2, or different methods, are never concatenated to claim success.

Methods:

- `Longest`: previously calculated longest gene-free interval
- `GenBank D-loop`: **only** `class=D-loop` features from `data/kingdom_genbank_annotations_split_new.csv` (multi-part annotations reconstructed)
- `GenBank CR`: the expanded GenBank set through **Tier 2 review**: explicit CR, Tier 1 review, AT-rich alternatives, and Tier 2 review; includes D-loop but excludes standalone replication origins, subfeatures, and remnants
- `Either`: logical OR of Longest, D-loop, and GenBank CR; not a union of their nucleotide intervals

The D-loop set is a subset of the expanded GenBank CR set by construction, so `Either` may be identical to `Longest OR GenBank CR`. The new D-loop column is designed to quantify the additional value of alternative GenBank annotations.

All ten original taxonomic ranks are included, and outputs use two denominators:
1. **All species** in the longest-gap rank input (species without MITOS2 OriH predictions count as unsuccessful in this denominator).
2. **MITOS2 OriH-predicted species** only.

**Run on McCleary**

```bash
cd /vast/palmer/scratch/lake_nicole/lt692/chrM_block_align
git pull origin main
mkdir -p logs
python3 test_compare_dloop_genbank_cr_longest_orih.py
sbatch run_compare_dloop_genbank_cr_longest_orih.slurm
```

**Output directory:** `longest_all_gene_desert_qc/genbank_Dloop_CR_comparison/`

- `orih_dloop_cr_all_species.png`, `.pdf`
- `orih_dloop_cr_orih_predicted.png`, `.pdf`
- `rank_100pct_containment.tsv`
- `species_100pct_containment.tsv`
- `summary.json`

Both figures separate title, subtitle, and legend and label the Y-axis with the relevant denominator. The original longest-gap results are reused rather than recomputed.

**Critical caveat:** The GenBank and MITOS2 reference accession, orientation and circular rotation have not been proven identical for all species; compare the source sequences before interpreting these provisional overlaps as validation of function. OriH is a computational prediction, not experimental replication initiation evidence.
