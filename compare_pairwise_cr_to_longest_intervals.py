#!/usr/bin/env python3
"""Map the human control region through pairwise aligned FASTAs into unrotated
species mtDNA coordinates, and compare with the UNFILTERED longest gene-free
interval plus MITOS2 OriH overlap status.

Coordinates:
- Human CR: NC_012920 positions 16024..16569 plus 1..576, 1-based inclusive.
- Pairwise species sequence: ROTATED; ungapped species base 1 is the rotated
  sequence start.
- Match_Position in kingdom_rotindex.csv: 1-based ORIGINAL species base that
  becomes rotated position 1.
- Other gene/origin annotations and longest-gap TSV: ORIGINAL unrotated genome.
- Each mapped base is transformed individually: original=((rot+match-2)%L)+1.
- A pairwise alignment is accepted only when its ungapped sequences exactly
  match the human reference and the independently rotated species FASTA.

Important: We compare mapped species BASES only. Human-gap insertion columns
are NOT inferred as aligned to the human control region, and the min/max of
mapped coordinates is NOT used to claim coverage across the circular origin.

Output for every rank:
- detailed per species TSV
- summary counts by OriH status, mapped-CR support
- overall summary for ranks
- QC file of missing/mismatched alignments and rotation index
"""
import argparse
import collections
import csv
import re
from pathlib import Path

CR_END = 576
CR_START = 16024
HUMAN_LENGTH = 16569
CR_LENGTH = (HUMAN_LENGTH - CR_START + 1) + CR_END
RANKS = ("kingdom", "phylum", "subphylum", "superclass", "class",
         "superorder", "order", "suborder", "infraorder", "parvorder")

def norm(x):
    return re.sub(r"[^a-z0-9]", "", str(x).casefold())

def fasta(path):
    name = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq).upper()

def write_tsv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def num(x):
    try:
        return int(x)
    except (ValueError, TypeError):
        return None

def original_pos(rotated_pos, match_position, length):
    """Convert 1-based ROTATED base to 1-based ORIGINAL base."""
    if not (1 <= rotated_pos <= length and 1 <= match_position <= length):
        raise ValueError("Invalid rotated position / Match_Position")
    return ((rotated_pos + match_position - 2) % length) + 1

def in_gap(pos, start, end, length):
    """Gap stored as 1-based inclusive [start, end], allowed to wrap."""
    if not (1 <= pos <= length and 1 <= start <= length and
            1 <= end <= length):
        return False
    return (start <= pos <= end) if start <= end else (
        pos >= start or pos <= end
    )

def rotate(seq, match_pos):
    if not 1 <= match_pos <= len(seq):
        raise ValueError("Match_Position is not 1-based and in sequence")
    k = match_pos - 1
    return seq[k:] + seq[:k]

def read_rotindex(path):
    data = {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        header = set(reader.fieldnames or [])
        required = {"Organism", "Sequence_Length", "Match_Position"}
        if not required <= header:
            raise ValueError(f"rotindex missing columns {required-header}")
        for row in reader:
            key = norm(row["Organism"])
            if not key:
                continue
            if key in data:
                raise ValueError(f"Duplicate rotation species: {row['Organism']}")
            data[key] = dict(length=num(row["Sequence_Length"]),
                             match=num(row["Match_Position"]),
                             organism=row["Organism"])
    return data

def read_longest(path):
    data = collections.defaultdict(dict)
    with open(path, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = (row["rank"], norm(row["species"]))
            if key in data:
                raise ValueError(f"Duplicate longest-gap record: {key}")
            data[key] = row
    return data

def mapped_human_cr(human_aln, species_aln, match_position, length):
    """Map each human control-region base to ORIGINAL species coordinates.

    Returns only species bases aligned to human CR nucleotides. A species
    insertion aligned to '-' in human is deliberately not included.
    """
    hpos = spos = 0
    mapped = set()
    h_cr_present = 0
    for h, s in zip(human_aln, species_aln):
        if h != "-":
            hpos += 1
        if s != "-":
            spos += 1
        if h != "-" and (hpos <= CR_END or hpos >= CR_START):
            h_cr_present += 1
            if s != "-":
                mapped.add(original_pos(spos, match_position, length))
    return mapped, h_cr_present, hpos, spos

def mapped_runs(bases, length):
    """Count separated mapped-base runs on a CIRCULAR genome."""
    if not bases:
        return 0
    if len(bases) == length:
        return 1
    return sum(1 for p in bases if (p - 2) % length + 1 not in bases)

def read_pairwise(path):
    recs = list(fasta(path))
    if len(recs) != 2:
        raise ValueError("expected exactly two FASTA records")
    human = [(name, seq) for name, seq in recs if norm(name) == "homosapiens"]
    if len(human) != 1:
        raise ValueError("expected exactly one Homo_sapiens record")
    species = [(name, seq) for name, seq in recs if norm(name) != "homosapiens"]
    if len(species) != 1:
        raise ValueError("expected exactly one species record")
    if len(human[0][1]) != len(species[0][1]):
        raise ValueError("pairwise aligned sequence lengths differ")
    return human[0][1], species[0][1], species[0][0]

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pairwise-dir",
                   default="/home/ak3238/project/Cons/chrm/pwa_rotated/pwa_aligned")
    p.add_argument("--rotindex", default="data/kingdom_rotindex.csv")
    p.add_argument("--unrot-fasta-dir", default="data/unrot_fasta")
    p.add_argument("--human-fasta", default="data/Humans_mtdna.fasta")
    p.add_argument("--longest-tsv", default=(
        "longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv"))
    p.add_argument("--ranks", nargs="+", choices=RANKS, default=list(RANKS))
    p.add_argument("--min-cr-mapped-fraction", type=float, default=0.5,
                   help="Minimum fraction of 1122 human CR bases mapped for a usable comparison")
    p.add_argument("--out-dir",
                   default="longest_all_gene_desert_qc/pairwise_CR_validation")
    a = p.parse_args()
    out = Path(a.out_dir)
    if not 0 <= a.min_cr_mapped_fraction <= 1:
        p.error("--min-cr-mapped-fraction must be from 0 to 1")

    human = list(fasta(a.human_fasta))
    if len(human) != 1 or len(human[0][1]) != HUMAN_LENGTH:
        raise ValueError("Expected a single 16569-base human reference")
    ref = human[0][1]
    rot = read_rotindex(a.rotindex)
    longest = read_longest(a.longest_tsv)

    # Indexed by file stem, so no folder scan inside each taxonomic rank.
    pairwise = {}
    for file in sorted(Path(a.pairwise_dir).glob("*.faa")):
        key = norm(file.stem)
        if key in pairwise:
            raise ValueError(f"Conflicting pairwise names: {pairwise[key]}, {file}")
        pairwise[key] = file
    if not pairwise:
        raise ValueError(f"No .faa alignments in {a.pairwise_dir}")

    rank_species = collections.defaultdict(list)
    seqs = {}
    lengths = {}
    for rank in a.ranks:
        path = Path(a.unrot_fasta_dir) / f"{rank}_chrm.faa"
        for name, seq in fasta(path):
            key = norm(name)
            if key == "homosapiens":
                continue
            if key in seqs and seqs[key] != seq:
                raise ValueError(f"Different unrotated sequences across ranks: {name}")
            if key in [k for k, _, _ in rank_species[rank]]:
                raise ValueError(f"Duplicate rank species: {rank}/{name}")
            seqs[key] = seq
            lengths[key] = len(seq)
            rank_species[rank].append((key, name, len(seq)))
    print("Unique species:", len(seqs), "pairwise FAA files:", len(pairwise),
          flush=True)

    mapped = {}
    for index, (key, original) in enumerate(seqs.items(), 1):
        row = dict(species_key=key, pairing_status="",
                   n_human_CR_mapped_bp=0, n_human_CR_total_bp=CR_LENGTH,
                   mapped_CR_fraction=0.0,
                   n_unique_mapped_species_bp=0, n_mapped_runs=0,
                   mapped_bases=None, species_alignment_label="")
        file = pairwise.get(key)
        ri = rot.get(key)
        if file is None:
            row["pairing_status"] = "missing_pairwise_file"
        elif ri is None:
            row["pairing_status"] = "missing_rotindex"
        elif ri["length"] != len(original):
            row["pairing_status"] = "rotindex_length_mismatch"
        elif ri["match"] is None or not 1 <= ri["match"] <= len(original):
            row["pairing_status"] = "invalid_1based_match_position"
        else:
            try:
                h_aln, s_aln, label = read_pairwise(file)
                h_ungapped = h_aln.replace("-", "")
                s_ungapped = s_aln.replace("-", "")
                if h_ungapped != ref:
                    raise ValueError("human_ungapped_not_exact_reference")
                if s_ungapped != rotate(original, ri["match"]):
                    raise ValueError("species_ungapped_not_exact_rotation")
                bases, n_h_cr, h_total, s_total = mapped_human_cr(
                    h_aln, s_aln, ri["match"], len(original))
                if h_total != HUMAN_LENGTH or s_total != len(original):
                    raise ValueError("unexpected_ungapped_lengths")
                if n_h_cr != CR_LENGTH:
                    raise ValueError("human_CR_reference_incomplete")
                row.update(pairing_status="validated", mapped_bases=bases,
                           species_alignment_label=label,
                           n_human_CR_mapped_bp=len(bases),
                           mapped_CR_fraction=round(len(bases)/CR_LENGTH, 6),
                           n_unique_mapped_species_bp=len(bases),
                           n_mapped_runs=mapped_runs(bases, len(original)))
            except (ValueError, OSError, IndexError) as e:
                row["pairing_status"] = f"invalid_alignment:{str(e)}"
        mapped[key] = row
        if index % 2000 == 0:
            print("Processed pairwise:", index, "/", len(seqs), flush=True)

    per_species = []
    per_rank = []
    by_status = []
    for rank in a.ranks:
        rank_rows = []
        for key, name, length in rank_species[rank]:
            g = longest.get((rank, key))
            m = mapped[key]
            result = dict(rank=rank, species=name,
                          genome_bp=length,
                          pairwise_status=m["pairing_status"],
                          n_human_CR_mapped_bp=m["n_human_CR_mapped_bp"],
                          human_CR_bp=CR_LENGTH,
                          CR_mapping_fraction=m["mapped_CR_fraction"],
                          n_mapped_CR_runs=m["n_mapped_runs"],
                          mapped_CR_usable=0,
                          candidate_defined=0, longest_gap_bp="",
                          gap_start_1based="", gap_end_1based="",
                          OriH_status="unknown",
                          CR_overlap_bp="", CR_overlap_any="",
                          mapped_CR_coverage_by_gap="",
                          longest_gap_fraction_covered_by_CR="",
                          support_class="unavailable")
            if g is not None:
                defined = num(g.get("candidate_defined")) == 1
                result.update(candidate_defined=int(defined),
                              longest_gap_bp=g.get("longest_gap_bp", ""),
                              gap_start_1based=g.get("gap_start_1based", ""),
                              gap_end_1based=g.get("gap_end_1based", ""))
                if num(g.get("mitos2_gff_available")) != 1:
                    result["OriH_status"] = "missing_MITOS2_GFF"
                elif (num(g.get("n_OriH_features")) or 0) == 0:
                    result["OriH_status"] = "no_OriH_prediction"
                elif defined and num(g.get("OriH_overlaps")) == 1:
                    result["OriH_status"] = "OriH_overlaps_longest"
                elif defined:
                    result["OriH_status"] = "OriH_predicted_elsewhere"
                else:
                    result["OriH_status"] = "OriH_predicted_no_longest_interval"
            if (g is not None and result["candidate_defined"] and
                    m["pairing_status"] == "validated"):
                usable = m["mapped_CR_fraction"] >= a.min_cr_mapped_fraction
                result["mapped_CR_usable"] = int(usable)
                start = num(g.get("gap_start_1based"))
                end = num(g.get("gap_end_1based"))
                if start is not None and end is not None:
                    bases = m["mapped_bases"]
                    n_hit = sum(in_gap(pos, start, end, length) for pos in bases)
                    result.update(
                        CR_overlap_bp=n_hit, CR_overlap_any=int(n_hit > 0),
                        mapped_CR_coverage_by_gap=round(n_hit/len(bases), 6)
                        if bases else "",
                        longest_gap_fraction_covered_by_CR=round(
                            n_hit/num(g["longest_gap_bp"]), 6)
                        if num(g.get("longest_gap_bp")) else "",
                    )
                    if usable:
                        ori = result["OriH_status"]
                        if ori == "OriH_overlaps_longest":
                            result["support_class"] = ("OriH+_CR+" if n_hit else "OriH+_CR-")
                        elif ori == "OriH_predicted_elsewhere":
                            result["support_class"] = ("OriH_elsewhere_CR+" if n_hit else "OriH_elsewhere_CR-")
                        elif ori == "no_OriH_prediction":
                            result["support_class"] = ("OriH_missing_CR+" if n_hit else "OriH_missing_CR-")
                        else:
                            result["support_class"] = ("OriH_unknown_CR+" if n_hit else "OriH_unknown_CR-")
                else:
                    result["pairwise_status"] = "invalid_gap_coordinates"
            rank_rows.append(result)
            per_species.append(result)

        def count(predicate):
            return sum(bool(predicate(x)) for x in rank_rows)
        valid = [x for x in rank_rows if x["mapped_CR_usable"]]
        per_rank.append(dict(
            rank=rank, n_total_species=len(rank_rows),
            n_pairwise_validated=count(lambda x: x["pairwise_status"] == "validated"),
            n_mapped_CR_usable=len(valid),
            n_longest_interval=count(lambda x: x["candidate_defined"] == 1),
            n_longest_with_CR_overlap=count(lambda x: x["mapped_CR_usable"] and x["CR_overlap_any"] == 1),
            n_longest_with_CR_full_mapped_containment=count(
                lambda x: x["mapped_CR_usable"] and x["mapped_CR_coverage_by_gap"] == 1),
            n_OriH_overlaps_longest_and_CR_overlap=count(lambda x: x["support_class"] == "OriH+_CR+"),
            n_OriH_overlaps_longest_but_not_CR=count(lambda x: x["support_class"] == "OriH+_CR-"),
            n_OriH_predicted_elsewhere_but_CR_overlap=count(lambda x: x["support_class"] == "OriH_elsewhere_CR+"),
            n_OriH_predicted_elsewhere_and_not_CR=count(lambda x: x["support_class"] == "OriH_elsewhere_CR-"),
            n_no_OriH_prediction_but_CR_overlap=count(lambda x: x["support_class"] == "OriH_missing_CR+"),
            n_no_OriH_prediction_and_not_CR=count(lambda x: x["support_class"] == "OriH_missing_CR-")
        ))
        for status, group in sorted(collections.Counter(x["support_class"] for x in rank_rows).items()):
            by_status.append(dict(rank=rank, support_class=status, n_species=group))
        print(rank, len(rank_rows), "mapped usable:", len(valid),
              "CR overlaps:", per_rank[-1]["n_longest_with_CR_overlap"],
              flush=True)

    species_fields = [
        "rank", "species", "genome_bp", "pairwise_status",
        "n_human_CR_mapped_bp", "human_CR_bp", "CR_mapping_fraction",
        "n_mapped_CR_runs", "mapped_CR_usable", "candidate_defined",
        "longest_gap_bp", "gap_start_1based", "gap_end_1based", "OriH_status",
        "CR_overlap_bp", "CR_overlap_any", "mapped_CR_coverage_by_gap",
        "longest_gap_fraction_covered_by_CR", "support_class"]
    write_tsv(out/"all_ranks_species_pairwise_CR_vs_longest.tsv",
              per_species, species_fields)
    write_tsv(out/"all_ranks_summary_pairwise_CR_vs_longest.tsv",
              per_rank, list(per_rank[0]))
    write_tsv(out/"all_ranks_support_categories.tsv",
              by_status, ["rank", "support_class", "n_species"])
    qc = [dict(species_key=k, pairwise_status=v["pairing_status"],
               pairwise_file=str(pairwise.get(k,"")),
               rotindex_Match_Position=rot.get(k, {}).get("match",""),
               rotindex_Sequence_Length=rot.get(k, {}).get("length",""),
               unrot_genome_length=lengths[k],
               n_human_CR_mapped_bp=v["n_human_CR_mapped_bp"],
               mapped_CR_fraction=v["mapped_CR_fraction"],
               n_mapped_runs=v["n_mapped_runs"])
          for k,v in mapped.items()]
    write_tsv(out/"pairwise_mapping_QC.tsv", qc, list(qc[0]))
    print("Wrote results to", out.resolve(), flush=True)
    print("Caution: full-length exact-match QC is conservative; review failures in"
          " pairwise_mapping_QC.tsv. A valid human CR mapping does not itself"
          " demonstrate control-region homology in distant taxa.")

if __name__ == "__main__":
    main()
