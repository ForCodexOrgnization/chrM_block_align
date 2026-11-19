#!/usr/bin/env python3
"""Trace human gene-overlap MSA columns to species mtDNA coordinates.

Read-only QC for the selected tRNAscan > GenBank > MITOS2 fallback pad100
pipeline. Rebuilds the original selected annotations, NOT alignments.
Each coordinate is 1-based on the source species FASTA. A human-aligned
column containing a species gap has no genomic coordinate. Gene A and B
need not be adjacent or on the same strand in distant taxa.
"""
import argparse
import collections
import csv
import json
import re
from pathlib import Path

import chrM_block_align as base
from calculate_overlap_aware_ci import PATTERN
from human_genbank import load_human_genbank
from run_genbank_trnascan_mitos_fallback import (
    build_feature_sets, load_aliases, load_genbank_csv,
)

BASES = set("ACGT")


def gene_footprint(feature, length):
    return {i % length for i in range(feature.start, feature.end)}


def physical_overlap(first, second, length):
    return bool(gene_footprint(first, length) & gene_footprint(second, length))


def source_coordinates(feature, length, padding, left, right, human_strand):
    """Per-base coordinates in the *final* (human-forward) aligned orientation.

    Padded builder extracts a forward genomic interval, reverse-complements
    for negative species strand, and finally reverses negative HUMAN blocks.
    Each of those two strand flips reverses nucleotide-to-position order.
    """
    beg = feature.start - (padding if left else 0)
    end = feature.end + (padding if right else 0)
    coords = [i % length + 1 for i in range(beg, end)]
    if (feature.strand == "-") != (human_strand == "-"):
        coords.reverse()
    return coords


def project_columns(aln_path, gene, human_feature, species_features, sequences,
                    human_seq, target_positions, allowed_flanks, padding=100):
    """Retrieve calls and 1-based species positions at selected human sites."""
    match = PATTERN.fullmatch(aln_path.name)
    if not match or match.group(1) != gene or match.group(3) != str(padding):
        raise ValueError(f"Unexpected cached gene alignment: {aln_path}")
    aligned = {base.norm(k):v for k,v in base.fasta(aln_path)}
    if len(aligned) != len(set(aligned)):
        raise ValueError(f"Duplicate species in {aln_path}")
    if len({len(v) for v in aligned.values()}) != 1:
        raise ValueError(f"Unequal sequence lengths in {aln_path}")
    if "homosapiens" not in aligned:
        raise ValueError(f"Human record absent in {aln_path}")
    if human_feature.strand == "-":
        aligned = {k:base.rc(v) for k,v in aligned.items()}
    hbeg, hend = human_feature.start-padding, human_feature.end+padding
    hcoords = [i % len(human_seq)+1 for i in range(hbeg,hend)]
    if aligned["homosapiens"].replace("-", "") != base.circular(human_seq,hbeg,hend):
        raise ValueError(f"Human reference mismatch: {aln_path}")
    columns = {}
    j = 0
    for col, nt in enumerate(aligned["homosapiens"]):
        if nt == "-":
            continue
        if hcoords[j] in target_positions:
            columns[hcoords[j]] = col
        j += 1
    if j != len(hcoords) or set(columns) != set(target_positions):
        raise ValueError(f"Human coordinate projection incomplete in {aln_path}")

    result = {pos:{} for pos in target_positions}
    for species, aligned_seq in aligned.items():
        if species == "homosapiens":
            continue
        feature = species_features.get(species, {}).get(gene)
        if feature is None or species not in sequences:
            raise ValueError(f"Missing selected feature/FASTA for {species} {gene}")
        seq = sequences[species]
        left, right = allowed_flanks(species, gene)
        positions = source_coordinates(
            feature, len(seq), padding, left, right, human_feature.strand
        )
        expected = "".join(seq[p-1] for p in positions)
        # Orientation after all reversals follows the actual species strand
        # transformed into the human gene's reference orientation.
        if feature.strand != human_feature.strand:
            expected = base.rc(expected)
        if aligned_seq.replace("-", "") != expected:
            raise ValueError(f"MSA-to-FASTA validation failed for {species} {gene}")
        non_gap = 0
        col_to_coord = {}
        for col, ch in enumerate(aligned_seq):
            if ch != "-":
                if col in columns.values():
                    if non_gap >= len(positions):
                        raise ValueError("Coordinate overrun")
                    col_to_coord[col] = positions[non_gap]
                non_gap += 1
        if non_gap != len(positions):
            raise ValueError(f"Species length mismatch for {species} {gene}")
        for pos,col in columns.items():
            result[pos][species] = (aligned_seq[col], col_to_coord.get(col))
    return result, set(aligned) - {"homosapiens"}


def classify(feature_a, feature_b, genome_length, call_a, coord_a,
             call_b, coord_b):
    if feature_a is None or feature_b is None:
        return "missing_annotation"
    if call_a is None or call_b is None or call_a not in BASES or call_b not in BASES:
        return "not_jointly_callable"
    if coord_a is None or coord_b is None:
        raise ValueError("Callable nucleotide without genomic coordinate")
    if coord_a == coord_b:
        return "same_coordinate"
    if physical_overlap(feature_a, feature_b, genome_length):
        return "different_coordinate_with_gene_overlap"
    return "no_physical_gene_overlap"


def run(args):
    rank_dir=Path(args.rank_dir)
    rank=rank_dir.name
    fasta_path=Path(args.fasta_dir)/f"{rank}_chrm.faa"
    out=Path(args.out_dir or rank_dir/"coordinate_overlap_qc")
    out.mkdir(parents=True,exist_ok=True)
    seqs={}
    display={}
    for name,seq in base.fasta(fasta_path):
        key=base.norm(name)
        if key in seqs:
            raise ValueError(f"Duplicate normalized species {key}")
        seqs[key]=seq
        display[key]=name
    human=seqs["homosapiens"]
    human_genes=load_human_genbank(args.human_genbank,human,base.Feature,base.EXPECTED)
    order=sorted(human_genes,key=lambda g:human_genes[g].start % len(human))

    aliases=load_aliases(args.gene_aliases)
    reference=[(g,human_genes[g].strand) for g in order]
    gb,_,_,_=load_genbank_csv(args.genbank_annotations,aliases,reference)
    trna=base.read_trnas(args.trnascan_csv)
    with open(args.mitos2_summary,newline="",encoding="utf-8-sig") as fh:
        mitos={base.norm(row["species"]):row["result_gff"]
               for row in csv.DictReader(fh,delimiter="\t")
               if row.get("species") and row.get("result_gff")}
    selected,_,_=build_feature_sets(rank,fasta_path,gb,trna,mitos)

    # Restore the exact neighbour predicate used when the cached MSAs were made.
    from chrM_block_align_padded import neighbours
    def allowed_flanks(species,gene):
        fs=selected[species]
        return neighbours(fs,gene,human_genes,order,len(seqs[species]))

    pair=tuple(args.genes)
    if len(pair)!=2 or any(g not in human_genes for g in pair):
        raise ValueError("Specify exactly two annotated human genes")
    positions=sorted(set(args.positions))
    for p in positions:
        if not (1<=p<=len(human)):
            raise ValueError(f"Invalid human coordinate: {p}")
        if not all(human_genes[g].start <= p-1 < human_genes[g].end for g in pair):
            raise ValueError(f"Position {p} is not in both human gene cores")
    block_dir=rank_dir/"gene_alignments"
    calls={}
    present={}
    for gene in pair:
        files=[p for p in block_dir.glob("*.aln.fa")
               if (m:=PATTERN.fullmatch(p.name)) and m.group(1)==gene and m.group(3)==str(args.padding)]
        if len(files)!=1:
            raise ValueError(f"Expected one {gene} pad{args.padding} MSA, found: {files}")
        calls[gene],present[gene]=project_columns(
            files[0],gene,human_genes[gene],selected,seqs,human,
            positions,allowed_flanks,args.padding
        )

    a,b=pair
    rows=[]
    counts=collections.Counter()
    all_species=sorted(set(present[a])|set(present[b]))
    for pos in positions:
        for s in all_species:
            feat_a=selected.get(s,{}).get(a)
            feat_b=selected.get(s,{}).get(b)
            va,ca=calls[a][pos].get(s,(None,None))
            vb,cb=calls[b][pos].get(s,(None,None))
            status=classify(feat_a,feat_b,len(seqs[s]),va,ca,vb,cb)
            counts[(pos,status)]+=1
            rows.append(dict(
                Position=pos,Species=display[s],species_key=s,
                Gene_A=a,Gene_B=b,Base_A=va or "",Base_B=vb or "",
                Coord_A=ca or "",Coord_B=cb or "",
                Coord_same=int(ca==cb) if ca is not None and cb is not None else "",
                Physical_gene_overlap=(int(physical_overlap(feat_a,feat_b,len(seqs[s])))
                                       if feat_a is not None and feat_b is not None else ""),
                Gene_A_start=(feat_a.start+1 if feat_a else ""),
                Gene_A_end=(feat_a.end if feat_a else ""),
                Gene_B_start=(feat_b.start+1 if feat_b else ""),
                Gene_B_end=(feat_b.end if feat_b else ""),
                Status=status
            ))
    def write(path,rs):
        with open(path,"w",newline="") as fh:
            w=csv.DictWriter(fh,fieldnames=list(rs[0]))
            w.writeheader();w.writerows(rs)
    write(out/"species_coordinate_qc.tsv",rows) if False else None
    # CSV is easier to filter by position or species.
    write(out/"species_coordinate_qc.csv",rows)
    summary=[]
    for pos in positions:
        row={"Position":pos,"Gene_A":a,"Gene_B":b,"n_species":len(all_species)}
        for category in ("same_coordinate","different_coordinate_with_gene_overlap",
                         "no_physical_gene_overlap","missing_annotation",
                         "not_jointly_callable"):
            row[category]=counts[(pos,category)]
        summary.append(row)
    write(out/"coordinate_status_summary.csv",summary)
    (out/"summary.json").write_text(json.dumps({
        "rank":rank,"genes":list(pair),"positions":positions,
        "n_species_in_either_MSA":len(all_species),
        "caveat":"Coordinates are 1-based in source FASTA; comparisons are for matched species only. A gap has no coordinate. Physical overlap alone does not establish aligned-site orthology."
    },indent=2)+"\n")
    for row in summary:
        print(row)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rank-dir",required=True)
    p.add_argument("--genes",nargs=2,default=["ATP6","ATP8"])
    p.add_argument("--positions",nargs="+",type=int,default=[8527,8530,8531,8563])
    p.add_argument("--padding",type=int,default=100)
    p.add_argument("--fasta-dir",default="data/unrot_fasta")
    p.add_argument("--human-genbank",default="data/sequence.gb")
    p.add_argument("--genbank-annotations",default="data/kingdom_genes_v2.csv")
    p.add_argument("--gene-aliases",default="data/genefiltfinal5_2026_09_24.csv")
    p.add_argument("--trnascan-csv",default="data/consolidated_trnascan-2.csv")
    p.add_argument("--mitos2-summary",default="data/mitos2_species_summary.tsv")
    p.add_argument("--out-dir")
    args=p.parse_args()
    run(args)

if __name__=="__main__":
    main()
