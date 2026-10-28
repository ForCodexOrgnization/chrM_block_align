#!/usr/bin/env python3
"""Per-human-coordinate conservation from a reference-projected stitched MSA.

CI_MSA reproduces ci_freq_counts_calc.py: percentage of all species' characters
equal to the human reference base (including gaps/ambiguous IUPAC in denominator).
Base frequencies reproduce that script: A/C/G/T/- over those five characters.
Additional callable-only CI excludes gaps and ambiguous bases.
"""
import argparse
import collections
import csv
import json
from pathlib import Path

IUPAC=set("ACGTRYSWKMBDHVN-")
BASES="ACGT-"

def read_fasta(path):
    label=None;parts=[]
    with open(path) as handle:
        for line in handle:
            line=line.strip()
            if not line:continue
            if line.startswith(">"):
                if label is not None:yield label,"".join(parts).upper()
                label=line[1:].strip();parts=[]
            else:parts.append(line)
    if label is not None:yield label,"".join(parts).upper()

def calculate(msa_path,out_dir,map_path=None,expected_length=16569):
    records=list(read_fasta(msa_path))
    if not records:raise ValueError("Empty alignment")
    names=[n for n,_ in records]
    if len(set(names))!=len(names):raise ValueError("Duplicate FASTA headers")
    human=[(n,s) for n,s in records if n.split("|")[0].strip().lower()=="homo sapiens"]
    if len(human)!=1:raise ValueError("Expected exactly one Homo sapiens in FASTA")
    ref=human[0][1]
    length=len(ref)
    if length!=expected_length:raise ValueError(f"Reference length {length} != {expected_length}")
    if set(ref)-set("ACGTN"):
        raise ValueError("Unexpected characters in human reference")
    for name,seq in records:
        if len(seq)!=length:raise ValueError(f"{name}: length {len(seq)} != {length}")
        if set(seq)-IUPAC:
            raise ValueError(f"{name}: contains invalid symbols {set(seq)-IUPAC}")
    if map_path:
        with open(map_path,newline="") as handle:
            rows=list(csv.DictReader(handle,delimiter="\t"))
        positions=[int(r["human_position"]) for r in rows]
        if positions!=list(range(1,length+1)):
            raise ValueError("Coordinate map must be contiguous 1..reference length")
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    ci_file=out/"conservation_ci.csv"
    freq_file=out/"nucleotide_frequencies.csv"
    n=len(records)
    total_callable=0
    with ci_file.open("w",newline="") as ci_h, freq_file.open("w",newline="") as freq_h:
        ci=csv.writer(ci_h);freq=csv.writer(freq_h)
        ci.writerow(["Position","Homo sapiens","CI_MSA","CI_callable","n_species",
                     "n_callable","n_gap","n_ambiguous","n_match"])
        freq.writerow(["Position","a","c","g","t","-","n_frequency_valid",
                       "count_a","count_c","count_g","count_t","count_gap","count_ambiguous"])
        seqs=[seq for _,seq in records]
        for p,base in enumerate(ref):
            counts=collections.Counter(seq[p] for seq in seqs)
            call=sum(counts[b] for b in "ACGT")
            match=counts[base] if base in "ACGT" else 0
            ci_all=100*match/n
            ci_call=100*match/call if call else ""
            valid=sum(counts[b] for b in BASES)
            total_callable+=call
            ci.writerow([p+1,base,f"{ci_all:.8f}",
                         f"{ci_call:.8f}" if call else "",n,call,
                         counts["-"],n-call-counts["-"],match])
            freq.writerow([p+1]+[f"{100*counts[b]/valid:.8f}" if valid else "0" for b in BASES]
                          +[valid]+[counts[b] for b in BASES]+[n-valid])
    report={"input":str(msa_path),"reference_length":length,"n_species":n,
            "ci_definition":"CI_MSA = 100 * reference-base matches / all species (gap and IUPAC ambiguity in denominator)",
            "ci_callable_definition":"CI_callable = 100 * reference-base matches / ACGT species; empty if none",
            "frequency_denominator":"A/C/G/T/gap characters only, excluding ambiguous IUPAC bases",
            "mean_callable_fraction":total_callable/(n*length),
            "outputs":[str(ci_file),str(freq_file)]}
    (out/"conservation_summary.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    return report

def main():
    p=argparse.ArgumentParser(description="Calculate reference-anchored CI and nucleotide frequencies")
    p.add_argument("--msa",required=True,help="Rank stitched FASTA")
    p.add_argument("--out-dir",required=True)
    p.add_argument("--coordinate-map")
    p.add_argument("--expected-length",type=int,default=16569)
    args=p.parse_args()
    calculate(args.msa,args.out_dir,args.coordinate_map,args.expected_length)

if __name__=="__main__":main()
