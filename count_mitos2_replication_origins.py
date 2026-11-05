#!/usr/bin/env python3
"""Count MITOS2 OriH and OriH-b predictions per taxonomic rank.

Uses the rank FASTAs and MITOS2 result.gff paths already used by the mtDNA
block-alignment pipeline. Does not use stitched alignment or GenBank calls.

Output:
  <output-dir>/mitos2_origin_summary.tsv
  <output-dir>/mitos2_origin_by_species.tsv
  <output-dir>/mitos2_origin_feature_examples.tsv

We distinguish OriH-b from OriH (do not count OriH-b as OriH).
Counts are based on MITOS2 GFF annotations, not validated origins.
"""
import argparse
import csv
import re
from pathlib import Path
from collections import Counter

from chrM_block_align import fasta, norm

RANKS = ("kingdom","phylum","subphylum","superclass","class",
         "superorder","order","suborder","infraorder","parvorder")
# Names in MITOS2 may be OriH-b, OriH_b, OriH b, Ori-b, or OH-b.
B_RE = re.compile(r"(?<![a-z0-9])(?:ori[\s_-]*h[\s_-]*b|ori[\s_-]*b|oh[\s_-]*b)(?![a-z0-9])",re.I)
# OriH, OH, heavy-strand origin (do not mistake OriH-b for OriH).
H_RE = re.compile(r"(?<![a-z0-9])(?:ori[\s_-]*h|oh)(?![a-z0-9])",re.I)
L_RE = re.compile(r"(?<![a-z0-9])(?:ori[\s_-]*l|ol)(?![a-z0-9])",re.I)

def read_summary(path):
    with open(path,encoding="utf-8-sig",newline="") as fh:
        rows=list(csv.DictReader(fh,delimiter="\t"))
    result={}
    for row in rows:
        key=norm(row["species"])
        gff=(row.get("result_gff") or "").strip()
        if not key or not gff:continue
        if key in result and result[key]!=gff:
            raise ValueError(f"Conflicting MITOS2 GFF paths for {row['species']}")
        result[key]=gff
    return result

def classify_gff(path):
    hits=set()
    examples=[]
    with open(path,encoding="utf-8",errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):continue
            fields=line.rstrip("\n").split("\t")
            if len(fields)<9:continue
            # GFF feature type and attributes; exclude source/seqid to avoid
            # matching species names or filenames accidentally.
            description=fields[2]+";"+fields[8]
            description=B_RE.sub(" ",description)
            # Determine B using ORIGINAL description to avoid collision with H.
            original=fields[2]+";"+fields[8]
            types=[]
            if B_RE.search(original):types.append("OriH-b")
            if H_RE.search(description):types.append("OriH")
            if L_RE.search(original):types.append("OriL")
            for typ in types:
                hits.add(typ)
                if len(examples)<20:
                    examples.append((typ,fields[2],fields[3],fields[4],fields[8]))
    return hits,examples

def write_tsv(path,columns,rows):
    with open(path,"w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=columns,delimiter="\t")
        w.writeheader()
        w.writerows(rows)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mitos2-summary",default="data/mitos2_species_summary.tsv")
    ap.add_argument("--fasta-dir",default="data/unrot_fasta")
    ap.add_argument("--output-dir",default="mitos2_origin_qc")
    ap.add_argument("--ranks",nargs="+",choices=RANKS,default=list(RANKS))
    a=ap.parse_args()
    out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    paths=read_summary(a.mitos2_summary)
    members={}
    for rank in a.ranks:
        records=list(fasta(Path(a.fasta_dir)/f"{rank}_chrm.faa"))
        keys=[norm(name) for name,_ in records if norm(name)!="homosapiens"]
        if len(keys)!=len(set(keys)):
            raise ValueError(f"Duplicate normalized species in {rank} FASTA")
        members[rank]=[(name,norm(name)) for name,_ in records if norm(name)!="homosapiens"]
    allkeys={key for pairs in members.values() for _,key in pairs}
    cache={};examples=[]
    for key in sorted(allkeys):
        gff=paths.get(key)
        if not gff:
            cache[key]=("missing_summary",set())
        elif not Path(gff).is_file():
            cache[key]=("missing_gff",set())
        else:
            hits,ex=classify_gff(gff)
            cache[key]=("read",hits)
            for typ,feature,start,end,attrs in ex:
                examples.append(dict(species_key=key,annotation=typ,feature_type=feature,
                                     start=start,end=end,attributes=attrs))
    details=[];summary=[]
    for rank,pairs in members.items():
        for name,key in pairs:
            status,hits=cache[key]
            details.append(dict(rank=rank,species=name,gff_status=status,
                                has_OriH=int("OriH" in hits),has_OriH_b=int("OriH-b" in hits),
                                has_both=int("OriH" in hits and "OriH-b" in hits),
                                has_either=int(bool({"OriH","OriH-b"}&hits)),
                                has_OriL=int("OriL" in hits),
                                mitos2_gff=paths.get(key,"")))
        group=[r for r in details if r["rank"]==rank]
        n=len(group);read=sum(r["gff_status"]=="read" for r in group)
        c=lambda k:sum(r[k] for r in group)
        both=c("has_both")
        summary.append(dict(rank=rank,n_species=n,n_gff_read=read,
            n_missing_gff_or_summary=n-read,
            n_OriH=c("has_OriH"),n_OriH_b=c("has_OriH_b"),
            n_both_OriH_and_OriH_b=both,
            n_either_OriH_or_OriH_b=c("has_either"),
            n_OriL=c("has_OriL"),
            pct_both_of_all=round(100*both/n,2) if n else "",
            pct_both_of_gff_read=round(100*both/read,2) if read else ""))
    write_tsv(out/"mitos2_origin_summary.tsv",list(summary[0]),summary)
    write_tsv(out/"mitos2_origin_by_species.tsv",list(details[0]),details)
    write_tsv(out/"mitos2_origin_feature_examples.tsv",
              ["species_key","annotation","feature_type","start","end","attributes"],examples)
    print("\t".join(["Rank","Species","GFF read","OriH","OriH-b","Both","Either","OriL"]))
    for x in summary:
        print("\t".join(str(x[k]) for k in (
          "rank","n_species","n_gff_read","n_OriH","n_OriH_b",
          "n_both_OriH_and_OriH_b","n_either_OriH_or_OriH_b","n_OriL")))
    print("Results:",out)
    print("Review mitos2_origin_feature_examples.tsv to confirm the exact MITOS2 GFF labels.")

if __name__=="__main__":
    main()
