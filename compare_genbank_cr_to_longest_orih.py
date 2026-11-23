#!/usr/bin/env python3
"""Compare 100% OriH containment by longest gene desert vs GenBank control region.

For each rank, denominators:
  OriH-predicted species: all species with at least one MITOS2 parent OriH
  All species: every species in the longest-gene-desert input rank

A method succeeds iff ONE of its candidate regions contains the ENTIRE
union of annotated OriH bases. CR1/CR2 cannot combine partial coverage;
Neither can longest plus CR combine partial coverage. Missing annotations
are failures for all-species/predicted-species success percentages, but are
tracked separately in per-species status. Input coordinates are assumed to
refer to the same unrotated mtDNA; without provenance checks results are
provisional. The reference longest-gap TSV may have been generated with old
GenBank gene annotations; this tool does NOT silently recompute longest gaps.
"""
import argparse
import collections
import csv
import json
from pathlib import Path

from analyze_longest_all_gene_desert import spans, norm
from compare_pairwise_cr_to_longest_intervals import (
    RANKS, read_orih_positions, read_longest, read_mitos2_paths, in_gap)
from control_region_orih_region_level import parse_dloop, parse_misc

FULL_TYPES = {"full_control_region", "putative_control_region"}
CATEGORIES = ("Longest", "GenBank_CR", "Both", "Either")


def write_tsv(path, rows, names):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer=csv.DictWriter(fh, fieldnames=names, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def deduplicate(regions):
    """Exact physical region duplicates across D-loop and misc_feature."""
    data={}
    priorities={"full_control_region":2,"putative_control_region":1}
    for region in regions:
        if region["classification"] not in FULL_TYPES:
            continue
        key=(region["species"], region["acc"], tuple(sorted(set(region["intervals"]))))
        if key not in data:
            item=dict(region)
            item["sources"]={region["source"]}
            data[key]=item
        else:
            obj=data[key]
            obj["sources"].add(region["source"])
            if priorities[region["classification"]]>priorities[obj["classification"]]:
                obj["classification"]=region["classification"]
            if obj["cr_label"]=="CR_unspecified" and region["cr_label"]!="CR_unspecified":
                obj["cr_label"]=region["cr_label"]
    grouped=collections.defaultdict(list)
    for region in data.values():
        grouped[region["species"]].append(region)
    return grouped


def coordinates_to_set(intervals, length):
    positions=set()
    for a,b in intervals:
        pieces=spans(a,b,length)
        if not pieces:
            return None
        for start,end in pieces:
            positions.update(range(start+1,end+1))
    return positions


def evaluate_regions(regions, ori_positions, length):
    """Do not join distinct CR copies to obtain 100% containment."""
    best=0
    full=[]
    if ori_positions is None:
        return None,[],0
    for region in regions:
        positions=coordinates_to_set(region["intervals"],length)
        if positions is None:
            continue
        bp=len(positions & ori_positions)
        if bp>best:
            best=bp
        if bp==len(ori_positions) and ori_positions:
            full.append(region)
    return bool(full),full,best


def evaluate_longest(row, ori_positions, length):
    if not row or row.get("candidate_defined")!="1" or not ori_positions:
        return False,0
    try:
        start=int(row["gap_start_1based"]);end=int(row["gap_end_1based"])
    except (ValueError,TypeError):
        return False,0
    if not (1<=start<=length and 0<=end<=length):
        return False,0
    bp=sum(bool(in_gap(p,start,end,length)) for p in ori_positions)
    return bp==len(ori_positions),bp


def compare(longest_tsv, dloop, misc, mitos_summary, out_dir, ranks=RANKS):
    longest=read_longest(longest_tsv)
    gff_paths=read_mitos2_paths(mitos_summary)
    dloop_regions,invalid=parse_dloop(dloop)
    regions=deduplicate(dloop_regions+parse_misc(misc))
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    ori_cache={}
    results=[]
    counts=[]
    for rank in ranks:
        rank_rows=[]
        for (r,species),gap in longest.items():
            if r!=rank:
                continue
            name=gap.get("species",species)
            length=int(gap.get("sequence_length") or 0)
            candidate=regions.get(species,[])
            accessions={x["acc"] for x in candidate if x["acc"]}
            n_regions=len(candidate)
            ori=None
            ori_status="missing_gff"
            gff_value=gff_paths.get(species)
            if not length:
                ori_status="missing_genome_length"
            elif gff_value:
                gff=Path(gff_value)
                if not gff.is_file():
                    ori_status="gff_not_found"
                else:
                    cachekey=(str(gff.resolve()),length)
                    if cachekey not in ori_cache:
                        try:
                            intervals,positions=read_orih_positions(gff,length)
                            ori_cache[cachekey]=("predicted" if intervals else "no_OriH",positions)
                        except (OSError,ValueError) as err:
                            ori_cache[cachekey]=("invalid_gff:"+str(err),None)
                    ori_status,ori=ori_cache[cachekey]
            long_success=cr_success=both=either=0
            long_bp=cr_bp=0
            matched=[]
            qc="provisional_coordinate_match"
            if len(accessions)>1:
                qc="multiple_genbank_accessions"
            if ori_status=="predicted":
                long_success,long_bp=evaluate_longest(gap,ori,length)
                cr_success,matched,cr_bp=evaluate_regions(candidate,ori,length)
                # Do not score ambiguous GenBank accessions as a valid CR hit.
                if len(accessions)>1:
                    cr_success=False;matched=[]
                both=int(long_success and cr_success)
                either=int(long_success or cr_success)
            rank_rows.append(dict(
                rank=rank,species=name,organism_key=species,
                genome_bp=length,ori_status=ori_status,
                ori_bp=len(ori) if ori is not None else "",
                longest_candidate_defined=gap.get("candidate_defined",""),
                n_genbank_control_regions=n_regions,
                n_genbank_accessions=len(accessions),
                coordinate_qc=qc if ori_status=="predicted" else "not_evaluable",
                Longest_full=int(bool(long_success)),
                GenBank_CR_full=int(bool(cr_success)),
                Both_full=both,Either_full=either,
                Longest_overlap_bp=long_bp if ori is not None else "",
                GenBank_CR_max_single_region_overlap_bp=cr_bp if ori is not None else "",
                GenBank_CR_containing_region_labels=";".join(
                    sorted({x["cr_label"] for x in matched})),
                GenBank_CR_containing_region_intervals=";".join(
                    ",".join(f"{a}-{b}" for a,b in x["intervals"])
                    for x in matched)))
        if not rank_rows:
            raise ValueError("No species found for rank: "+rank)
        for denominator in ("OriH_predicted_species","All_species"):
            subset=([r for r in rank_rows if r["ori_status"]=="predicted"]
                    if denominator=="OriH_predicted_species" else rank_rows)
            count=len(subset)
            totals={category:sum(r[category+"_full"] for r in subset)
                    for category in CATEGORIES}
            counts.append(dict(rank=rank,denominator=denominator,
                n_denominator=count,
                n_all_species=len(rank_rows),
                n_OriH_predicted=sum(r["ori_status"]=="predicted" for r in rank_rows),
                **{"n_"+c:totals[c] for c in CATEGORIES},
                **{"pct_"+c:round(100*totals[c]/count,4) if count else ""
                   for c in CATEGORIES}))
        results.extend(rank_rows)
    columns=list(results[0])
    write_tsv(out/"species_100pct_containment.tsv",results,columns)
    write_tsv(out/"rank_100pct_containment.tsv",counts,list(counts[0]))
    summary=dict(
        ranks=list(ranks),
        n_rank_species_rows=len(results),
        n_invalid_split_dloop=len(invalid),
        denominators=["OriH_predicted_species","All_species"],
        definition="All OriH bases must lie within ONE candidate region; candidates cannot concatenate partial coverage.",
        categories=list(CATEGORIES),
        caveat="GenBank/MITOS2 assembly and circular coordinate origin must be verified separately; these raw positional comparisons are provisional.",
        longest_source=str(longest_tsv))
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
    return results,counts


def plot(counts,out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    out=Path(out_dir)
    colors=("#2778B5","#F28B30","#57A05B","#D62929")
    labels=("Longest ∩ OriH","GenBank CR ∩ OriH","Both","Either")
    for denom,subtitle in (("OriH_predicted_species","OriH-predicted species"),
                           ("All_species","All species")):
        rows={r["rank"]:r for r in counts if r["denominator"]==denom}
        ranks=[rank for rank in RANKS if rank in rows]
        x=np.arange(len(ranks))
        fig,ax=plt.subplots(figsize=(12.8,5))
        w=.19
        for i,(c,label,color) in enumerate(zip(CATEGORIES,labels,colors)):
            y=[float(rows[r]["pct_"+c]) for r in ranks]
            ax.bar(x+(i-1.5)*w,y,w,label=label,color=color)
        ax.set_ylim(0,100)
        ax.set_yticks(range(0,101,20))
        ax.set_yticklabels([f"{i}%" for i in range(0,101,20)])
        ax.set_xticks(x)
        ax.set_xticklabels(ranks,rotation=35,ha="right")
        ax.set_xlabel("Taxonomic rank")
        ax.set_ylabel("Species (%)")
        ax.grid(axis="y",alpha=.25,linestyle="--")
        ax.set_axisbelow(True)
        ax.set_title("OriH containment across taxonomic ranks\n"
                     "100% OriH coverage | "+subtitle,loc="left",fontweight="bold")
        ax.legend(loc="upper center",bbox_to_anchor=(.5,1.16),ncol=4,
                  frameon=False,fontsize=9)
        fig.tight_layout()
        name=("orih_containment_genbank_predicted" if denom=="OriH_predicted_species"
              else "orih_containment_genbank_all_species")
        fig.savefig(out/(name+".png"),dpi=220,bbox_inches="tight")
        fig.savefig(out/(name+".pdf"),bbox_inches="tight")
        plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--longest-tsv",default="longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv")
    p.add_argument("--dloop",default="data/kingdom_genbank_annotations_split_new.csv")
    p.add_argument("--misc",default="data/genbank_misc_ftr.csv")
    p.add_argument("--mitos-summary",default="data/mitos2_species_summary.tsv")
    p.add_argument("--out-dir",default="longest_all_gene_desert_qc/genbank_CR_comparison")
    p.add_argument("--ranks",nargs="+",choices=RANKS,default=list(RANKS))
    a=p.parse_args()
    _,counts=compare(a.longest_tsv,a.dloop,a.misc,a.mitos_summary,a.out_dir,a.ranks)
    plot(counts,a.out_dir)


if __name__=="__main__":
    main()
