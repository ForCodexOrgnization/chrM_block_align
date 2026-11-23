#!/usr/bin/env python3
"""Kingdom cumulative CR assignment and independent 100% MITOS2 OriH support.

Candidates are defined WITHOUT OriH; the validation uses existing (unrotated)
MITOS2 OriH predictions. Annotation candidates at differing coordinates are
independent: no concatenation of CR1 and CR2 to reach 100%. Reports missing
predictions and unverified coordinate provenance explicitly.
"""
import argparse
import collections
import csv
import json
from pathlib import Path
from audit_genbank_control_region_vocabulary import classify, norm
from control_region_orih_region_level import parse_dloop, parse_misc, read_fasta_manifest, verify_species_coordinate
from compare_genbank_cr_to_longest_orih import (
    evaluate_regions, evaluate_longest, deduplicate, coordinates_to_set)
from compare_pairwise_cr_to_longest_intervals import read_longest,read_mitos2_paths,read_orih_positions

LEVELS=("tier1_explicit","tier1_review","tier2_at_rich",
        "tier2_review","tier3_longest")


def candidate_tier(region):
    """Classify original annotation evidence, including putative D-loop."""
    feature="D-loop" if region["source"]=="D-loop" else ""
    category,tier=classify(region["note"],feature)
    if tier=="tier1":
        return "tier1_explicit"
    if tier=="tier1_review":
        return "tier1_review"
    if tier=="tier2":
        return "tier2_at_rich"
    if tier=="tier2_review":
        return "tier2_review"
    return None


def gather(dloop,misc):
    annotated,invalid=parse_dloop(dloop)
    annotated+=parse_misc(misc)
    # Keep complete and putative source annotations, but exclude motifs and
    # remnants. Reclassify using vocabulary, not old broad regex.
    buckets=collections.defaultdict(list)
    seen={}
    tier_priority={v:i for i,v in enumerate(LEVELS)}
    for reg in annotated:
        tier=candidate_tier(reg)
        if tier is None:
            continue
        identity=(reg["species"],reg["acc"],tuple(sorted(set(reg["intervals"]))))
        if identity not in seen:
            candidate=dict(reg,tier=tier,annotation_sources={reg["source"]})
            seen[identity]=candidate
        else:
            candidate=seen[identity]
            candidate["annotation_sources"].add(reg["source"])
            if tier_priority[tier]<tier_priority[candidate["tier"]]:
                candidate["tier"]=tier
            if candidate["cr_label"]=="CR_unspecified" and reg["cr_label"]!="CR_unspecified":
                candidate["cr_label"]=reg["cr_label"]
    for reg in seen.values():
        buckets[reg["species"]].append(reg)
    return buckets,invalid


def write_tsv(path,records,columns):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter="\t")
        w.writeheader()
        w.writerows(records)


def run(longest_tsv,dloop,misc,mitos_summary,out_dir,manifest_path=None,
        allow_unverified=True):
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    gaps=read_longest(longest_tsv)
    by_species,invalid=gather(dloop,misc)
    gff_paths=read_mitos2_paths(mitos_summary)
    manifest=read_fasta_manifest(manifest_path)
    ori_cache={}
    species_rows=[]
    candidate_rows=[]
    # Kingdom only: the same species must not be duplicated by rank.
    for (rank,species),gap in gaps.items():
        if rank!="kingdom":
            continue
        display=gap.get("species",species)
        length=int(gap.get("sequence_length") or 0)
        current=by_species.get(species,[])
        accessions={r["acc"] for r in current if r["acc"]}
        assembly=verify_species_coordinate(species,accessions,manifest)
        ori_status="missing_gff"
        ori=None
        gff_path=gff_paths.get(species)
        if not length:
            ori_status="invalid_genome_length"
        elif gff_path:
            path=Path(gff_path)
            if not path.is_file():
                ori_status="gff_not_found"
            else:
                k=(str(path.resolve()),length)
                if k not in ori_cache:
                    try:
                        intervals,positions=read_orih_positions(path,length)
                        ori_cache[k]=("predicted" if intervals else "no_OriH",positions)
                    except (ValueError,OSError) as err:
                        ori_cache[k]=("invalid_gff:"+str(err),None)
                ori_status,ori=ori_cache[k]
        has_longest=gap.get("candidate_defined")=="1"
        assignments={tier:[] for tier in LEVELS}
        for r in current:
            assignments[r["tier"]].append(r)
        if has_longest:
            try:
                start,end=int(gap["gap_start_1based"]),int(gap["gap_end_1based"])
                if not 1<=start<=length or not 0<=end<=length:
                    raise ValueError("invalid longest coordinates")
                assignments["tier3_longest"]=[dict(
                    species=species,display=display,acc="",source="longest_gene_free",
                    tier="tier3_longest",cr_label="longest_gene_free",
                    intervals=[(start,end if end else length)],
                    note="longest gene-free interval")]
            except (TypeError,ValueError):
                has_longest=False
        assigned=False
        selected_tier="unresolved"
        selected_interval=""
        for tier in LEVELS:
            if assignments[tier] and not assigned:
                selected_tier=tier
                selected_interval=";".join(
                    ",".join(f"{a}-{b}" for a,b in r["intervals"])
                    for r in assignments[tier])
                assigned=True
        tier_hits={}
        for tier in LEVELS:
            candidates=assignments[tier]
            if ori_status!="predicted" or not allow_unverified and assembly!="identical":
                tier_hits[tier]=None
            elif not candidates:
                tier_hits[tier]=False
            elif tier=="tier3_longest":
                tier_hits[tier],_=evaluate_longest(gap,ori,length)
            elif len(accessions)>1:
                tier_hits[tier]=None
            else:
                tier_hits[tier],_,_=evaluate_regions(candidates,ori,length)
            for i,r in enumerate(candidates):
                if tier=="tier3_longest":
                    full=tier_hits[tier]
                    overlap=""
                    if ori_status=="predicted":
                        _,overlap=evaluate_longest(gap,ori,length)
                elif ori_status=="predicted" and length:
                    full,_,overlap=evaluate_regions([r],ori,length)
                else:
                    full=None;overlap=""
                candidate_rows.append(dict(
                    species=display,organism_key=species,tier=tier,
                    candidate_index=i+1,accession=r["acc"],
                    region_label=r["cr_label"],
                    source=";".join(sorted(r.get("annotation_sources",{r["source"]}))),
                    intervals=";".join(f"{a}-{b}" for a,b in r["intervals"]),
                    note=r.get("note",""),OriH_status=ori_status,
                    OriH_total_bp=len(ori) if ori is not None else "",
                    overlap_bp=overlap,full_containment=(
                        "" if full is None else int(full)),
                    coordinate_verification=assembly))
        cumulative_assignment={}
        cumulative_containment={}
        covered=False
        for tier in LEVELS:
            covered=covered or bool(assignments[tier])
            cumulative_assignment[tier]=int(covered)
            prior=[tier_hits[t] for t in LEVELS[:LEVELS.index(tier)+1]]
            cumulative_containment[tier]=int(any(x is True for x in prior)) if ori_status=="predicted" else ""
        species_rows.append(dict(
            species=display,organism_key=species,genome_bp=length,
            OriH_status=ori_status,OriH_bp=len(ori) if ori is not None else "",
            coordinate_verification=assembly,
            n_genbank_accessions=len(accessions),
            best_annotation_tier=next((t for t in LEVELS[:-1] if assignments[t]),"none"),
            selected_tier=selected_tier,selected_intervals=selected_interval,
            n_candidates=sum(map(len,assignments.values())),
            **{f"has_{t}":int(bool(assignments[t])) for t in LEVELS},
            **{f"contained_{t}":"" if tier_hits[t] is None else int(tier_hits[t]) for t in LEVELS},
            **{f"cumulative_assigned_{t}":cumulative_assignment[t] for t in LEVELS},
            **{f"cumulative_contained_{t}":cumulative_containment[t] for t in LEVELS}))
    if not species_rows:
        raise ValueError("No kingdom species in longest-gap TSV")
    summary=[]
    count=len(species_rows)
    predicted=[r for r in species_rows if r["OriH_status"]=="predicted"]
    n_pred=len(predicted)
    previous_assign=previous_cont=0
    for tier in LEVELS:
        c=f"cumulative_assigned_{tier}"
        h=f"cumulative_contained_{tier}"
        n_assign=sum(r[c]==1 for r in species_rows)
        n_cont=sum(r[h]==1 for r in predicted)
        summary.append(dict(
            tier=tier,n_all_species=count,n_orih_predicted=n_pred,
            n_cumulative_assigned=n_assign,n_newly_assigned=n_assign-previous_assign,
            assignment_pct=round(100*n_assign/count,3),
            n_cumulative_orih_full=n_cont,n_new_orih_full=n_cont-previous_cont,
            containment_pct_of_predicted=round(100*n_cont/n_pred,3) if n_pred else "",
            n_unresolved=count-n_assign))
        previous_assign=n_assign;previous_cont=n_cont
    write_tsv(out/"kingdom_candidates.tsv",candidate_rows,
              list(candidate_rows[0]) if candidate_rows else ["species"])
    write_tsv(out/"kingdom_species_tier_recovery.tsv",species_rows,list(species_rows[0]))
    write_tsv(out/"kingdom_cumulative_recovery.tsv",summary,list(summary[0]))
    report=dict(n_kingdom_species=count,n_orih_predicted=n_pred,
                n_invalid_split_dloop=len(invalid),allow_unverified=allow_unverified,
                cumulative_recovery=summary,
                caveat="Positional OriH containment requires assembly-matched identical FASTAs. Without manifest, exploratory results are provisional and must not be called independent biological validation.")
    (out/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    return report


def plot_summary(rows,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels=["Tier 1 explicit","+ Tier 1 review","+ Tier 2 AT-rich",
            "+ Tier 2 review","+ Longest gap"]
    x=list(range(len(labels)))
    fig,ax=plt.subplots(figsize=(10,5))
    a=[r["assignment_pct"] for r in rows]
    b=[float(r["containment_pct_of_predicted"] or 0) for r in rows]
    ax.plot(x,a,"o-",label="Candidate assigned | All species")
    ax.plot(x,b,"s-",label="100% OriH containment | OriH-predicted species")
    ax.set_ylim(0,100);ax.set_yticks(range(0,101,20))
    ax.set_yticklabels([f"{v}%" for v in range(0,101,20)])
    ax.set_xticks(x);ax.set_xticklabels(labels,rotation=20,ha="right")
    ax.set_ylabel("Kingdom species (%)")
    ax.set_title("Cumulative mitochondrial control-region candidate recovery")
    ax.legend();ax.grid(axis="y",alpha=.25)
    fig.tight_layout()
    fig.savefig(Path(out)/"kingdom_cumulative_recovery.png",dpi=220)
    fig.savefig(Path(out)/"kingdom_cumulative_recovery.pdf")
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--longest-tsv",default="longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv")
    p.add_argument("--genbank",default="data/kingdom_genbank_annotations_split_new.csv")
    p.add_argument("--misc",default="data/genbank_misc_ftr.csv")
    p.add_argument("--mitos-summary",default="data/mitos2_species_summary.tsv")
    p.add_argument("--out-dir",default="control_region_vocabulary_audit/kingdom_recovery")
    p.add_argument("--assembly-manifest")
    p.add_argument("--require-verified-coordinates",action="store_true")
    a=p.parse_args()
    report=run(a.longest_tsv,a.genbank,a.misc,a.mitos_summary,a.out_dir,
               a.assembly_manifest,not a.require_verified_coordinates)
    plot_summary(report["cumulative_recovery"],a.out_dir)


if __name__=="__main__":
    main()
