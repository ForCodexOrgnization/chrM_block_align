#!/usr/bin/env python3
"""Audit Kingdom-wide GenBank mitochondrial control-region annotation vocabulary.

Read both current GenBank CSV files, retain original annotations and separate
whole control regions, plausible alternative names, motifs and remnants.
Outputs evidence-based, annotation-only tiers; does NOT claim OriH validation.
"""
import argparse
import collections
import csv
import json
import re
from pathlib import Path

def norm(s):
    return re.sub(r"[^a-z0-9]", "", (s or "").casefold())

def match(pat, s):
    return bool(re.search(pat, s, re.I))

# Order matters: subfeatures and remnants must not become full CR candidates.
REMNANT = r"\b(?:remnant|relict|pseudo[- ]?control|pseudoregion|degenerate|partial duplication)\b"
MOTIF = (r"\b(?:CSB[- _]?[0-9IVX]+|conserved sequence block|"
         r"coreTAS|TAS|termination.associated.sequence|SSB[- _]?[0-9IVX]+|"
         r"CoRe|control region conserved region|control region C.rich sequence|"
         r"C.rich sequence|mtTERM binding site)\b")
EXPLICIT = r"\b(?:control|contral|contorl|conrtol|cntrol|conrol|conrol)[-_ ]+region\b|\bD[-_ ]?loop\b"
AT_RICH = (r"\b(?:A\s*[+&/\-]\s*T|A\s*T|AT)[-_ ]*rich\b|"
           r"\bAT[-_ ]?(?:loop|region)\b|\bAT_rich\b|"
           r"\bA[+&/]T[-_ ]?region\b|\bA[-_ ]T[-_ ]rich\b|\bA[-_ ]T[-_ ]region\b")
NONCODING = r"\b(?:major|long|large)[-_ ]+(?:non[-_ ]?coding|intergenic)(?:[-_ ]+region)?\b"
ORIGIN = (r"\b(?:origin of replication|replication origin|"
          r"heavy strand origin|Ori[-_ ]?[HL]|OH|OL)\b")
BARE_CR = r"(?<![A-Z0-9])CR[-_ ]?(?:1|2|I|II)(?![A-Z0-9])"
TYPO = r"\b(?:contral|contorl|conrtol|cntrol)[-_ ]+region\b"


def classify(note, feature="", marked=False):
    """Return hierarchical category and normalized interpretation."""
    txt=(note or "").strip()
    if match(REMNANT,txt):
        return "remnant","not_full_control_region"
    if feature.casefold()=="d-loop":
        return ("putative_cr" if match(r"\\b(?:putative|possible|probable|predicted)\\b",txt)
                else "explicit_cr"),"tier1"
    # Treat a feature as a subfeature only when its primary description is
    # itself a short motif; a full control region note may mention CSB/TAS.
    first_clause=re.split(r"[;|]",txt,1)[0].strip()
    if match(MOTIF,first_clause) and not (match(EXPLICIT,first_clause) and
            match(r"\b(?:contains|including|with)\b",first_clause)):
        return "subfeature","not_full_control_region"
    if match(EXPLICIT,txt):
        return ("putative_cr" if match(r"\b(?:putative|possible|probable|predicted)\b",txt)
                else "explicit_cr"),"tier1"
    if match(BARE_CR,txt) and not match(r"\b(?:chloroplast|plastid)\b",txt):
        return "numbered_cr_abbreviation","tier1_review"
    if marked:
        return "source_dloop_like","tier1_review"
    if match(AT_RICH,txt) or match(r"\(\s*A\s*\+\s*T\s*\)[-_ ]*rich",txt):
        return "at_rich_alternative","tier2"
    if match(NONCODING,txt):
        return "long_noncoding_alternative","tier2_review"
    if match(ORIGIN,txt):
        return "replication_origin_only","origin_not_cr"
    return "unrelated","none"


def write(path,columns,rows):
    with path.open("w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=columns,delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def audit(genbank,misc,out_dir,kingdom_tsv=None):
    out=Path(out_dir);out.mkdir(parents=True,exist_ok=True)
    vocabulary=collections.defaultdict(lambda:dict(n=0,species=set(),accessions=set(),examples=set()))
    species_tiers=collections.defaultdict(set)
    candidates=[]
    seen=set()
    kingdom=None
    if kingdom_tsv:
        kingdom=set()
        with open(kingdom_tsv,newline="",encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh,delimiter="\t"):
                if (row.get("rank") or "").casefold()=="kingdom":
                    kingdom.add(norm(row.get("species","")))
    def record(source,row,label,feature=""):
        sp=row.get("Organism") or row.get("species") or ""
        nk=norm(sp)
        if not nk or (kingdom is not None and nk not in kingdom):
            return
        acc=(row.get("acc") or "").strip()
        marked=str(row.get("dloop_like","")).lower()=="true"
        category,tier=classify(label,feature,marked)
        raw=(row.get("Note") or "").strip()
        key=(source,category,raw)
        v=vocabulary[key]
        v["n"]+=1;v["species"].add(nk);v["accessions"].add(acc)
        if label and len(v["examples"])<4:
            v["examples"].add(label[:180])
        if tier in {"tier1","tier1_review","tier2","tier2_review"}:
            species_tiers[nk].add(tier)
            signature=(source,acc,sp,category,str(row.get("start","")),
                       str(row.get("end","")),str(row.get("location","")))
            if signature not in seen:
                seen.add(signature)
                candidates.append(dict(species=sp,accession=acc,source=source,
                    category=category,tier=tier,note=raw,
                    start=row.get("start",""),end=row.get("end",""),
                    location=row.get("location",""),
                    feature_index=row.get("feature_index",""),
                    dloop_like=row.get("dloop_like","")))
    with open(genbank,newline="",encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            if (row.get("class") or "").casefold()!="d-loop":
                continue
            label=" ".join(row.get(x,"") for x in (
                "Gene","G_Symbol","alt_name","alt_synonym","alt_Product","Note"))
            record("GenBank_D-loop",row,label,row.get("class",""))
    with open(misc,newline="",encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            if (row.get("misc_kind") or "").casefold()!="other":
                continue
            label=" ".join(row.get(x,"") for x in (
                "Note","Gene","G_Symbol","Product","other_qualifiers"))
            record("GenBank_misc_feature",row,label)
    vocab_rows=[]
    for (source,category,note),v in vocabulary.items():
        vocab_rows.append(dict(source=source,category=category,raw_note=note,
            annotation_rows=v["n"],n_species=len(v["species"]),
            n_accessions=len(v["accessions"]),example_text=" | ".join(sorted(v["examples"]))))
    vocab_rows.sort(key=lambda r:(r["source"],r["category"],-r["annotation_rows"],r["raw_note"]))
    write(out/"annotation_vocabulary.tsv",[
        "source","category","raw_note","annotation_rows","n_species","n_accessions",
        "example_text"],vocab_rows)
    write(out/"candidate_annotations.tsv",[
        "species","accession","source","category","tier","note",
        "start","end","location","feature_index","dloop_like"],candidates)
    rows=[]
    if kingdom is None:
        keys=set(species_tiers)
    else:
        keys=kingdom
    for species in sorted(keys):
        tiers=species_tiers.get(species,set())
        rows.append(dict(organism_key=species,
            has_tier1=int("tier1" in tiers),has_tier1_review=int("tier1_review" in tiers),
            has_tier2=int("tier2" in tiers),has_tier2_review=int("tier2_review" in tiers),
            best_annotation_tier=(
                "tier1" if "tier1" in tiers else
                "tier1_review" if "tier1_review" in tiers else
                "tier2" if "tier2" in tiers else
                "tier2_review" if "tier2_review" in tiers else "none")))
    write(out/"species_annotation_tiers.tsv",[
        "organism_key","has_tier1","has_tier1_review","has_tier2",
        "has_tier2_review","best_annotation_tier"],rows)
    histogram=collections.Counter(row["best_annotation_tier"] for row in rows)
    counts=collections.Counter(row["category"] for row in candidates)
    result=dict(n_species_scoped=len(rows),n_candidate_annotation_rows=len(candidates),
                species_tier_counts=dict(sorted(histogram.items())),
                candidate_category_counts=dict(sorted(counts.items())),
                typo_contral_region_rows=sum(v["n"] for (src,cat,note),v in vocabulary.items()
                    if match(TYPO,note)),
                validation_warning="Annotation nomenclature tiers are not biological proof of control-region function; OriH comparison requires matched assemblies.")
    (out/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genbank",default="data/kingdom_genbank_annotations_split_new.csv")
    p.add_argument("--misc",default="data/genbank_misc_ftr.csv")
    p.add_argument("--kingdom-tsv",default="longest_all_gene_desert_qc/longest_gene_free_interval_by_species.tsv")
    p.add_argument("--all-genbank-species",action="store_true",help="Disable Kingdom cohort filter")
    p.add_argument("--out-dir",default="control_region_vocabulary_audit")
    a=p.parse_args()
    audit(a.genbank,a.misc,a.out_dir,None if a.all_genbank_species else a.kingdom_tsv)


if __name__=="__main__":
    main()
