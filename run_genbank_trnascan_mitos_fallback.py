#!/usr/bin/env python3
"""GenBank-first 37 genes; tRNAscan tRNA and MITOS2 CDS/rRNA fallback pipeline.

Input annotation CSV requires: acc, Organism, Gene, start, end, strand.
Alias table: the 38-column genefiltfinal5_2026_09_24.csv (columns are
canonical gene identities; cells are alternate GenBank names).

Gene-only runs the existing baseline stitcher; pad50/pad100 run the
existing padded stitcher. Both receive the same merged annotation source,
with tRNAscan fallback for missing tRNAs and MITOS2 fallback for missing CDS/rRNA.
Ambiguous GenBank loci are left unresolved, not silently overridden.
Human GenBank coordinates always take priority.
"""
import argparse,csv,collections,logging,re
from pathlib import Path
from types import SimpleNamespace
import chrM_block_align as base
import chrM_block_align_padded as padded
from human_genbank import load_human_genbank
import resolve_gene_order_annotations as go

LOG=logging.getLogger("genbank_trnascan_mitos_fallback")
COLMAP={"nd1":"ND1","nd2":"ND2","nd3":"ND3","nd4":"ND4","nd4l":"ND4L",
        "nd5":"ND5","nd6":"ND6","cox1":"CO1","cox2":"CO2","cox3":"CO3",
        "cyb":"CYB","atp6":"ATP6","atp8":"ATP8","rnr1":"RNR1","rnr2":"RNR2"}
for _aa in "A R N D C Q E G H I L K M F P S T W Y V".split():
    COLMAP["MT_T"+_aa]="TRN"+_aa
for _s in ("L1","L2","S1","S2"):COLMAP["MT_T"+_s]="TRN"+_s

def normalize_name(value):
    return re.sub(r"[^a-z0-9]","",str(value).casefold())

def load_aliases(path):
    with open(path,newline="",encoding="utf-8-sig") as inp:
        rows=csv.DictReader(inp)
        alias=collections.defaultdict(set)
        for col,canonical in COLMAP.items():
            alias[normalize_name(col)].add(canonical)
            alias[normalize_name(canonical)].add(canonical)
        for row in rows:
            for col,canonical in COLMAP.items():
                raw=(row.get(col) or "").strip()
                if raw:alias[normalize_name(raw)].add(canonical)
    return {k:frozenset(v) for k,v in alias.items()}

def load_genbank_csv(path,aliases,reference):
    by_species=collections.defaultdict(list);unmapped=[];invalid=[];resolution=[]
    groups=collections.defaultdict(list)
    with open(path,newline="",encoding="utf-8-sig") as inp:
        reader=csv.DictReader(inp)
        if not {"start","end","strand"}.issubset(reader.fieldnames or []):
            raise ValueError("GenBank annotations require start,end,strand")
        if not ({"Organism","organism","species"} & set(reader.fieldnames or [])):
            raise ValueError("GenBank annotations require Organism/organism/species")
        for idx,row in enumerate(reader,2):
            species=row.get("Organism") or row.get("organism") or row.get("species") or ""
            species_key=base.norm(species)
            if not species_key:invalid.append((idx,"missing species"));continue
            hits=go.candidates(row,aliases)
            if not hits:
                unmapped.append((idx,species,row.get("Gene","")));continue
            try:start,end=int(row["start"]),int(row["end"])
            except (ValueError,TypeError):
                invalid.append((idx,"invalid coordinates"));continue
            strand=(row.get("strand") or "").strip()
            if strand not in ("+","-"):
                invalid.append((idx,"invalid strand"));continue
            groups[species_key].append((hits,start,end,strand,row.get("acc",""),idx,species))
    for species_key,entries in groups.items():
        accessions={e[4] for e in entries if e[4]}
        if len(accessions)>1:
            invalid.extend((e[5],"multiple_accessions") for e in entries)
            continue
        entries.sort(key=lambda e:e[1])
        choices=[(e[0],e[3],e[5]) for e in entries]
        resolved,details=go.resolve(choices,reference)
        for detail in details:
            e=entries[detail["index"]]
            resolution.append((species_key,e[6],e[5],detail["candidates"],
                               detail["status"],detail["assigned"],
                               detail.get("supported",""),detail.get("n_best_assignments","")))
        for e,gene in zip(entries,resolved):
            if gene is not None:
                by_species[species_key].append((gene,e[1],e[2],e[3],e[4],e[5]))
    return by_species,unmapped,invalid,resolution

def build_feature_sets(rank,fa,gb,trnas,mitos_paths,strict=False):
    """GenBank for all 37; tRNAscan for missing tRNAs; MITOS2 for missing CDS/rRNA.

    Missing means no unambiguous valid GenBank locus. An ambiguous GenBank
    locus is withheld rather than superseded by an alternate annotation.
    """
    merged={};audit=[];ambiguous=[]
    for name,seq in base.fasta(fa):
        key=base.norm(name);L=len(seq)
        if key=="homosapiens":continue
        gb_features=[]
        for gene,start,end,strand,accession,idx in gb.get(key,[]):
            if gene not in base.EXPECTED:continue
            if start<1 or end<1 or start>L or end>L:
                audit.append((rank,name,gene,"genbank","invalid_coordinates",accession,idx))
                continue
            a=start-1;b=end if end>=start else end+L
            if b<=a or b-a>L:
                audit.append((rank,name,gene,"genbank","invalid_interval",accession,idx))
                continue
            gb_features.append(base.Feature(gene,a,b,strand))
        chosen_gb,amb_gb=base.pick_features(gb_features,L)
        chosen_trnas,amb_trnas=base.pick_features(trnas.get(key,[]),L)
        chosen_trnas={g:f for g,f in chosen_trnas.items() if g in base.TRNA}
        mitos_gff=mitos_paths.get(key)
        mitos_features=(base.read_gff(mitos_gff,L)
                        if mitos_gff and Path(mitos_gff).is_file() else [])
        chosen_mitos,amb_mitos=base.pick_features(mitos_features,L)
        chosen_mitos={g:f for g,f in chosen_mitos.items() if g in (base.CDS|base.RRNA)}
        chosen=dict(chosen_gb)
        sources={g:"genbank" for g in chosen}
        for gene in sorted(base.EXPECTED):
            if gene in chosen or gene in amb_gb:
                continue
            if gene in base.TRNA and gene in chosen_trnas:
                chosen[gene]=chosen_trnas[gene]
                sources[gene]="trnascan_fallback"
            elif gene in (base.CDS|base.RRNA) and gene in chosen_mitos:
                chosen[gene]=chosen_mitos[gene]
                sources[gene]="mitos2_fallback"
        for gene in sorted(base.EXPECTED):
            if gene in sources:src=sources[gene];status="selected"
            elif gene in amb_gb:src="ambiguous_genbank";status="excluded"
            elif gene in base.TRNA and gene in amb_trnas:
                src="ambiguous_trnascan";status="excluded"
            elif gene in (base.CDS|base.RRNA) and gene in amb_mitos:
                src="ambiguous_mitos2";status="excluded"
            else:src="missing";status="excluded"
            audit.append((rank,name,gene,src,status,"",""))
        for gene,n in amb_gb.items():ambiguous.append((rank,name,gene,n,"genbank"))
        for gene,n in amb_trnas.items():
            if gene not in chosen_gb:ambiguous.append((rank,name,gene,n,"trnascan"))
        for gene,n in amb_mitos.items():
            if gene not in chosen_gb:ambiguous.append((rank,name,gene,n,"mitos2"))
        merged[key]=chosen
    n_gb=sum(row[3]=="genbank" and row[4]=="selected" for row in audit)
    counts=collections.Counter(row[3] for row in audit if row[4]=="selected")
    LOG.info("%s: %d GenBank features; tRNAscan fallback %d; MITOS2 fallback %d",
             rank,n_gb,counts["trnascan_fallback"],counts["mitos2_fallback"])
    if n_gb==0:raise ValueError(f"No GenBank genes matched {rank}; check annotations")
    if strict and any(r[4]=="excluded" for r in audit):
        raise ValueError("Strict annotation validation failed; see annotation_source_audit.tsv")
    return merged,audit,ambiguous

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genbank-annotations",required=True,
                   help="Export the kingdom_genes_v2 sheet as CSV with Organism/Gene/start/end/strand")
    p.add_argument("--gene-aliases",required=True,
                   help="genefiltfinal5_2026_09_24.csv synonym table")
    p.add_argument("--fasta-dir",default="data/unrot_fasta")
    p.add_argument("--trnascan-csv",default="data/consolidated_trnascan-2.csv")
    p.add_argument("--mitos2-summary",default="data/mitos2_species_summary.tsv")
    p.add_argument("--human-fasta",default="data/Humans_mtdna.fasta")
    p.add_argument("--human-genbank",default="data/sequence.gb")
    p.add_argument("--output-dir",default="block_results_genbank_trnascan_mitos_fallback_v1")
    p.add_argument("--ranks",nargs="+",default=["class"])
    p.add_argument("--approaches",nargs="+",choices=["gene_only","pad50","pad100"],
                   default=["gene_only","pad50","pad100"])
    p.add_argument("--resume-incomplete",action="store_true",help="Reuse validated existing MAFFT caches in an incomplete rank; never overwrite completed stitched output")
    p.add_argument("--threads",type=int,default=4)
    p.add_argument("--mafft",default="mafft")
    p.add_argument("--max-flank-bp",type=int,default=2500)
    p.add_argument("--flank-batch-size",type=int,default=100)
    a=p.parse_args()
    logging.basicConfig(level=logging.INFO,format="%(asctime)s %(levelname)s %(message)s")
    outroot=Path(a.output_dir).resolve()
    if outroot in {Path(x).resolve() for x in ("block_results_v6","block_results_padded_v1",
            "block_results_genbank_trnascan_v1","block_results_genbank_mitos_v1")}:
        p.error("Do not use existing output folders")
    aliases=load_aliases(a.gene_aliases)
    human=list(base.fasta(a.human_fasta))[0][1]
    reference_features=load_human_genbank(a.human_genbank,human,base.Feature,base.EXPECTED)
    reference=[(g,reference_features[g].strand) for g in sorted(reference_features,key=lambda g:reference_features[g].start)]
    gb,unmapped,invalid,resolution=load_genbank_csv(a.genbank_annotations,aliases,reference)
    if not gb:raise ValueError("No recognized GenBank annotation rows. Check aliases and Organism names.")
    trnas=base.read_trnas(a.trnascan_csv)
    with open(a.mitos2_summary,newline="",encoding="utf-8-sig") as fh:
        mitos={base.norm(r["species"]):r["result_gff"] for r in csv.DictReader(fh,delimiter="\t")
               if r.get("species") and r.get("result_gff")}
    orig_base_read=base.read_gff
    orig_padded_read=padded.read_gff
    orig_base_trnas=base.read_trnas
    orig_padded_trnas=padded.read_trnas
    for rank in a.ranks:
        fa=Path(a.fasta_dir)/(rank+"_chrm.faa")
        if not fa.exists():raise FileNotFoundError(fa)
        selected,audit,ambiguous=build_feature_sets(rank,fa,gb,trnas,mitos)
        audit_dir=outroot/"annotation_qc"/rank;audit_dir.mkdir(parents=True,exist_ok=True)
        base.write_tsv(audit_dir/"annotation_source_audit.tsv",audit,
                       ["rank","species","gene","source","status","accession","row"])
        base.write_tsv(audit_dir/"ambiguous_gene_annotations.tsv",ambiguous,
                       ["rank","species","gene","n_loci","source"])
        base.write_tsv(audit_dir/"unmapped_gene_names.tsv",unmapped,
                       ["row","species","raw_gene"])
        base.write_tsv(audit_dir/"invalid_rows.tsv",invalid,["row","reason"])
        base.write_tsv(audit_dir/"gene_order_resolutions.tsv",resolution,
                       ["organism_key","species","row","candidates","status","assigned","supported","n_best_assignments"])
        # The existing run() code checks that each nonhuman species has a
        # 'GFF' path. Use the actual FASTA as an existing path and substitute
        # the merged in-memory feature map at read_gff(). The original modules
        # and their on-disk files remain unchanged.
        # Put a species key in a lightweight proxy path and intercept its
        # content; selected annotations are already prepared above.
        proxy_dir=audit_dir/"feature_proxies";proxy_dir.mkdir(exist_ok=True)
        proxy_paths={}
        for species in selected:
            file=proxy_dir/(species+".features")
            file.touch(exist_ok=True)
            proxy_paths[species]=str(file)
        def read_selected(path,length):
            key=Path(path).stem
            return list(selected[key].values())
        base.read_gff=read_selected
        padded.read_gff=read_selected
        empty=collections.defaultdict(list)
        base.read_trnas=lambda path:empty
        padded.read_trnas=lambda path:empty
        for approach in a.approaches:
            target=outroot/approach/rank
            if target.exists() and any(target.iterdir()):
                stitched=target/(rank+"_stitched_aln.fasta")
                if not a.resume_incomplete or stitched.exists():
                    raise FileExistsError(f"Refusing to overwrite existing results: {target}; use --resume-incomplete only when no stitched alignment exists")
                LOG.info("Resuming incomplete %s/%s using validated MAFFT caches",approach,rank)
            args=SimpleNamespace(output_dir=str(outroot/approach),
                 human_fasta=a.human_fasta,human_genbank=a.human_genbank,
                 threads=a.threads,mafft=a.mafft,
                 max_flank_sequences=300,max_flank_bp=a.max_flank_bp,
                 flank_batch_size=a.flank_batch_size,keep_intermediate=True,
                 padding=50 if approach=="pad50" else 100)
            if approach=="gene_only":
                base.run(rank,fa,proxy_paths,empty,args)
            else:
                padded.run(rank,fa,proxy_paths,empty,args)
        base.read_gff=orig_base_read;padded.read_gff=orig_padded_read
        base.read_trnas=orig_base_trnas;padded.read_trnas=orig_padded_trnas

if __name__=="__main__":main()
