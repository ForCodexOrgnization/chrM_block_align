#!/usr/bin/env python3
"""Isolated GenBank-first + MITOS2-fallback chrM block pipeline.

Input annotation CSV requires: acc, Organism, Gene, start, end, strand.
Alias table: the 38-column genefiltfinal5_2026_09_24.csv (columns are
canonical gene identities; cells are alternate GenBank names).

Gene-only runs the existing baseline stitcher; pad50/pad100 run the
existing padded stitcher. Both receive the same merged annotation source,
which never reads tRNAscan. Human GenBank coordinates always take priority.
"""
import argparse,csv,collections,logging,re
from pathlib import Path
from types import SimpleNamespace
import chrM_block_align as base
import chrM_block_align_padded as padded

LOG=logging.getLogger("genbank_mitos")
COLMAP={"nd1":"ND1","nd2":"ND2","nd3":"ND3","nd4":"ND4","nd4l":"ND4L",
        "nd5":"ND5","nd6":"ND6","cox1":"CO1","cox2":"CO2","cox3":"CO3",
        "cyb":"CYB","atp6":"ATP6","atp8":"ATP8","rnr1":"RNR1","rnr2":"RNR2"}
for _aa in "A R N D C Q E G H I L K M F P S T W Y V".split():
    COLMAP["MT_T"+_aa]="TRN"+_aa
for _s in ("L1","L2","S1","S2"):COLMAP["MT_T"+_s]="TRN"+_s

def normalize_name(value):
    return re.sub(r"[^a-z0-9]","",str(value).casefold())

def load_aliases(path):
    with open(path,newline="") as inp:
        rows=csv.DictReader(inp)
        unknown=set(rows.fieldnames or [])-set(COLMAP)-{"D-loop"}
        if unknown:raise ValueError(f"Unknown alias columns: {sorted(unknown)}")
        alias=collections.defaultdict(set)
        for col,canonical in COLMAP.items():
            alias[normalize_name(col)].add(canonical)
            alias[normalize_name(canonical)].add(canonical)
        for row in rows:
            for col,canonical in COLMAP.items():
                raw=(row.get(col) or "").strip()
                if raw:alias[normalize_name(raw)].add(canonical)
    conflicts={k:sorted(v) for k,v in alias.items() if len(v)>1}
    if conflicts:raise ValueError(f"Ambiguous gene aliases: {conflicts}")
    return {k:next(iter(v)) for k,v in alias.items()}

def load_genbank_csv(path,aliases):
    by_species=collections.defaultdict(list);unmapped=[];invalid=[]
    with open(path,newline="") as inp:
        reader=csv.DictReader(inp)
        required={"Gene","start","end","strand"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("GenBank annotation CSV must contain Gene,start,end,strand "
                             "(plus Organism or an accession-to-species mapping). "
                             "The uploaded synonym CSV is NOT the coordinate file.")
        if not ({"Organism","organism","species"} & set(reader.fieldnames or [])):
            raise ValueError("GenBank coordinate CSV needs Organism/organism/species to match rank FASTA headers.")
        for idx,row in enumerate(reader,2):
            species=row.get("Organism") or row.get("organism") or row.get("species") or ""
            key=base.norm(species)
            if not key:invalid.append((idx,"missing species"));continue
            raw=(row.get("Gene") or "").strip()
            gene=aliases.get(normalize_name(raw))
            if gene is None:
                unmapped.append((idx,species,raw));continue
            try:
                start=int(float(row["start"]));end=int(float(row["end"]))
            except (TypeError,ValueError):
                invalid.append((idx,"invalid coordinates"));continue
            strand=(row.get("strand") or "").strip()
            if strand not in ("+","-"):
                invalid.append((idx,"invalid strand"));continue
            # GenBank coordinate export is 1-based inclusive; circular coordinates
            # are normalized to forward 0-based intervals per species below.
            by_species[key].append((gene,start,end,strand,row.get("acc",""),idx))
    return by_species,unmapped,invalid

def build_feature_sets(rank,fa,gb,mitos_paths,strict=False):
    merged={};audit=[];ambiguous=[]
    for name,seq in base.fasta(fa):
        key=base.norm(name);L=len(seq)
        if key=="homosapiens":continue
        genbank=[]
        for gene,start,end,strand,accession,idx in gb.get(key,[]):
            if start<1 or end<1 or start>L or end>L:
                audit.append((rank,name,gene,"genbank","invalid_coordinates",accession,idx))
                continue
            a=start-1;b=end if end>=start else end+L
            if b<=a or b-a>L:
                audit.append((rank,name,gene,"genbank","invalid_interval",accession,idx))
                continue
            genbank.append(base.Feature(gene,a,b,strand))
        chosen_gb,amb_gb=base.pick_features(genbank,L)
        mitos=[]
        gff=mitos_paths.get(key)
        if gff and Path(gff).exists():
            mitos=base.read_gff(gff,L)
        chosen_mitos,amb_mitos=base.pick_features(mitos,L)
        chosen=dict(chosen_gb)
        for gene,feature in chosen_mitos.items():
            if gene not in chosen and gene not in amb_gb:
                chosen[gene]=feature
        for gene in base.EXPECTED:
            if gene in chosen_gb:src="genbank"
            elif gene in chosen:src="mitos2_fallback"
            elif gene in amb_gb:src="ambiguous_genbank"
            elif gene in amb_mitos:src="ambiguous_mitos2"
            else:src="missing"
            audit.append((rank,name,gene,src,"selected" if gene in chosen else "excluded","", ""))
        for gene,n in amb_gb.items():ambiguous.append((rank,name,gene,n,"genbank"))
        for gene,n in amb_mitos.items():ambiguous.append((rank,name,gene,n,"mitos2"))
        merged[key]=chosen
    matched=sum(1 for row in audit if row[3]=="genbank")
    LOG.info("%s: %d genes assigned from GenBank; %d species with annotation records",rank,matched,len(gb))
    if matched==0:raise ValueError(f"No GenBank genes matched rank {rank}. Check species names and annotation coordinate file.")
    if strict and any(r[3] in ("missing","ambiguous_genbank") for r in audit):
        raise ValueError("Strict annotation validation failed; see annotation_source_audit.tsv")
    return merged,audit,ambiguous

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genbank-annotations",required=True,
                   help="Export the kingdom_genes_v2 sheet as CSV with Organism/Gene/start/end/strand")
    p.add_argument("--gene-aliases",required=True,
                   help="genefiltfinal5_2026_09_24.csv synonym table")
    p.add_argument("--fasta-dir",default="data/unrot_fasta")
    p.add_argument("--mitos2-summary",default="data/mitos2_species_summary.tsv")
    p.add_argument("--human-fasta",default="data/Humans_mtdna.fasta")
    p.add_argument("--human-genbank",default="data/sequence.gb")
    p.add_argument("--output-dir",default="block_results_genbank_mitos_v1")
    p.add_argument("--ranks",nargs="+",default=["class"])
    p.add_argument("--approaches",nargs="+",choices=["gene_only","pad50","pad100"],
                   default=["gene_only","pad50","pad100"])
    p.add_argument("--threads",type=int,default=4)
    p.add_argument("--mafft",default="mafft")
    p.add_argument("--max-flank-bp",type=int,default=2500)
    p.add_argument("--flank-batch-size",type=int,default=100)
    a=p.parse_args()
    logging.basicConfig(level=logging.INFO,format="%(asctime)s %(levelname)s %(message)s")
    outroot=Path(a.output_dir).resolve()
    if outroot in {Path("block_results_v6").resolve(),Path("block_results_padded_v1").resolve()}:
        p.error("Do not use existing output folders")
    aliases=load_aliases(a.gene_aliases)
    gb,unmapped,invalid=load_genbank_csv(a.genbank_annotations,aliases)
    if not gb:raise ValueError("No recognized GenBank annotation rows. Check aliases and Organism names.")
    with open(a.mitos2_summary,newline="") as h:
        mitos={base.norm(r["species"]):r["result_gff"] for r in csv.DictReader(h,delimiter="\t") if r.get("result_gff")}
    orig_base_read=base.read_gff
    orig_padded_read=padded.read_gff
    orig_base_trnas=base.read_trnas
    orig_padded_trnas=padded.read_trnas
    for rank in a.ranks:
        fa=Path(a.fasta_dir)/(rank+"_chrm.faa")
        if not fa.exists():raise FileNotFoundError(fa)
        selected,audit,ambiguous=build_feature_sets(rank,fa,gb,mitos)
        audit_dir=outroot/"annotation_qc"/rank;audit_dir.mkdir(parents=True,exist_ok=True)
        base.write_tsv(audit_dir/"annotation_source_audit.tsv",audit,
                       ["rank","species","gene","source","status","accession","row"])
        base.write_tsv(audit_dir/"ambiguous_gene_annotations.tsv",ambiguous,
                       ["rank","species","gene","n_loci","source"])
        base.write_tsv(audit_dir/"unmapped_gene_names.tsv",unmapped,
                       ["row","species","raw_gene"])
        base.write_tsv(audit_dir/"invalid_rows.tsv",invalid,["row","reason"])
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
                raise FileExistsError(f"Refusing to overwrite existing results: {target}")
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
