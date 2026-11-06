#!/usr/bin/env python3
"""Find the longest circular interval without ANY annotated gene, per species.

Union of ALL GenBank gene rows (not just human homologs), MITOS2 gene/CDS/
rRNA/tRNA/ncRNA features, and tRNAscan-SE genes. OriH is validation ONLY.
One chosen interval per species per rank; no length threshold. Coordinates
are in the original unrotated FASTA. Missing annotations remain a limitation.
"""
import argparse,csv,re,collections
from pathlib import Path

RANKS=("kingdom","phylum","subphylum","superclass","class","superorder",
       "order","suborder","infraorder","parvorder")
TYPES={"gene","cds","trna","rrna","ncrna","tmrna","snrna","snorna","lncrna","misc_rna"}
NON_GENE=re.compile(r"d[-_ ]?loop|control[-_ ]?region|origin[_ -]?of[_ -]?replication|replication[_ -]?origin|\bori[-_ ]?[hl]\b|\boh\b|\bol\b",re.I)

def norm(x):return re.sub("[^a-z0-9]","",str(x).casefold())
def fasta(path):
    name=None;parts=[]
    with open(path) as fh:
        for line in fh:
            line=line.strip()
            if line.startswith(">"):
                if name is not None:yield name,"".join(parts)
                name=line[1:].strip();parts=[]
            elif line:parts.append(line)
    if name is not None:yield name,"".join(parts)

def number(x):
    try:return int(x)
    except (TypeError,ValueError):return None

def spans(a,b,L):
    """Convert 1-based inclusive circular coordinates to half-open pieces."""
    if a is None or b is None or not (1<=a<=L and 1<=b<=L):return []
    return [(a-1,b)] if b>=a else [(a-1,L),(0,b)]

def union(v):
    result=[]
    for a,b in sorted(v):
        if not result or a>result[-1][1]:result.append([a,b])
        else:result[-1][1]=max(result[-1][1],b)
    return result

def longest_gap(annotations,L):
    occupied=union(annotations)
    if not occupied:return None,0
    gaps=[]
    for i,(_,b) in enumerate(occupied):
        a=occupied[(i+1)%len(occupied)][0]+(L if i==len(occupied)-1 else 0)
        if a>b:gaps.append((b,a))
    return (min(gaps,key=lambda g:(-(g[1]-g[0]),g[0]%L)),len(occupied)) if gaps else (None,len(occupied))

def overlap(a,b,c,d,L):
    """Intervals a:b and c:d may each wrap around the circular origin."""
    def pieces(x,y):
        n=y-x
        if n<0 or n>L:return []
        x%=L
        return [(x,min(x+n,L))]+([(0,x+n-L)] if x+n>L else [])
    return sum(max(0,min(y,v)-max(x,u)) for x,y in pieces(a,b)
               for u,v in pieces(c,d))

def gff_annotations(path,L):
    gene=[];oh=[];count=0
    with open(path,errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):continue
            x=line.rstrip("\n").split("\t")
            if len(x)<9:continue
            type_=x[2].lower()
            attrs={k.casefold():v for item in x[8].split(";") if "=" in item
                   for k,v in [item.split("=",1)]}
            name=next((attrs[k] for k in ("name","gene_id","id") if attrs.get(k)),"")
            pieces=spans(number(x[3]),number(x[4]),L)
            if not pieces:continue
            if type_ in ("origin_of_replication","rep_origin"):
                labels=[re.sub("[^A-Z0-9]","",attrs.get(k,"").upper())
                        for k in ("name","gene_id","id","product")]
                if any(label in {"OH","ORIH","GENEOH","TRANSCRIPTOH","HEAVYSTRANDORIGIN"} for label in labels):
                    oh.append(pieces)
            if type_ in TYPES and not NON_GENE.search(name):
                gene.extend(pieces);count+=1
    return gene,oh,count

def genbank_genes(path):
    found=collections.defaultdict(list)
    with open(path,newline="",encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            k=norm(r.get("Organism") or r.get("organism") or r.get("species") or "")
            if not k:continue
            label=" ".join(str(r.get(x) or "") for x in ("Gene","G_Symbol","G_locus_tag","alt_name","alt_Product"))
            typ=str(r.get("feature") or r.get("Feature") or r.get("type") or "").casefold().strip()
            if NON_GENE.search(label) or (typ and typ not in TYPES):continue
            a,b=number(r.get("start")),number(r.get("end"))
            if a is not None and b is not None:found[k].append((a,b))
    return found

def trnascan_genes(path):
    result=collections.defaultdict(list)
    with open(path,newline="",encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            k=norm(r.get("Organism",""));a=number(r.get("Begin"));b=number(r.get("End"))
            if k and a is not None and b is not None:result[k].append((min(a,b),max(a,b)))
    return result

def gff_paths(path):
    d={}
    with open(path,newline="",encoding="utf-8-sig") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            k=norm(r.get("species",""));v=(r.get("result_gff") or "").strip()
            if k and v:
                if k in d and d[k]!=v:raise ValueError("Conflicting GFF paths: "+k)
                d[k]=v
    return d

def write(path,rows):
    if not rows:raise ValueError("No output rows")
    with open(path,"w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
        writer.writeheader();writer.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta-dir",default="data/unrot_fasta")
    p.add_argument("--genbank-annotations",default="data/kingdom_genes_v2.csv")
    p.add_argument("--trnascan-csv",default="data/consolidated_trnascan-2.csv")
    p.add_argument("--mitos2-summary",default="data/mitos2_species_summary.tsv")
    p.add_argument("--ranks",nargs="+",choices=RANKS,default=list(RANKS))
    p.add_argument("--out-dir",default="longest_all_gene_desert_qc")
    args=p.parse_args();out=Path(args.out_dir);out.mkdir(parents=True,exist_ok=True)
    gb=genbank_genes(args.genbank_annotations);tr=trnascan_genes(args.trnascan_csv)
    gff=gff_paths(args.mitos2_summary)
    cache={};rows=[];summaries=[]
    for rank in args.ranks:
        rankrows=[];seen=set()
        for species,seq in fasta(Path(args.fasta_dir)/(rank+"_chrm.faa")):
            k=norm(species)
            if k=="homosapiens":continue
            if k in seen:raise ValueError(f"{rank}: duplicate species {species}")
            seen.add(k);L=len(seq)
            if (k,L) not in cache:
                file=gff.get(k)
                cache[k,L]=gff_annotations(file,L) if file and Path(file).is_file() else ([],[],0)
            gg,origins,ng=cache[k,L]
            occupied=list(gg);n_gb=n_tr=0
            for a,b in gb.get(k,[]):
                z=spans(a,b,L)
                if z:occupied.extend(z);n_gb+=1
            for a,b in tr.get(k,[]):
                z=spans(a,b,L)
                if z:occupied.extend(z);n_tr+=1
            candidate,nblocks=longest_gap(occupied,L)
            gff_ok=bool(gff.get(k) and Path(gff[k]).is_file())
            oh_size=sum(b-a for group in origins for a,b in group)
            oh_overlap=0;oh_hits=0
            if candidate:
                for group in origins:
                    this=sum(overlap(*candidate,a,b,L) for a,b in group)
                    oh_overlap+=this
                    oh_hits+=int(this>0)
            row=dict(rank=rank,species=species,sequence_length=L,
                n_genbank_gene_features=n_gb,n_mitos2_gene_features=ng,
                n_trnascan_gene_features=n_tr,n_merged_occupied_blocks=nblocks,
                candidate_defined=int(candidate is not None),
                longest_gap_bp=candidate[1]-candidate[0] if candidate else "",
                gap_start_1based=candidate[0]%L+1 if candidate else "",
                gap_end_1based=candidate[1]%L if candidate else "",
                wraps_origin=int(candidate[1]>L) if candidate else "",
                mitos2_gff_available=int(gff_ok),
                n_OriH_features=len(origins) if gff_ok else "",
                OriH_total_annotated_bp=oh_size if gff_ok else "",
                OriH_overlap_bp=oh_overlap if candidate and gff_ok else "",
                OriH_overlaps=int(oh_overlap>0) if candidate and gff_ok else "",
                OriH_fully_contained=int(oh_size>0 and oh_overlap==oh_size) if candidate and gff_ok else "",
                OriH_coverage_fraction=round(oh_overlap/oh_size,6) if candidate and oh_size else "",
                n_OriH_features_overlapping=oh_hits if candidate and gff_ok else "",
                status=("no_annotated_genes" if nblocks==0 else
                        "no_intergenic_gap" if candidate is None else
                        "mitos2_gff_missing" if not gff_ok else
                        "no_OriH_prediction" if not origins else
                        "OriH_fully_contained" if oh_overlap==oh_size else
                        "OriH_partial_overlap" if oh_overlap else "OriH_elsewhere"))
            rows.append(row);rankrows.append(row)
        valid=[r for r in rankrows if r["candidate_defined"]]
        predicted=[r for r in valid if r["mitos2_gff_available"] and r["n_OriH_features"]>0]
        hit=sum(r["OriH_overlaps"]==1 for r in predicted)
        contained=sum(r["OriH_fully_contained"]==1 for r in predicted)
        summaries.append(dict(rank=rank,n_species=len(rankrows),
          n_longest_candidate=len(valid),n_longest_gap_ge300=sum(r["longest_gap_bp"]>=300 for r in valid),
          n_gff_available=sum(r["mitos2_gff_available"] for r in rankrows),
          n_with_OriH_predictions=len(predicted),n_longest_overlaps_OriH=hit,
          n_longest_fully_contains_OriH=contained,
          overlap_rate_among_OriH_predicted=round(100*hit/len(predicted),2) if predicted else "",
          containment_rate_among_OriH_predicted=round(100*contained/len(predicted),2) if predicted else ""))
        print(rank,len(rankrows),"species;",hit,"OriH overlaps",flush=True)
    write(out/"longest_gene_free_interval_by_species.tsv",rows)
    write(out/"all_ranks_summary.tsv",summaries)
    print("Output:",out.resolve())

if __name__=="__main__":main()
