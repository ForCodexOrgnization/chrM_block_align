#!/usr/bin/env python3
"""Full-gene mitochondrial MAFFT + deduplicated human-reference stitching."""
import argparse
import collections
import csv
import json
import logging
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

COMP=str.maketrans("ACGTNacgtn-","TGCANtgcan-")
CDS=set("ND1 ND2 ND3 ND4 ND4L ND5 ND6 CO1 CO2 CO3 ATP6 ATP8 CYB".split())
RRNA={"RNR1","RNR2"}
TRNA=set("TRNA TRNR TRNN TRND TRNC TRNQ TRNE TRNG TRNH TRNI TRNK TRNM TRNF TRNP TRNT TRNW TRNY TRNV TRNL1 TRNL2 TRNS1 TRNS2".split())
EXPECTED=CDS|RRNA|TRNA
AA=dict(zip("Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val".lower().split(),
            list("ARNDCQEGHILKMFPSTWYV")))
LOG=logging.getLogger("chrM")
def norm(s):return re.sub("[^a-z0-9]","",s.lower())
def rc(s):return s.translate(COMP)[::-1]
def fasta(path):
    name=None;chunks=[]
    with open(path) as inp:
        for line in inp:
            s=line.strip()
            if not s:continue
            if s.startswith(">"):
                if name is not None:yield name,"".join(chunks).upper()
                name=s[1:].strip();chunks=[]
            else:chunks.append(s)
        if name is not None:yield name,"".join(chunks).upper()
def save_fasta(path,records):
    with open(path,"w") as out:
        for name,seq in records:
            out.write(">"+name+"\n")
            for start in range(0,len(seq),80):out.write(seq[start:start+80]+"\n")
def canonical_mitos(raw):
    s=re.sub(r"^(gene_|transcript_|mt[_-])","",raw.lower())
    s=re.sub(r"_\d+$","",s).split("(")[0].replace("_","").replace("-","")
    if s in {"rrns","12s","rnr1"}:return "RNR1"
    if s in {"rrnl","16s","rnr2"}:return "RNR2"
    if s in {"cob","cytb","cyb"}:return "CYB"
    if re.fullmatch(r"nad(4l|[1-6])",s):return "ND"+s[3:].upper()
    if re.fullmatch(r"cox[123]",s):return "CO"+s[3:]
    if s in ("atp6","atp8"):return s.upper()
    return None
def canonical_trna(typ,anti):
    typ=typ.strip().lower();anti=anti.strip().upper().replace("U","T")
    if typ=="leu":return {"TAA":"TRNL1","TAG":"TRNL2"}.get(anti)
    if typ=="ser":return {"TGA":"TRNS1","GCT":"TRNS2"}.get(anti)
    aa=AA.get(typ)
    gene=("TRN"+aa) if aa else None
    return gene if gene in TRNA else None
@dataclass
class Feature:
    gene:str
    start:int # 0-based, may extend across circular origin
    end:int # exclusive
    strand:str
def read_gff(path,length):
    features=[]
    with open(path) as inp:
        for line in inp:
            if line.startswith("#"):continue
            parts=line.rstrip("\n").split("\t")
            if len(parts)<9 or parts[2].lower() not in ("gene","rrna"):continue
            try:a,b=int(parts[3])-1,int(parts[4])
            except ValueError:continue
            if a<0 or a>=length or b<=a or b-a>length:continue
            attrs=dict(pair.split("=",1) for pair in parts[8].split(";") if "=" in pair)
            gene=canonical_mitos(attrs.get("Name",attrs.get("gene_id",attrs.get("ID",""))))
            if gene and ((parts[2].lower()=="gene" and gene in CDS) or (parts[2].lower()=="rrna" and gene in RRNA)):
                features.append(Feature(gene,a,b,parts[6]))
    return features
def read_trnas(path):
    by=collections.defaultdict(list)
    with open(path,newline="") as inp:
        for row in csv.DictReader(inp):
            try:a,b=int(row["Begin"]),int(row["End"])
            except (KeyError,ValueError,TypeError):continue
            gene=canonical_trna(row["Type"],row["AntiCodon"])
            if gene:by[norm(row["Organism"])].append(Feature(gene,min(a,b)-1,max(a,b),"+" if a<b else "-"))
    return by
def pick_features(items,L):
    by=collections.defaultdict(list)
    for f in items:
        if 0<=f.start<L and f.start<f.end and f.end-f.start<=L:
            by[f.gene].append(f)
    selected={};ambiguous={}
    for gene,group in by.items():
        loci=[]
        for f in group:
            fs={i%L for i in range(f.start,f.end)}
            matching=[g for g in loci if g.strand==f.strand and len(fs & {i%L for i in range(g.start,g.end)})>=.8*min(len(fs),g.end-g.start)]
            if matching:
                if f.end-f.start>matching[0].end-matching[0].start:loci[loci.index(matching[0])]=f
            else:loci.append(f)
        if len(loci)==1:selected[gene]=loci[0]
        else:ambiguous[gene]=len(loci)
    return selected,ambiguous
def circular(seq,a,b):
    return "".join(seq[i%len(seq)] for i in range(a,b))

def align(records,folder,name,threads,mafft):
    if not records:return {}
    if len(records)==1:return dict(records)
    src=folder/(name+".fa");dest=folder/(name+".aln.fa")
    save_fasta(src,records)
    cmd=[mafft,"--retree","2","--maxiterate","0","--thread",str(threads),str(src)]
    with open(dest,"w") as fd:
        result=subprocess.run(cmd,stdout=fd,stderr=subprocess.PIPE,text=True)
    if result.returncode:raise RuntimeError("MAFFT failed for "+name+": "+result.stderr[-1500:])
    rows=dict(fasta(dest))
    if set(rows)!={x[0] for x in records} or len(set(map(len,rows.values())))!=1:
        raise ValueError("Invalid MAFFT output "+name)
    return rows
def write_tsv(path,rows,columns):
    with open(path,"w",newline="") as inp:
        w=csv.DictWriter(inp,delimiter="\t",fieldnames=columns)
        w.writeheader();w.writerows(rows)
def run(rank,fa,gffs,trnas,args):
    out=Path(args.output_dir)/rank;out.mkdir(parents=True,exist_ok=True)
    intermediate=out/"gene_alignments";intermediate.mkdir(exist_ok=True)
    ref=list(fasta(args.human_fasta))
    if len(ref)!=1 or len(ref[0][1])!=16569:raise ValueError("Need one 16569-bp human FASTA")
    human=ref[0][1];L=len(human)
    samples={};skipped=[];ambiguous_rows=[]
    for name,seq in fasta(fa):
        key=norm(name)
        path=gffs.get(key)
        if not path or not Path(path).exists():
            skipped.append((name,"missing_gff"));continue
        all_features=read_gff(path,len(seq))+trnas.get(key,[])
        features,amb=pick_features(all_features,len(seq))
        for gene,count in amb.items():ambiguous_rows.append((name,gene,count))
        if key=="homosapiens":
            if len(seq)!=16569:
                raise ValueError(f"Human rank FASTA is {len(seq)} bp, expected 16569. Replace the Homo sapiens record in this rank FASTA; no 16565-bp liftover is performed.")
            if seq!=human:
                raise ValueError("Human rank FASTA differs from --human-fasta. Use the same 16569-bp reference sequence in both.")
            if set(features)!=EXPECTED or amb:
                raise ValueError("Human gene set not unique: "+repr(sorted(EXPECTED-set(features)))+" "+repr(amb))
            # Annotations must be generated for this exact 16569-bp sequence.
            region_lengths=set()
            with open(path) as gff_in:
                for line in gff_in:
                    p=line.rstrip("\n").split("\t")
                    if len(p)>=5 and p[2].lower()=="region":
                        try:region_lengths.add(int(p[4])-int(p[3])+1)
                        except ValueError:pass
            if region_lengths and region_lengths!={16569}:
                raise ValueError(f"Human MITOS2 GFF region length {sorted(region_lengths)} != 16569; supply the 16569-bp human GFF.")
            for f in features.values():
                if f.end>16569:
                    raise ValueError(f"Human annotation exceeds 16569 bp: {f.gene} {f.start+1}-{f.end}")
        samples[key]=(name,seq,features)
    if "homosapiens" not in samples:raise ValueError("Homo sapiens not found")
    keys=["homosapiens"]+[key for key in samples if key!="homosapiens"]
    human_features=samples["homosapiens"][2]
    order=sorted(human_features,key=lambda gene:human_features[gene].start%L)
    # Full gene alignments are intentionally overlapping. Stitched reference bases
    # are emitted only once, and duplicate calls are audited.
    primary={};duplicates=[];blocks=[]
    def incorporate(aligned,positions,block,segment):
        reference=aligned["homosapiens"]
        pos_iter=iter(positions)
        for i,ch in enumerate(reference):
            if ch=="-":continue # Reference-insertion columns not part of human-coordinate MSA.
            try:pos=next(pos_iter)
            except StopIteration:raise ValueError("Reference alignment longer than block: "+block)
            value={key:aligned.get(key,"-"*len(reference))[i] for key in keys}
            if pos in primary:
                prev=primary[pos][0]
                disagree=sum(1 for key in keys if key in aligned and primary[pos][1].get(key,"-")!=value[key])
                duplicates.append((pos,prev,block,disagree))
            else:primary[pos]=(block,value,segment)
        if next(pos_iter,None) is not None:raise ValueError("Reference alignment shorter than block: "+block)
    for i,gene in enumerate(order):
        core=[]
        for key in keys:
            seq=samples[key][1];features=samples[key][2]
            f=features.get(gene)
            if f:
                s=circular(seq,f.start,f.end)
                core.append((key,rc(s) if f.strand=="-" else s))
        aligned=align(core,intermediate,f"{i:02d}_{gene}_core",args.threads,args.mafft)
        f=human_features[gene]
        if f.strand=="-":aligned={k:rc(v) for k,v in aligned.items()}
        coords=[j%L+1 for j in range(f.start,f.end)]
        if aligned["homosapiens"].replace("-","")!=circular(human,f.start,f.end):
            raise ValueError("Human reference changed in gene "+gene)
        incorporate(aligned,coords,gene,"gene")
        blocks.append((gene,"gene",len(core),"aligned"))
    # Intergenic regions: the union complement of gene features on human chrM.
    # Only accept other species if both bounding genes are physically adjacent.
    for i,gene in enumerate(order):
        prev=order[i-1];current_f=human_features[gene];previous_f=human_features[prev]
        prev_end=previous_f.end
        if i==0:prev_end-=L
        start=prev_end;stop=current_f.start
        if stop<=start:continue
        # Only positions still not covered by any gene are assigned to IGR.
        gap_coords=[j%L+1 for j in range(start,stop)]
        if not any(p not in primary for p in gap_coords):continue
        records=[]
        for key in keys:
            _,seq,features=samples[key];fs=sorted(features.values(),key=lambda f:f.start%len(seq))
            a=features.get(prev);b=features.get(gene)
            if a is None or b is None:continue
            ia=fs.index(a);ib=fs.index(b)
            if (ia+1)%len(fs)!=ib:continue
            local_end=a.end if ia!=len(fs)-1 else a.end-len(seq)
            if b.start<=local_end:continue
            if a.strand!=previous_f.strand or b.strand!=current_f.strand:continue
            records.append((key,circular(seq,local_end,b.start)))
        if not any(k=="homosapiens" for k,_ in records):raise ValueError("Missing human IGR")
        guarded=(len(records)>args.max_flank_sequences or max(len(s) for _,s in records)>args.max_flank_bp)
        if guarded:records=[(k,s) for k,s in records if k=="homosapiens"]
        aligned=align(records,intermediate,f"{i:02d}_{gene}_intergenic",args.threads,args.mafft)
        if aligned["homosapiens"].replace("-","")!=circular(human,start,stop):
            raise ValueError("Human intergenic mismatch: "+gene)
        incorporate(aligned,gap_coords,gene,"intergenic")
        blocks.append((gene,"intergenic",len(records),"guarded" if guarded else "aligned"))
    if set(primary)!=set(range(1,L+1)):
        raise ValueError(f"Human coverage {len(primary)}/{L}; missing positions {sorted(set(range(1,L+1))-set(primary))[:10]}")
    alignment={key:"".join(primary[i][1][key] for i in range(1,L+1)) for key in keys}
    if alignment["homosapiens"]!=human:raise ValueError("Human sequence mismatch")
    save_fasta(out/(rank+"_stitched_aln.fasta"),[(samples[key][0],alignment[key]) for key in keys])
    write_tsv(out/"human_coordinate_map.tsv",
              [(i,i,primary[i][0],primary[i][2]) for i in range(1,L+1)],
              ["alignment_column","human_position","primary_block","segment"])
    write_tsv(out/"overlap_dedup.tsv",duplicates,
              ["human_position","primary_block","duplicate_block","n_discordant_species"])
    write_tsv(out/"block_summary.tsv",blocks,["gene","segment","n_sequences","status"])
    write_tsv(out/"ambiguous_genes.tsv",ambiguous_rows,["species","gene","n_loci"])
    write_tsv(out/"skipped_species.tsv",skipped,["species","reason"])
    summary=dict(rank=rank,n_species=len(keys),reference_length=L,n_overlap_columns_deduplicated=len(duplicates))
    (out/"run_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    if not args.keep_intermediate:
        import shutil;shutil.rmtree(intermediate)
    LOG.info("%s completed: %s",rank,summary)
def main():
    p=argparse.ArgumentParser()
    for x in ("fasta-dir","annotation-summary","trnascan-csv","human-fasta","output-dir"):
        p.add_argument("--"+x,required=True)
    p.add_argument("--ranks",nargs="+",default=None)
    p.add_argument("--threads",type=int,default=4)
    p.add_argument("--mafft",default="mafft")
    p.add_argument("--max-flank-sequences",type=int,default=300)
    p.add_argument("--max-flank-bp",type=int,default=2500)
    p.add_argument("--keep-intermediate",action="store_true")
    args=p.parse_args()
    logging.basicConfig(level=logging.INFO,format="%(asctime)s %(levelname)s %(message)s")
    with open(args.annotation_summary) as inp:
        gffs={norm(r["species"]):r["result_gff"] for r in csv.DictReader(inp,delimiter="\t") if r.get("result_gff")}
    trnas=read_trnas(args.trnascan_csv)
    paths=list(Path(args.fasta_dir).glob("*_chrm.faa"))
    if not paths:raise ValueError("No *_chrm.faa FASTA files found")
    for path in paths:
        rank=path.name.removesuffix("_chrm.faa")
        if args.ranks is None or rank in args.ranks:run(rank,path,gffs,trnas,args)
if __name__=="__main__":main()
