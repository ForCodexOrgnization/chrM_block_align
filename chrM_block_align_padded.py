#!/usr/bin/env python3
"""Experimental padded gene-block MSA; separate from the running baseline."""
import argparse,csv,json,logging
from pathlib import Path
from collections import defaultdict
from chrM_block_align import Feature,EXPECTED,fasta,save_fasta,norm,rc,circular,read_gff,read_trnas,pick_features,align,write_tsv
from human_genbank import load_human_genbank

def ordered(fs,L):return sorted(fs,key=lambda g:fs[g].start%L)
def neighbours(fs,g,reference,order,L):
    if g not in fs:return (False,False)
    a=ordered(fs,L);i=a.index(g);j=order.index(g)
    left,right=order[j-1],order[(j+1)%len(order)]
    return (a[i-1]==left and left in fs and fs[left].strand==reference[left].strand,
            a[(i+1)%len(a)]==right and right in fs and fs[right].strand==reference[right].strand)

def run(rank,fa,gffs,trnas,args):
    out=Path(args.output_dir)/rank;out.mkdir(parents=True,exist_ok=True)
    msadir=out/"gene_alignments";msadir.mkdir(exist_ok=True)
    ref=list(fasta(args.human_fasta))
    if len(ref)!=1 or len(ref[0][1])!=16569:raise ValueError("Need one canonical 16569-bp human FASTA")
    human=ref[0][1];L=len(human)
    genes=load_human_genbank(args.human_genbank,human,Feature,EXPECTED)
    order=ordered(genes,L)
    samples={};skipped=[]
    for name,seq in fasta(fa):
        k=norm(name)
        if k in samples:raise ValueError("Duplicate species: "+k)
        if k=="homosapiens":
            if seq!=human:raise ValueError("Human reference mismatch")
            fs=genes
        else:
            gf=gffs.get(k)
            if not gf or not Path(gf).exists():
                skipped.append((name,"missing_gff"));continue
            fs,_=pick_features(read_gff(gf,len(seq))+trnas.get(k,[]),len(seq))
        samples[k]=(name,seq,fs)
    if "homosapiens" not in samples:raise ValueError("Missing human")
    keys=["homosapiens"]+[k for k in samples if k!="homosapiens"]
    proposed=defaultdict(list);flankqc=[];blocks=[]
    for idx,g in enumerate(order):
        h=genes[g];beg=h.start-args.padding;end=h.end+args.padding
        coords=[i%L+1 for i in range(beg,end)]
        records=[];allowed={}
        for k in keys:
            _,seq,fs=samples[k]
            if g not in fs:continue
            f=fs[g];left,right=neighbours(fs,g,genes,order,len(seq))
            if k=="homosapiens":left=right=True
            raw=circular(seq,f.start-(args.padding if left else 0),f.end+(args.padding if right else 0))
            records.append((k,rc(raw) if f.strand=="-" else raw))
            allowed[k]=(left,right)
            flankqc.append((g,k,int(left),int(right)))
        aligned=align(records,msadir,f"{idx:02d}_{g}_pad{args.padding}",args.threads,args.mafft)
        if h.strand=="-":aligned={k:rc(v) for k,v in aligned.items()}
        if aligned["homosapiens"].replace("-","")!=circular(human,beg,end):
            raise ValueError("Human block mismatch "+g)
        iterator=iter(coords)
        for col,nt in enumerate(aligned["homosapiens"]):
            if nt=="-":continue
            pos=next(iterator)
            left_dist=(h.start-(pos-1))%L
            right_dist=((pos-1)-(h.end-1))%L
            core=(pos-1-h.start)%L < h.end-h.start
            is_left=not core and left_dist<=args.padding
            dist=0 if core else min(left_dist,right_dist)
            calls={}
            for k in keys:
                v=aligned[k][col] if k in aligned else "-"
                if k!="homosapiens" and not core and k in allowed:
                    left,right=allowed[k]
                    if (is_left and not left) or (not is_left and not right):v="-"
                calls[k]=v
            proposed[pos].append((0 if core else 1,dist,idx,g,"gene" if core else "padding",calls))
        blocks.append((g,len(records),len(coords)))
    igrqc=[]
    for idx,g in enumerate(order):
        prev=order[idx-1];a=genes[prev];b=genes[g]
        start=a.end-(L if idx==0 else 0);end=b.start
        if end<=start:continue
        coords=[i%L+1 for i in range(start,end)]
        control=(prev,g)==("TRNP","TRNF")
        if not control and all(p in proposed for p in coords):continue
        records=[]
        for k in keys:
            _,seq,fs=samples[k]
            if prev not in fs or g not in fs:continue
            phys=ordered(fs,len(seq));i=phys.index(prev)
            if phys[(i+1)%len(phys)]!=g:continue
            x,y=fs[prev],fs[g]
            if x.strand!=a.strand or y.strand!=b.strand:continue
            gap=(y.start-x.end)%len(seq)
            raw=circular(human,start,end) if k=="homosapiens" else circular(seq,x.end,x.end+gap)
            if k!="homosapiens" and len(raw)>args.max_flank_bp:continue
            records.append((k,raw))
        if not records or records[0][0]!="homosapiens":raise ValueError("No human IGR")
        aligned=align(records,msadir,f"{idx:02d}_{g}_igr",args.threads,args.mafft,
                      incremental=True,batch_size=args.flank_batch_size)
        if aligned["homosapiens"].replace("-","")!=circular(human,start,end):
            raise ValueError("Human IGR mismatch")
        iterator=iter(coords)
        for col,ch in enumerate(aligned["homosapiens"]):
            if ch=="-":continue
            p=next(iterator)
            if p not in proposed or control:
                calls={k:aligned[k][col] if k in aligned else "-" for k in keys}
                proposed[p].append((2,0,idx,g,"control_igr" if control else "intergenic",calls))
        igrqc.append((prev,g,len(records),len(coords),int(control)))
    if len(proposed)!=L:raise ValueError(f"Human coordinate coverage: {len(proposed)}/{L}")
    owners={p:min(proposed[p],key=lambda v:v[:3]) for p in range(1,L+1)}
    seqs={k:"".join(owners[p][5][k] for p in range(1,L+1)) for k in keys}
    if seqs["homosapiens"]!=human:raise ValueError("Stitched human mismatch")
    save_fasta(out/f"{rank}_stitched_aln.fasta",[(samples[k][0],seqs[k]) for k in keys])
    write_tsv(out/"human_coordinate_map.tsv",[(p,p,owners[p][3],owners[p][4]) for p in range(1,L+1)],
              ["alignment_column","human_position","primary_block","segment"])
    write_tsv(out/"stitch_decisions.tsv",[(p,owners[p][3],owners[p][4],len(proposed[p]),owners[p][1]) for p in range(1,L+1)],
              ["human_position","primary_block","segment","n_candidates","anchor_distance"])
    conflicts=[]
    for p in range(1,L+1):
        winner=owners[p]
        for alternative in proposed[p]:
            if alternative is winner:continue
            n=sum(winner[5][k]!=alternative[5][k] and winner[5][k] in "ACGT" and alternative[5][k] in "ACGT" for k in keys)
            conflicts.append((p,winner[3],alternative[3],n))
    write_tsv(out/"overlap_conflicts.tsv",conflicts,
              ["human_position","primary_block","alternative_block","n_discordant_callable_species"])
    write_tsv(out/"flank_context_qc.tsv",flankqc,["gene","species","left_supported","right_supported"])
    write_tsv(out/"padded_block_summary.tsv",blocks,["gene","n_sequences","human_block_bp"])
    write_tsv(out/"intergenic_qc.tsv",igrqc,["left_gene","right_gene","n_sequences","human_gap_bp","control_region"])
    write_tsv(out/"skipped_species.tsv",skipped,["species","reason"])
    (out/"run_summary.json").write_text(json.dumps({"rank":rank,"padding":args.padding,
         "n_species":len(keys),"reference_length":L,"n_overlap_candidates":len(conflicts)},indent=2)+"\n")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--fasta-dir",default="data/unrot_fasta")
    p.add_argument("--annotation-summary",default="data/mitos2_species_summary.tsv")
    p.add_argument("--trnascan-csv",default="data/consolidated_trnascan-2.csv")
    p.add_argument("--human-fasta",default="data/Humans_mtdna.fasta")
    p.add_argument("--human-genbank",default="data/sequence.gb")
    p.add_argument("--output-dir",default="block_results_padded_v1")
    p.add_argument("--ranks",nargs="+",default=["class"])
    p.add_argument("--padding",type=int,choices=[50,100],default=50)
    p.add_argument("--threads",type=int,default=4)
    p.add_argument("--mafft",default="mafft")
    p.add_argument("--max-flank-bp",type=int,default=2500)
    p.add_argument("--flank-batch-size",type=int,default=100)
    a=p.parse_args()
    if Path(a.output_dir).resolve()==Path("block_results_v6").resolve():p.error("Cannot overwrite baseline output")
    with open(a.annotation_summary) as f:
        gffs={norm(row["species"]):row["result_gff"] for row in csv.DictReader(f,delimiter="\t") if row.get("result_gff")}
    trnas=read_trnas(a.trnascan_csv)
    for rank in a.ranks:run(rank,Path(a.fasta_dir)/f"{rank}_chrm.faa",gffs,trnas,a)
if __name__=="__main__":main()
