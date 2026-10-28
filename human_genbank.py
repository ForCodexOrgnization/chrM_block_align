"""Read the 37 canonical human mitochondrial genes directly from GenBank."""
from collections import Counter

def load_human_genbank(path, reference_sequence, Feature, expected):
    try:
        from Bio import SeqIO
        from Bio.SeqFeature import CompoundLocation
    except ImportError as exc:
        raise RuntimeError("Biopython is required for --human-genbank. Install: python3 -m pip install biopython") from exc

    record=SeqIO.read(path,"genbank")
    seq=str(record.seq).upper()
    if len(seq)!=16569:
        raise ValueError(f"Human GenBank sequence is {len(seq)} bp; expected 16569")
    if seq!=reference_sequence.upper():
        mismatches=sum(a!=b for a,b in zip(seq,reference_sequence.upper()))
        raise ValueError(f"Human GenBank sequence differs from --human-fasta at {mismatches} positions")
    features={}
    counts=Counter()
    types={"CDS","rRNA","tRNA"}
    for item in record.features:
        if item.type not in types:
            continue
        names=item.qualifiers.get("gene",[])
        if not names:
            continue
        gene=names[0].upper().replace("MT-","")
        if gene not in expected:
            continue
        if isinstance(item.location,CompoundLocation):
            raise ValueError(f"Compound human gene location requires review: {gene}")
        start,end=int(item.location.start),int(item.location.end)
        if not (0<=start<end<=len(seq)):
            raise ValueError(f"Out-of-range GenBank location: {gene}: {start}-{end}")
        strand={1:"+",-1:"-"}.get(item.location.strand)
        if strand is None:
            raise ValueError(f"Unknown human gene strand: {gene}")
        counts[gene]+=1
        if counts[gene]!=1:
            raise ValueError(f"Duplicate human GenBank annotation: {gene}")
        features[gene]=Feature(gene,start,end,strand)
    if set(features)!=expected:
        raise ValueError("Human GenBank gene set mismatch. Missing: "
                         +repr(sorted(expected-set(features)))
                         +"; unexpected: "+repr(sorted(set(features)-expected)))
    print(f"Human GenBank verified: accession={record.id}, bp={len(seq)}, genes={len(features)}")
    return features
