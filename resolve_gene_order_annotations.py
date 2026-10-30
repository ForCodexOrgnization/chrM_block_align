#!/usr/bin/env python3
"""Conservative gene-order-assisted disambiguation for GenBank feature labels.

Reference order comes from canonical human GenBank Feature coordinates.
Only a UNIQUE exact circular gene-order solution may resolve an ambiguous name.
Gene order can be reversed for diagnostic purposes, but no FASTA or annotation
coordinates are rotated or reverse-complemented by this module.
"""
from collections import Counter
from itertools import product
import json
import re

FIELDS=("Gene","G_Symbol","G_locus_tag","G_Sequence","G_tag",
        "alt_name","alt_synonym","alt_Product")

def key(s):
    return str(s or "").strip().casefold()

def normalized(s):
    return re.sub(r"[^a-z0-9]","",key(s))

def candidates(row,aliases):
    hits=set()
    for field in FIELDS:
        value=row.get(field,"")
        if key(value) in ("","???","na","nan","none","n/a","-"):continue
        hits.update(aliases.get(normalized(value),()))
        k=normalized(value)
        if k in ("trnaleu","trnl","mttl"):
            hits.update(("TRNL1","TRNL2"))
        elif k in ("trnaser","trns","mtts"):
            hits.update(("TRNS1","TRNS2"))
    return hits

def signature(seq):
    # Canonical orientation-invariant circular representation for comparison.
    seq=tuple(seq)
    if not seq:return ()
    return min(seq[i:]+seq[:i] for i in range(len(seq)))

def equivalent(reference,observed,allow_reverse=True):
    if len(reference)!=len(observed):return False
    variants=[observed]
    if allow_reverse:
        variants.append([(gene,"+" if strand=="-" else "-") for gene,strand in reversed(observed)])
    return any(signature(reference)==signature(seq) for seq in variants)

def contexts(seq):
    result={}
    if len(seq)<3:return result
    for i,(name,_) in enumerate(seq):
        result.setdefault(name,Counter())[(seq[i-1][0],seq[(i+1)%len(seq)][0])]+=1
    return result

def score(reference,observed):
    # For QC only: remove unequal-copy gene groups before comparing context.
    refcounts=Counter(g for g,_ in reference)
    obcounts=Counter(g for g,_ in observed)
    keep={g for g in refcounts if refcounts[g]==obcounts[g]}
    r=[(g,s) for g,s in reference if g in keep]
    q=[(g,s) for g,s in observed if g in keep]
    cr,cq=contexts(r),contexts(q)
    bad=sum(cr.get(g)!=cq.get(g) for g in keep)
    strand_ref={g:Counter(s for n,s in r if n==g) for g in keep}
    strand_obs={g:Counter(s for n,s in q if n==g) for g in keep}
    bad_strands=sum(strand_ref[g]!=strand_obs[g] for g in keep)
    copy=sum(abs(refcounts[g]-obcounts[g]) for g in set(refcounts)|set(obcounts))
    return (copy,bad,bad_strands)

def resolve(rows,reference,max_candidates=4096):
    """rows: sorted (candidate_set, strand, metadata) for one organism.
    Return (resolved gene or None per source row, detailed QC).
    """
    options=[tuple(sorted(c)) for c,_,_ in rows]
    variable=[i for i,c in enumerate(options) if len(c)>1]
    choices=1
    for i in variable:choices*=len(options[i])
    original=[next(iter(c)) if len(c)==1 else None for c in options]
    detail=[]
    if not variable:return original,detail
    if choices>max_candidates:
        return original,[{"index":i,"candidates":json.dumps(options[i]),
               "status":"search_limit_exceeded","assigned":""} for i in variable]
    best=None;best_assignments=[]
    # One row maps to one gene identity. Do not force equality of copy counts
    # without checking the entire annotated sequence.
    for selection in product(*(options[i] for i in variable)):
        labels=original.copy()
        for i,name in zip(variable,selection):labels[i]=name
        observed=[(g,rows[i][1]) for i,g in enumerate(labels) if g is not None]
        is_exact=equivalent(reference,observed)
        cost=(0 if is_exact else 1,*score(reference,observed))
        if best is None or cost<best:
            best=cost;best_assignments=[labels]
        elif cost==best:best_assignments.append(labels)
    # Do not automatically accept inexact best-fit assignments. Treat their
    # score as diagnostic because distant taxa can truly rearrange genes.
    for i in variable:
        supported=sorted({trial[i] for trial in best_assignments})
        if best[0]==0 and len(supported)==1:
            original[i]=supported[0];status="unique_exact_order"
        elif best[0]==0:status="ambiguous_exact_order"
        elif len(supported)==1:status="inferred_best_fit_not_accepted"
        else:status="ambiguous_best_fit"
        detail.append({"index":i,"candidates":json.dumps(options[i]),
            "supported":json.dumps(supported),"assigned":original[i] or "",
            "status":status,"best_score":json.dumps(best),
            "n_best_assignments":len(best_assignments)})
    return original,detail
