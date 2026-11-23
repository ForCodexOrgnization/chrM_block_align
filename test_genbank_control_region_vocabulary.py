#!/usr/bin/env python3
"""Unit tests for Kingdom control-region synonym classification."""
from audit_genbank_control_region_vocabulary import classify

def verify(text,expected):
    result=classify(text)[0]
    assert result==expected,(text,result,expected)

def main():
    for s in ("control region","Control Region","contral region",
              "control region 1","control region II","D-loop",
              "A+T-rich control region"):
        verify(s,"explicit_cr")
    verify("putative control region 2","putative_cr")
    verify("control region CSB-1","subfeature")
    verify("control region coreTAS sequence","subfeature")
    verify("control region SSB-3","subfeature")
    verify("remnant control region","remnant")
    verify("A+T rich region","at_rich_alternative")
    verify("AT-loop","at_rich_alternative")
    verify("A-T rich region","at_rich_alternative")
    verify("major non-coding region","long_noncoding_alternative")
    verify("origin of replication","replication_origin_only")
    verify("chloroplast-derived region","unrelated")
    print("PASS GenBank control region vocabulary tests")

if __name__=="__main__":
    main()
