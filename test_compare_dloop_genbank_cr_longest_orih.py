#!/usr/bin/env python3
"""Synthetic regression tests for independent D-loop and GenBank-CR scoring."""
from compare_dloop_genbank_cr_longest_orih import score_methods
from compare_genbank_cr_to_longest_orih import deduplicate


def region(source,note,a,b):
    return dict(species="test",acc="NC_TEST",source=source,note=note,
                intervals=[(a,b)],cr_label="CR_unspecified",
                classification="full_control_region")


def main():
    ori=set(range(10,21))
    gap={"candidate_defined":"1","gap_start_1based":"30","gap_end_1based":"60"}
    dloop=deduplicate([region("D-loop","D-loop",10,14)])
    all_cr=deduplicate([region("D-loop","D-loop",10,14),
                        region("misc_feature","AT-loop",10,20)])
    flags,_=score_methods(gap,ori,100,dloop["test"],all_cr["test"])
    assert flags==dict(Longest=0,GenBank_Dloop=0,GenBank_CR=1,Either=1),flags

    # CR1 and CR2 must not be combined into one 100%-containing region.
    dloop=deduplicate([region("D-loop","D-loop",10,15),
                       region("D-loop","D-loop",16,20)])
    flags,_=score_methods(gap,ori,100,dloop["test"],dloop["test"])
    assert flags==dict(Longest=0,GenBank_Dloop=0,GenBank_CR=0,Either=0),flags

    # Longest-only success must show up in Either.
    gap={"candidate_defined":"1","gap_start_1based":"10","gap_end_1based":"20"}
    flags,_=score_methods(gap,ori,100,[],[])
    assert flags==dict(Longest=1,GenBank_Dloop=0,GenBank_CR=0,Either=1),flags

    # Without a predicted OriH, no method can be counted as successful.
    flags,_=score_methods(gap,None,100,[],[])
    assert sum(flags.values())==0
    print("PASS Longest / D-loop-only / GenBank CR / Either scoring")


if __name__=="__main__":
    main()
