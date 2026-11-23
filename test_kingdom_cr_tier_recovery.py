#!/usr/bin/env python3
"""Pure synthetic checks for Kingdom incremental CR recovery."""
from audit_genbank_control_region_vocabulary import classify
from analyze_kingdom_cr_tier_recovery import candidate_tier
from compare_genbank_cr_to_longest_orih import evaluate_regions


def region(source,text,coords):
    return dict(source=source,note=text,intervals=coords,species="test",
                acc="NC_1",cr_label="CR_unspecified")


def main():
    assert classify("conrol region")[0]=="explicit_cr"
    assert classify("control region1")[0]=="explicit_cr"
    assert classify("control region2")[0]=="explicit_cr"
    assert classify("putative control region","D-loop")[0]=="putative_cr"
    assert classify("(A+T)-rich region; putative")[0]=="at_rich_alternative"
    assert classify("non-coding region")[1]=="tier2_review"
    assert classify("control region CSB-1")[0]=="subfeature"
    assert candidate_tier(region("misc_feature","contral region",[(4,25)]))=="tier1_explicit"
    assert candidate_tier(region("misc_feature","AT-loop",[(4,25)]))=="tier2_at_rich"
    assert candidate_tier(region("misc_feature","CSB-2",[(4,25)])) is None
    a=region("D-loop","control region 1",[(1,5)])
    b=region("D-loop","control region 2",[(6,10)])
    assert evaluate_regions([a,b],set(range(1,11)),100)[0] is False
    print("PASS Kingdom tier recovery vocabulary and independent CR copies")


if __name__=="__main__":
    main()
