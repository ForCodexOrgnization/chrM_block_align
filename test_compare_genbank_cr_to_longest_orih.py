#!/usr/bin/env python3
"""Regression tests for strict single-interval GenBank CR containment."""
from compare_genbank_cr_to_longest_orih import (
    coordinates_to_set, evaluate_longest, evaluate_regions, deduplicate)


def region(a,b,label="CR_unspecified",source="D-loop"):
    return dict(species="test",acc="NC_TEST",intervals=[(a,b)],
                classification="full_control_region",cr_label=label,
                source=source)


def test_distinct_regions_do_not_combine():
    ori=set(range(10,21))
    cr1=region(10,15,"CR1")
    cr2=region(16,20,"CR2")
    success,matching,best=evaluate_regions([cr1,cr2],ori,100)
    assert not success and matching==[] and best==6
    both=region(8,22,"CR2")
    success,matching,best=evaluate_regions([cr1,both],ori,100)
    assert success and matching==[both] and best==11


def test_exact_dedup_and_separate_cr():
    first=region(10,20,"CR1")
    duplicate=region(10,20,"CR_unspecified",source="misc_feature")
    second=region(21,40,"CR2")
    grouped=deduplicate([first,duplicate,second])
    assert len(grouped["test"])==2
    assert grouped["test"][0]["sources"]=={"D-loop","misc_feature"}
    assert {r["cr_label"] for r in grouped["test"]}=={"CR1","CR2"}


def test_circular_containment():
    ori={98,99,100,1,2}
    cr=dict(region(98,2),intervals=[(98,100),(1,2)])
    positions=coordinates_to_set(cr["intervals"],100)
    assert positions==ori
    success,matching,_=evaluate_regions([cr],ori,100)
    assert success and len(matching)==1
    gap=dict(candidate_defined="1",gap_start_1based="98",gap_end_1based="2")
    success,bp=evaluate_longest(gap,ori,100)
    assert success and bp==5
    assert not evaluate_longest(gap,ori|{3},100)[0]


def test_subregion_not_considered_full():
    r=region(10,20)
    r["classification"]="subfeature"
    assert deduplicate([r])=={}


def test_expanded_genbank_vocabulary():
    from compare_genbank_cr_to_longest_orih import deduplicate
    def entry(note, loc, source="misc_feature"):
        return dict(species="speciesx", acc="NC_TEST", source=source,
                    note=note, intervals=[loc], cr_label="CR_unspecified",
                    classification="expanded_candidate")
    regions=deduplicate([
        entry("contral region", (1, 100)),
        entry("AT-loop", (200, 300)),
        entry("non-coding region", (400, 500)),
        entry("control region CSB-1", (600, 620)),
        entry("remnant control region", (700, 720)),
        entry("origin of replication", (800, 850)),
    ])
    assert len(regions["speciesx"]) == 3
    assert {r["vocabulary_tier"] for r in regions["speciesx"]} == {
        "tier1", "tier2", "tier2_review"
    }


if __name__=="__main__":
    test_distinct_regions_do_not_combine()
    test_exact_dedup_and_separate_cr()
    test_circular_containment()
    test_subregion_not_considered_full()
    test_expanded_genbank_vocabulary()
    print("PASS strict 100% OriH containment tests")
