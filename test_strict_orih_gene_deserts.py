#!/usr/bin/env python3
"""Small synthetic checks for strict OriH containment and circular gaps.

Run: python -m unittest test_strict_orih_gene_deserts.py
"""
import unittest
from analyze_strict_orih_gene_deserts import circular_gaps, select_interval
from analyze_longest_all_gene_desert import spans

class TestStrictOriginIntervals(unittest.TestCase):
    def test_circular_wrap_gap(self):
        # Genes occupy 101..300 and 501..900 (1-based inclusive).
        # Longest gap wraps 901..1000 + 1..100 = 200 bp.
        occupied=spans(101,300,1000)+spans(501,900,1000)
        gaps=circular_gaps(occupied,1000)
        self.assertIn((900,1100),gaps)
        origins=[spans(950,50,1000)]
        _,winner,count,length,_=select_interval(occupied,origins,1000)
        self.assertEqual(winner,(900,1100))
        self.assertEqual(count,1)
        self.assertEqual(length,101)

    def test_full_containment_not_one_bp_overlap(self):
        occupied=spans(101,200,1000)+spans(401,900,1000)
        # Gap 201-400 overlaps OH 350-450 but does not contain it.
        origins=[spans(350,450,1000)]
        _,winner,count,_,_=select_interval(occupied,origins,1000)
        self.assertIsNone(winner)
        self.assertEqual(count,0)

    def test_choose_longest_only_after_full_containment(self):
        occupied=spans(101,200,1000)+spans(301,700,1000)
        # The long wraparound gap (~400bp) does not contain OH.
        # The short internal gap (100bp) fully contains OH.
        origins=[spans(220,240,1000)]
        gaps,winner,count,_,_=select_interval(occupied,origins,1000)
        self.assertEqual(len(gaps),2)
        self.assertEqual(winner,(200,300))
        self.assertEqual(count,1)

    def test_multiple_origins_must_fit_in_same_interval(self):
        occupied=spans(201,300,1000)+spans(601,700,1000)
        origins=[spans(400,450,1000),spans(800,850,1000)]
        _,winner,count,_,_=select_interval(occupied,origins,1000)
        self.assertIsNone(winner)
        self.assertEqual(count,0)

if __name__=="__main__":
    unittest.main()
