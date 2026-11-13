#!/usr/bin/env python3
"""Regression tests for circular overlap and coordinate conversion."""
import unittest
from compare_pairwise_cr_to_longest_intervals import in_gap, original_pos

class TestCircularIntervals(unittest.TestCase):
    def test_end_zero_means_last_base(self):
        # Legacy longest-gap TSV writes 0 for a gap ending at the last bp.
        self.assertTrue(in_gap(16000, 15001, 0, 16000))
        self.assertTrue(in_gap(15001, 15001, 0, 16000))
        self.assertFalse(in_gap(1, 15001, 0, 16000))

    def test_wraparound_interval(self):
        self.assertTrue(in_gap(16000, 15900, 400, 16000))
        self.assertTrue(in_gap(100, 15900, 400, 16000))
        self.assertFalse(in_gap(500, 15900, 400, 16000))

    def test_normal_interval(self):
        self.assertTrue(in_gap(100, 100, 120, 16000))
        self.assertFalse(in_gap(121, 100, 120, 16000))

    def test_one_based_rotation(self):
        # When the rotated reference begins at original position 4001.
        self.assertEqual(original_pos(1, 4001, 16000), 4001)
        self.assertEqual(original_pos(12000, 4001, 16000), 16000)
        self.assertEqual(original_pos(12001, 4001, 16000), 1)

if __name__ == "__main__":
    unittest.main()
