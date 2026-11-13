#!/usr/bin/env python3
"""Regression tests for circular overlap and coordinate conversion."""
import unittest
from compare_pairwise_cr_to_longest_intervals import (
    in_gap, original_pos, cr_mapping_group, overlap_label, ori_group,
    ori_overlap_value
)

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


class TestIndependentOverlapClassification(unittest.TestCase):
    def test_below_threshold_remains_valid_mapping(self):
        self.assertEqual(cr_mapping_group("validated", 0.30, 0.50),
                         "below_threshold")
        self.assertEqual(cr_mapping_group("validated", 0.50, 0.50),
                         "at_or_above_threshold")
        self.assertEqual(cr_mapping_group("invalid_alignment", 0.0, 0.50),
                         "no_valid_mapping")

    def test_overlap_is_tri_state(self):
        self.assertEqual(overlap_label(1), "yes")
        self.assertEqual(overlap_label(0), "no")
        self.assertEqual(overlap_label(""), "unavailable")

    def test_oriH_is_independent_of_CR_mapping_fraction(self):
        row = dict(mitos2_gff_available="1", n_OriH_features="1",
                   OriH_overlaps="1", OriH_overlap_bp="140")
        self.assertEqual(ori_group(row, True), "OriH_overlaps_longest")
        self.assertEqual(ori_overlap_value(row, True), 1)
        row["OriH_overlaps"] = "0"
        self.assertEqual(ori_overlap_value(row, True), 0)

    def test_missing_OriH_prediction_not_equivalent_to_non_overlap(self):
        row = dict(mitos2_gff_available="1", n_OriH_features="0",
                   OriH_overlaps="0")
        self.assertEqual(ori_group(row, True), "no_OriH_prediction")
        self.assertEqual(ori_overlap_value(row, True), "")

    def test_no_longest_interval_is_unavailable(self):
        row = dict(mitos2_gff_available="1", n_OriH_features="1",
                   OriH_overlaps="0")
        self.assertEqual(ori_overlap_value(row, False), "")

if __name__ == "__main__":
    unittest.main()
