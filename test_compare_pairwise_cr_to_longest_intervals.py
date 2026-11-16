#!/usr/bin/env python3
"""Regression tests for circular overlap and coordinate conversion."""
import unittest
import tempfile
from pathlib import Path
from compare_pairwise_cr_to_longest_intervals import (
    in_gap, original_pos, cr_mapping_group, overlap_label, ori_group,
    ori_overlap_value, read_orih_positions, orih_guided_strategy,
    orih_coverage_flags, longest_orih_coverage, strategy_summary,
    mapped_human_cr, gap_column_in_human_cr
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


class TestDirectCRvsOriH(unittest.TestCase):
    def test_oriH_only_parent_feature(self):
        gff = (
            "mtDNA\tMITOS2\torigin_of_replication\t95\t5\t.\t+\t.\tID=OH;Name=OH\n"
            "mtDNA\tMITOS2\texon\t96\t98\t.\t+\t.\tID=OH-b;Parent=OH\n"
            "mtDNA\tMITOS2\torigin_of_replication\t20\t25\t.\t+\t.\tID=OL;Name=OL\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            filename = Path(tmp) / "test.gff"
            filename.write_text(gff)
            intervals, bases = read_orih_positions(filename, 100)
        self.assertEqual(intervals, [(95, 5)])
        self.assertEqual(bases, set(range(95, 101)) | set(range(1, 6)))
        self.assertEqual(len({1, 2, 40} & bases), 2)

    def test_two_features_and_union(self):
        gff = (
            "chrM\tx\torigin_of_replication\t10\t20\t.\t+\t.\tName=OH\n"
            "chrM\tx\trep_origin\t18\t30\t.\t+\t.\tName=OriH\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            filename = Path(tmp) / "test.gff"
            filename.write_text(gff)
            intervals, bases = read_orih_positions(filename, 100)
        self.assertEqual(len(intervals), 2)
        self.assertEqual(len(bases), 21)

    def test_fallback_order_and_intersection(self):
        self.assertEqual(orih_guided_strategy(1, 1), {
            "longest_first": "longest_interval",
            "CR_first": "pairwise_CR", "selected": 1, "both": 1})
        self.assertEqual(orih_guided_strategy(1, 0)["CR_first"],
                         "longest_interval")
        self.assertEqual(orih_guided_strategy(0, 1)["longest_first"],
                         "pairwise_CR")
        self.assertEqual(orih_guided_strategy("", 1)["selected"], 1)
        self.assertEqual(orih_guided_strategy(0, "")["selected"], 0)
        self.assertEqual(orih_guided_strategy(0, "")["CR_first"], "none")



class TestOriHCoverageThresholds(unittest.TestCase):
    def test_exact_thresholds(self):
        self.assertEqual(orih_coverage_flags(1, 100)["any"], 1)
        self.assertEqual(orih_coverage_flags(79, 100)["ge80"], 0)
        self.assertEqual(orih_coverage_flags(80, 100)["ge80"], 1)
        self.assertEqual(orih_coverage_flags(80, 100)["full"], 0)
        self.assertEqual(orih_coverage_flags(100, 100)["full"], 1)
        self.assertEqual(orih_coverage_flags(0, 100)["any"], 0)

    def test_missing_and_zero_oriH_unavailable(self):
        for x, n in (("", 20), (0, ""), (0, 0)):
            flags = orih_coverage_flags(x, n)
            self.assertEqual(flags["any"], "")
            self.assertEqual(flags["ge80"], "")
            self.assertEqual(flags["full"], "")

    def test_circular_longest_interval(self):
        # An interval running 91..100 contains both OriH pieces 95..100 and 1..5
        ori = set(range(95, 101)) | set(range(1, 6))
        g = dict(candidate_defined="1", gap_start_1based="91",
                 gap_end_1based="5")
        self.assertEqual(longest_orih_coverage(g, 100, ori), 11)
        g["gap_start_1based"] = "97"
        self.assertEqual(longest_orih_coverage(g, 100, ori), 9)
        g["gap_end_1based"] = "0"
        self.assertEqual(longest_orih_coverage(g, 100, ori), 4)

    def test_threshold_fallback_and_both(self):
        def row(l, c):
            return {
                "Longest_OriH_any_exact": int(l > 0),
                "Longest_OriH_ge80": int(l >= 80),
                "Longest_OriH_full": int(l == 100),
                "CR_OriH_any_exact": int(c > 0),
                "CR_OriH_ge80": int(c >= 80),
                "CR_OriH_full": int(c == 100),
            }
        rows = [row(100, 100), row(90, 0), row(0, 80), row(5, 5)]
        full = strategy_summary(rows, "full")
        self.assertEqual(full["longest_first_total"], 1)
        self.assertEqual(full["both_methods"], 1)
        eighty = strategy_summary(rows, "ge80")
        self.assertEqual(eighty["longest_first_primary"], 2)
        self.assertEqual(eighty["longest_first_CR_fallback"], 1)
        self.assertEqual(eighty["both_methods"], 1)
        any_result = strategy_summary(rows, "any")
        self.assertEqual(any_result["longest_first_total"], 4)
        self.assertEqual(any_result["both_methods"], 2)


class TestHumanCRInsertionMapping(unittest.TestCase):
    def test_gap_assignment_boundaries(self):
        self.assertTrue(gap_column_in_human_cr(0))
        self.assertTrue(gap_column_in_human_cr(1))
        self.assertTrue(gap_column_in_human_cr(575))
        self.assertFalse(gap_column_in_human_cr(576))
        self.assertFalse(gap_column_in_human_cr(16023))
        self.assertTrue(gap_column_in_human_cr(16024))
        self.assertTrue(gap_column_in_human_cr(16569))

    def test_internal_species_insertions_retained(self):
        # Human 1,2,'-','-',3,4; two species bases are insertions.
        expanded, anchor, n, h, sp = mapped_human_cr(
            "AC--GT", "ACTGGT", 1, 6)
        self.assertEqual(anchor, {1, 2, 5, 6})
        self.assertEqual(expanded, {1, 2, 3, 4, 5, 6})
        self.assertEqual((n, h, sp), (4, 4, 6))

    def test_species_deletion_does_not_create_bases(self):
        expanded, anchor, n, h, sp = mapped_human_cr(
            "ACGT", "A-GT", 1, 3)
        self.assertEqual(expanded, {1, 2, 3})
        self.assertEqual(anchor, {1, 2, 3})
        self.assertEqual((n, h, sp), (4, 4, 3))

    def test_trailing_and_leading_reference_gaps_within_circular_cr(self):
        expanded, anchor, n, h, sp = mapped_human_cr(
            "--ACGT--", "TTACGTTT", 1, 8)
        self.assertEqual(len(expanded), 8)
        self.assertEqual(len(anchor), 4)

    def test_reference_boundary_gap_outside_cr(self):
        # 576 CR human positions, followed by a human gap insertion
        # between reference 576 and 577, which must be excluded.
        human = "A" * 576 + "--" + "C"
        species = "A" * 576 + "GG" + "C"
        expanded, anchor, n, h, sp = mapped_human_cr(
            human, species, 1, len(species))
        self.assertEqual(len(anchor), 576)
        self.assertEqual(expanded, anchor)
        self.assertEqual((n, h, sp), (576, 577, len(species)))

    def test_ref_start_boundary_gap_outside_cr(self):
        human = "A" * 16023 + "--" + "C"
        species = "A" * 16023 + "GG" + "C"
        expanded, anchor, n, h, sp = mapped_human_cr(
            human, species, 1, len(species))
        self.assertNotIn(16024, expanded)
        self.assertNotIn(16025, expanded)
        self.assertIn(16026, expanded)


if __name__ == "__main__":
    unittest.main()
