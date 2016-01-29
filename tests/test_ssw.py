#!/usr/bin/env python

import unittest
import ssw

class TestSSW(unittest.TestCase):
    def test_mixed_case(self):
        reference = "GTGCGATGTGCGATGAGATC"
        query = reference.lower()
        aligner = ssw.Aligner()
        al = aligner.align(query, reference)
        self.assertEqual(al.match_count, len(reference))
        self.assertEqual(al.mismatch_count, 0)
        self.assertEqual(al.insertion_count, 0)
        self.assertEqual(al.deletion_count, 0)
        self.assertEqual(al.cigar, '20M')

    def test_perfect_alignment(self):
        reference = "GTGCGATGTGCGATGAGATC"
        query = reference
        aligner = ssw.Aligner()
        al = aligner.align(query, reference)
        self.assertEqual(al.match_count, len(reference))
        self.assertEqual(al.mismatch_count, 0)
        self.assertEqual(al.insertion_count, 0)
        self.assertEqual(al.deletion_count, 0)
        self.assertEqual(al.cigar, '20M')

    def test_rc_alignment(self):
        reference = "GTGCGATGTGCGATGAGATC"
        query = "GATCTCATCGCACATCGCAC"
        aligner = ssw.Aligner()
        al = aligner.align(query, reference)
        self.assertEqual(al.match_count, 20)
        self.assertEqual(al.mismatch_count, 0)
        self.assertEqual(al.insertion_count, 0)
        self.assertEqual(al.deletion_count, 0)
        self.assertEqual(al.cigar, "20M")

    def test_insertion(self):
        reference = "GTGCGATGTGCGATGAGATC"
        query = reference[:10] + 'A' + reference[10:]
        aligner = ssw.Aligner()
        al = aligner.align(query, reference)
        self.assertEqual(al.match_count, 20)
        self.assertEqual(al.mismatch_count, 0)
        self.assertEqual(al.insertion_count, 1)
        self.assertEqual(al.deletion_count, 0)
        self.assertEqual(al.cigar, "10M1I10M")

    def test_mismatch(self):
        reference = "GTGCGATGTGCGATGAGATC"
        query = reference[:9] + 'A' + reference[10:]
        aligner = ssw.Aligner()
        al = aligner.align(query, reference)
        self.assertEqual(al.match_count, 19)
        self.assertEqual(al.mismatch_count, 1)
        self.assertEqual(al.insertion_count, 0)
        self.assertEqual(al.deletion_count, 0)
        self.assertEqual(al.cigar, "20M")

    def test_deletion(self):
        reference = "GTGCGATGTGCGATGAGATC"
        query = reference[:10] + reference[11:]
        aligner = ssw.Aligner()
        al = aligner.align(query, reference)
        self.assertEqual(al.match_count, 19)
        self.assertEqual(al.mismatch_count, 0)
        self.assertEqual(al.insertion_count, 0)
        self.assertEqual(al.deletion_count, 1)
        self.assertEqual(al.cigar, "10M1D9M")

    def test_issue_1(self):
        # https://github.com/vishnubob/ssw/issues/1
        reference = "CCC" + "AGCT" * 10
        query = "AGGT" * 10
        aligner = ssw.Aligner()
        alignment = aligner.align(query, reference)
        (r_line, m_line, q_line) = alignment.alignment
        self.assertEqual(r_line, "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG")
        self.assertEqual(m_line, "||*|||*|||*|||*|||*|||*|||*|||*|||*|||")
        self.assertEqual(q_line, "AGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAG")

    def test_issue_3(self):
        # https://github.com/vishnubob/ssw/issues/3
        self.assertRaises(ValueError, ssw.Aligner, gap_open=1, gap_extend=2)
        self.assertRaises(ValueError, ssw.Aligner, gap_open=1, gap_extend=1)

if __name__ == '__main__':
    unittest.main()
