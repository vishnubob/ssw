#!/usr/bin/env python

import unittest
import ssw

class TestSSW(unittest.TestCase):
    def test_mixed_case(self):
        reference = "GTGCGATGTGCGATGAGATC"
        query = reference.lower()
        aligner = ssw.Aligner()
        al = aligner.align(query, reference)
        self.assertEquals(al.match_count, len(reference))
        self.assertEquals(al.mismatch_count, 0)
        self.assertEquals(al.insertion_count, 0)
        self.assertEquals(al.deletion_count, 0)
        self.assertEquals(al.cigar, '20M')

    def test_perfect_alignment(self):
        reference = "GTGCGATGTGCGATGAGATC"
        query = reference
        aligner = ssw.Aligner()
        al = aligner.align(query, reference)
        self.assertEquals(al.match_count, len(reference))
        self.assertEquals(al.mismatch_count, 0)
        self.assertEquals(al.insertion_count, 0)
        self.assertEquals(al.deletion_count, 0)
        self.assertEquals(al.cigar, '20M')

    def test_rc_alignment(self):
        reference = "GTGCGATGTGCGATGAGATC"
        query = "GATCTCATCGCACATCGCAC"
        aligner = ssw.Aligner()
        al = aligner.align(query, reference)
        self.assertEquals(al.match_count, 20)
        self.assertEquals(al.mismatch_count, 0)
        self.assertEquals(al.insertion_count, 0)
        self.assertEquals(al.deletion_count, 0)
        self.assertEquals(al.cigar, "20M")

    def test_insertion(self):
        reference = "GTGCGATGTGCGATGAGATC"
        query = reference[:10] + 'A' + reference[10:]
        aligner = ssw.Aligner()
        al = aligner.align(query, reference)
        self.assertEquals(al.match_count, 20)
        self.assertEquals(al.mismatch_count, 0)
        self.assertEquals(al.insertion_count, 1)
        self.assertEquals(al.deletion_count, 0)
        self.assertEquals(al.cigar, "10M1I10M")

    def test_mismatch(self):
        reference = "GTGCGATGTGCGATGAGATC"
        query = reference[:9] + 'A' + reference[10:]
        aligner = ssw.Aligner()
        al = aligner.align(query, reference)
        self.assertEquals(al.match_count, 19)
        self.assertEquals(al.mismatch_count, 1)
        self.assertEquals(al.insertion_count, 0)
        self.assertEquals(al.deletion_count, 0)
        self.assertEquals(al.cigar, "20M")

    def test_deletion(self):
        reference = "GTGCGATGTGCGATGAGATC"
        query = reference[:10] + reference[11:]
        aligner = ssw.Aligner()
        al = aligner.align(query, reference)
        self.assertEquals(al.match_count, 19)
        self.assertEquals(al.mismatch_count, 0)
        self.assertEquals(al.insertion_count, 0)
        self.assertEquals(al.deletion_count, 1)
        self.assertEquals(al.cigar, "10M1D9M")

if __name__ == '__main__':
    unittest.main()
