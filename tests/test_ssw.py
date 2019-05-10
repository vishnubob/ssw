#!/usr/bin/env python

import unittest
import pickle
import ssw
from ssw import iupac

class TestIUPAC(unittest.TestCase):
    def test_degen_revcomp(self):
        test_seq = "AGTCMRYSWKBDHVN"
        test_seq_rc = "NBDHVMWSRYKGACT"
        _test_seq_rc = iupac.nucleotide_reverse_complement(test_seq)
        self.assertEqual(test_seq_rc, _test_seq_rc)
        _test_seq = iupac.nucleotide_reverse_complement(test_seq_rc)
        self.assertEqual(test_seq, _test_seq)

class TestPickle(unittest.TestCase):
    def test_alignment_pickle(self):
        reference = "GTGCGATGTGCGATGAGATC"
        query = reference
        aligner = ssw.Aligner()
        al = aligner.align(query, reference)
        orig_dict = al.__dict__.copy()
        clone = pickle.loads(pickle.dumps(al))
        clone_dict = clone.__dict__.copy()
        for key in orig_dict.keys():
            self.assertIn(key, clone_dict)
            self.assertEqual(orig_dict[key], clone_dict[key])

class TestAlignment(unittest.TestCase):
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

    def test_coverage(self):
        bplen = 5
        reference = "GTGCGATGTGCGATGAGATC"
        query = reference[bplen:bplen*2]
        aligner = ssw.Aligner()
        al = aligner.align(query, reference)
        self.assertEqual(al.match_count, bplen)
        self.assertEqual(al.mismatch_count, 0)
        self.assertEqual(al.insertion_count, 0)
        self.assertEqual(al.deletion_count, 0)
        self.assertEqual(al.cigar, '%dM' % bplen)
        self.assertEqual(al.query_coverage, 1.0)
        self.assertEqual(al.reference_coverage, bplen / len(reference))

    def test_degen_alignment(self):
        # XXX: note, this fails if seq_1 and seq_2 are switched
        # see: https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/issues/63
        seq_1 = "AGCGATCACGT"
        seq_2 = "MRYSWKBDHVN"
        aligner = ssw.Aligner()
        al = aligner.align(seq_1, seq_2)
        self.assertEqual(al.match_count, len(seq_1))
        self.assertEqual(al.mismatch_count, 0)
        self.assertEqual(al.insertion_count, 0)
        self.assertEqual(al.deletion_count, 0)
        self.assertEqual(al.cigar, '%dM' % len(seq_1))

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
