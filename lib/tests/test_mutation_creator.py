import unittest
from .. import mutation_creator
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna

class TestMutationCreator(unittest.TestCase):

    def test_deletion(self):
        my_seq = MutableSeq("ACTCGTCGTC", generic_dna)
        mc = mutation_creator.Mutation_Creator()
        new_seq = mc.create_deletion(my_seq, start=5, end=10)
        self.assertEqual(new_seq, MutableSeq("ACTCG", generic_dna))

    def test_inversion(self):
        my_seq = MutableSeq("ACTCGTCGTC", generic_dna)
        mc = mutation_creator.Mutation_Creator()
        new_seq = mc.create_inversion(my_seq, start=3, end=7)
        self.assertEqual(new_seq, MutableSeq("ACTCTGCGTC", generic_dna))

    def test_snv(self):
        my_seq = MutableSeq("ACTCGTCGTC", generic_dna)
        mc = mutation_creator.Mutation_Creator()
        new_seq = mc.create_snv(my_seq, start=3, new_base='A')
        self.assertEqual(new_seq, MutableSeq("ACTAGTCGTC", generic_dna))

    def test_snv_works_on_N(self):
        my_seq = MutableSeq("AN", generic_dna)
        mc = mutation_creator.Mutation_Creator()
        new_seq = mc.create_snv(my_seq, start=1, new_base='A')
        self.assertEqual(new_seq, MutableSeq("AA", generic_dna))

    def test_snv_works_can_insert_N(self):
        my_seq = MutableSeq("AA", generic_dna)
        mc = mutation_creator.Mutation_Creator()
        new_seq = mc.create_snv(my_seq, start=1, new_base='N')
        self.assertEqual(new_seq, MutableSeq("AN", generic_dna))

    def test_insertion(self):
        my_seq = MutableSeq("ACTCGTCGTC", generic_dna)
        mc = mutation_creator.Mutation_Creator()
        new_seq = mc.create_insertion(my_seq, start=3, new_seq='AAAA')
        self.assertEqual(new_seq, MutableSeq("ACTAAAACGTCGTC", generic_dna))

    def test_inversion_larger_than_seq(self):
        my_seq = MutableSeq("CAAAT", generic_dna)
        mc = mutation_creator.Mutation_Creator()
        new_seq = mc.create_inversion(my_seq, start=1, end=7)
        self.assertEqual(new_seq, MutableSeq("CTAAA", generic_dna))

    def test_translocation(self):
        seq1 = MutableSeq("AAAAAAAAA", generic_dna)
        seq2 = MutableSeq("CCCC", generic_dna)
        mc = mutation_creator.Mutation_Creator()
        (new_seq1, new_seq2)  = mc.create_translocation(seq1, seq2, start1=3, start2=2, length1=5, length2=2)
        expected_seq1 = MutableSeq("AAACCA", generic_dna)
        expected_seq2 = MutableSeq("CCAAAAA", generic_dna)
        self.assertEqual(new_seq1, expected_seq1)
        self.assertEqual(new_seq2, expected_seq2)

    def test_translocation_handles_excessive_event_lengths(self):
        seq1 = MutableSeq("AAAAAAAAA", generic_dna)
        seq2 = MutableSeq("CCCC", generic_dna)
        mc = mutation_creator.Mutation_Creator()
        (new_seq1, new_seq2)  = mc.create_translocation(seq1, seq2, start1=3, start2=2, length1=5, length2=5)
        expected_seq1 = MutableSeq("AAACCA", generic_dna)
        expected_seq2 = MutableSeq("CCAAAAA", generic_dna)
        self.assertEqual(new_seq1, expected_seq1)
        self.assertEqual(new_seq2, expected_seq2)

    def test_translocation_handles_length_1_events(self):
        seq1 = MutableSeq("AAAAAAA", generic_dna)
        seq2 = MutableSeq("CCCC", generic_dna)
        mc = mutation_creator.Mutation_Creator()
        (new_seq1, new_seq2)  = mc.create_translocation(seq1, seq2, start1=3, start2=2, length1=1, length2=1)
        expected_seq1 = MutableSeq("AAACAAA", generic_dna)
        expected_seq2 = MutableSeq("CCAC", generic_dna)
        self.assertEqual(new_seq1, expected_seq1)
        self.assertEqual(new_seq2, expected_seq2)

