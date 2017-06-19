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

    def test_insertion(self):
        my_seq = MutableSeq("ACTCGTCGTC", generic_dna)
        mc = mutation_creator.Mutation_Creator()
        new_seq = mc.create_insertion(my_seq, start=3, new_seq='AAAA')
        self.assertEqual(new_seq, MutableSeq("ACTAAAACGTCGTC", generic_dna))
