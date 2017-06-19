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
