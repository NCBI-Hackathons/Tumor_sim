import unittest
from .. import mutation_orchestrator
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna

class TestMutationOrchestrator(unittest.TestCase):

    def test_deletion(self):
        my_seq = MutableSeq("ACTCGTCGTC", generic_dna)
        mutation_creator = create_mutations.Mutation_Creator()
        new_seq = mutation_creator.create_deletion(my_seq, start=5, end=10)
        self.assertEqual(new_seq, MutableSeq("ACTCG", generic_dna))
