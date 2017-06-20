import unittest
from .. import mutation_orchestrator
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna

class TestMutationOrchestrator(unittest.TestCase):

    def setUp(self):
        self.genome = {'chr1': MutableSeq("ACTCGTCGTC", generic_dna),
        'chr2': MutableSeq("ACTCGTCGTC", generic_dna)}
        self.mo = mutation_orchestrator.Mutation_Orchestrator()

    def test_deletion(self):
        return

    def test_pick_chromosome_always_pick_non_empty_chrom(self):
        self.genome['chr2'] = MutableSeq("", generic_dna)
        self.assertEqual(self.mo.pick_chromosome(self.genome), 'chr1')

    def test_generate_structural_variations(self):
        number = 8 
        self.mo.generate_structural_variations(self.genome, number)
