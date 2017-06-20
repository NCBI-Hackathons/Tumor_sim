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
        self.assertEqual(self.mo.pick_chromosomes(self.genome)[0], 'chr1')

    def test_generate_structural_variations(self):
        number = 8 
        self.mo.generate_structural_variations(self.genome, number)

    def test_get_location_on_sequence(self):
        self.genome = {'chr1': MutableSeq("NANNNNNNNNN", generic_dna)}
        location = self.mo.get_location_on_sequence(self.genome['chr1'], 'uniform')
        self.assertEqual(location, 1)

    def test_pick_chromosomes(self):
        return
