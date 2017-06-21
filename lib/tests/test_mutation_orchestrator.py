import unittest
from .. import mutation_orchestrator
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna

class Dummy_Mutation_Orchestrator(mutation_orchestrator.Mutation_Orchestrator):

    def get_event_length(self, p=1):
        return 2

    def pick_chromosomes(self, genome, number = 1):
        return ['chr2']

    def get_location_on_sequence(self, seq):
        return 2

class TestMutationOrchestrator(unittest.TestCase):

    def setUp(self):
        self.genome = {'chr1': MutableSeq("ACTCGTCGTC", generic_dna),
        'chr2': MutableSeq("ACTCGTCGTC", generic_dna)}
        self.mo = mutation_orchestrator.Mutation_Orchestrator()

    def test_duplication(self):
        self.mo = Dummy_Mutation_Orchestrator()
        self.mo.orchestrate_duplication(self.genome)
        self.mo.tracker.collapse_list(self.genome)
        expected_name = 'duplication (times 2)'
        realized_name = self.mo.tracker.log_data_frame['name'][0]
        self.assertEqual(expected_name, realized_name)
        expected_alt = 'TCTC'
        realized_alt = self.mo.tracker.log_data_frame['alt'][0]
        self.assertEqual(expected_alt, realized_alt)

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

    def test_get_location_on_sequence_at_end_of_Ns(self):
        self.genome = {'chr1': MutableSeq("NNNNNNNNNNA", generic_dna)}
        location = self.mo.get_location_on_sequence(self.genome['chr1'], 'uniform')
        self.assertEqual(location, 10)

    def test_orchestrate_translocation(self):
        return

    def test_pick_chromosomes(self):
        return
