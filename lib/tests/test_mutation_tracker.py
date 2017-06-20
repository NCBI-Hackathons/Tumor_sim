import unittest
from .. import mutation_tracker

class TestMutationTracker(unittest.TestCase):
    def setUp(self):
        self.genome = {'chr1': MutableSeq("ACTCGTCGTCGCAAC", generic_dna),
        'chr2': MutableSeq("ACTCGTCGTC", generic_dna)}
        self.tracker = mutation_tracker.Mutation_Tracker()

    def test_deletion(self):
        self.tracker.create_structure(self.genome)
        self.tracker.track_deletion('chr1', 11, 13)
        print(self.tracker.structure)
