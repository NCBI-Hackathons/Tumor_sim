import unittest
from .. import mutation_tracker
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna


class TestMutationTracker(unittest.TestCase):
    def setUp(self):
        self.genome = {'chr1': MutableSeq("ACTCGTCGTCGCAAC", generic_dna),
        'chr2': MutableSeq("ACTCGTCGTC", generic_dna)}
        self.tracker = mutation_tracker.Mutation_Tracker()

    def test_deletion(self):
        self.tracker.create_structure(self.genome)
        self.assertEqual(self.tracker.track_deletion('chr1', 11, 13), [{'chr1': (0, 10)}, {'chr1': (13, 15)}])

    def test_multiple_deletions_til_the_end(self):
        self.tracker.create_structure(self.genome)
        self.tracker.structure['chr1'] = [{'chr1': (0, 10)}, {'chr1': (13, 15)}]
        self.assertEqual(self.tracker.track_deletion('chr1', 11, 13), [{'chr1': (0, 10)}])

    def test_multiple_deletions(self):
        self.tracker.create_structure(self.genome)
        self.tracker.structure['chr1'] = [{'chr1': (0, 10)}, {'chr1': (13, 15)}]
        self.assertEqual(self.tracker.track_deletion('chr1', 10, 12), [{'chr1': (0, 9)}, {'chr1': (14, 15)}])

