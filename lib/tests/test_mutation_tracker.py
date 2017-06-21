import unittest
from .. import mutation_orchestrator
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna
import numpy as np
import pandas as pd

class TestMutationTracker(unittest.TestCase):
    def setUp(self):
        self.genome = {'chr1': MutableSeq("ACTCGTCGTC", generic_dna),
        'chr2': MutableSeq("ACTCGTCGTC", generic_dna)}
        self.mc = mutation_orchestrator.Mutation_Tracker()

    def test_deletion_adds_to_list(self):
        self.mc.create_deletion('chr1', start=5, end=10)
        self.assertEqual(self.mc.list, [['chr1', 5, 10, 'deletion', '-', 0]])

    def test_deletion_adds_to_function_dict(self):
        self.mc.create_deletion('chr1', start=5, end=10)
        new_genome = self.mc.collapse_list(self.genome)
        expected_chr1 = MutableSeq("ACTCG", generic_dna)
        self.assertEqual(new_genome['chr1'], expected_chr1)

    def test_two_deletion_adds_to_function_dict(self):
        self.mc.create_deletion('chr1', start=5, end=6)
        self.mc.create_deletion('chr1', start=8, end=10)
        new_genome = self.mc.collapse_list(self.genome)
        expected_chr1 = MutableSeq("ACTCGCG", generic_dna)
        self.assertEqual(new_genome['chr1'], expected_chr1)

    def test_two_deletion_added_to_two_chromosomes(self):
        self.mc.create_deletion('chr1', start=5, end=9)
        self.mc.create_deletion('chr2', start=1, end=10)
        new_genome = self.mc.collapse_list(self.genome)
        expected_chr1 = MutableSeq("ACTCGC", generic_dna)
        expected_chr2 = MutableSeq("A", generic_dna)
        self.assertEqual(new_genome['chr1'], expected_chr1)
        self.assertEqual(new_genome['chr2'], expected_chr2)

    def test_inversion(self):
        self.mc.create_inversion('chr1', start=1, end=3)
        new_genome = self.mc.collapse_list(self.genome)
        self.assertEqual(new_genome['chr1'], MutableSeq("ATCCGTCGTC", generic_dna))

    def test_inversion_outside_range(self):
        self.mc.create_inversion('chr1', start=6, end=15)
        new_genome = self.mc.collapse_list(self.genome)
        self.assertEqual(new_genome['chr1'], MutableSeq("ACTCGTCTGC", generic_dna))

    def test_insertion(self):
        self.mc.create_insertion('chr1', start=4, new_seq='GGAA')
        new_genome = self.mc.collapse_list(self.genome)
        self.assertEqual(new_genome['chr1'], MutableSeq("ACTCGGAAGTCGTC", generic_dna))

    # Start: ACTCGTCGTC
    # Inversion: ACTCGT - CTGC
    # Deletion: ACTC - TGC

    # Switch order:
    # Start: ACTCGTCGTC
    # Deletion: ACTC - GTC
    # Inversion:ACTCGTC
    # DIFFERENT RESULT!

    def test_multiple_structural_variations(self):
        self.mc.create_inversion('chr1', start=6, end=15)
        self.mc.create_deletion('chr1', start=4, end=7)
        self.mc.create_insertion('chr2', start=6, new_seq = 'AAA')
        new_genome = self.mc.collapse_list(self.genome)
        expected_genome = {'chr1': MutableSeq("ACTCGTCTGC", generic_dna),
                           'chr2': MutableSeq("ACTCGTAAACGTC", generic_dna)}
        self.assertEqual(new_genome, expected_genome)
        lists = [['chr2', 6, 7, 'insertion', 'AAA', 2],['chr1', 6, 15, 'inversion', '-', 0]]
        expected_df = pd.DataFrame(lists)
        expected_df.columns = ['chrom', 'start', 'end', 'name', 'ALT', 'uid']
        expected_df.index = [2,0]
        
        self.assertTrue(expected_df.equals(self.mc.log_data_frame))


    def tst_deletion_outside_range(self):
        my_seq = MutableSeq("AAAAAA", generic_dna)
        mc = mutation_creator.Mutation_Creator()
        new_seq = mc.create_deletion(my_seq, start=5, end=10)
        self.assertEqual(new_seq, MutableSeq("AAAAA", generic_dna))

    def tst_inversion(self):
        my_seq = MutableSeq("ACTCGTCGTC", generic_dna)
        mc = mutation_creator.Mutation_Creator()
        new_seq = mc.create_inversion(my_seq, start=3, end=7)
        self.assertEqual(new_seq, MutableSeq("ACTCTGCGTC", generic_dna))
