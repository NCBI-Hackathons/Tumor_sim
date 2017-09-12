import unittest
from .. import mutation_orchestrator
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna
import pandas as pd
import copy

class Dummy_Mutation_Orchestrator(mutation_orchestrator.Mutation_Orchestrator):

    def get_event_length(self, p=1):
        return 2

    def pick_chromosomes(self, genome, number = 1, replace=True):
        if replace:
            return ['chr2']
        else:
            return(['chr2', 'chr1'])

    def get_location_on_sequence(self, seq):
        return 2

class TestMutationOrchestrator(unittest.TestCase):

    def setUp(self):
        self.genome = {'chr1': MutableSeq("TCGTCGTC", generic_dna),
        'chr2': MutableSeq("ACTCGTCGTC", generic_dna)}
        self.mo = mutation_orchestrator.Mutation_Orchestrator()

    def test_duplication(self):
        self.mo = Dummy_Mutation_Orchestrator()
        self.mo.orchestrate_duplication(self.genome)
        original_genome = copy.deepcopy(self.genome)
        (returned_genome, bed) = self.mo.tracker.collapse_list(self.genome)
        expected_name = 'duplication (times 2)'
        realized_name = bed['name'][0]
        self.assertEqual(expected_name, realized_name)
        expected_alt = 'TCTC'
        realized_alt = bed['alt'][0]
        self.assertEqual(expected_alt, realized_alt)
        # Assert that the insertion increased the length of the genome at chr2 by 4
        self.assertEqual(len(original_genome['chr2']) + 4, len(returned_genome['chr2']))

    def test_pick_chromosome_always_pick_non_empty_chrom(self):
        self.genome['chr2'] = MutableSeq("", generic_dna)
        self.assertEqual(self.mo.pick_chromosomes(self.genome)[0], 'chr1')

    def test_generate_structural_variations(self):
        number = 8 
        s = self.mo.generate_structural_variations(self.genome, number)
        self.assertEqual(None, s)

    def test_generate_structural_variations(self):
        self.mo = Dummy_Mutation_Orchestrator()
        number = 8 
        s = self.mo.generate_indels(self.genome, number)
        self.assertEqual(None, s)
        (gg, bed) = self.mo.generate_fasta_and_bed(self.genome)

    def test_get_location_on_sequence(self):
        self.genome = {'chr1': MutableSeq("NANNNNNNNNN", generic_dna)}
        location = self.mo.get_location_on_sequence(self.genome['chr1'], 'uniform')
        self.assertEqual(location, 1)

    def test_get_location_on_sequence_at_end_of_Ns(self):
        self.genome = {'chr1': MutableSeq("NNNNNNNNNNA", generic_dna)}
        location = self.mo.get_location_on_sequence(self.genome['chr1'], 'uniform')
        self.assertEqual(location, 10)

    def test_bed_correct(self):
        lists = [['chr2', 6, 6, 'insertion', 'AAA', 2],['chr1', 6, 15, 'inversion', '-', 0]]
        df = pd.DataFrame(lists)
        df.columns = ['chrom', 'start', 'end', 'name', 'alt', 'uid']
        lists = [['chr2', 6, 7, 'insertion', 'AAA', 2],['chr1', 6, 15, 'inversion', '-', 0]]
        expected_df = pd.DataFrame(lists)
        expected_df.columns = ['chrom', 'start', 'end', 'name', 'alt', 'uid']
        new_bed = self.mo.bed_correct(df)
        self.assertTrue(expected_df.equals(new_bed))

    def test_orchestrate_deletion(self):
        self.mo = Dummy_Mutation_Orchestrator()
        self.mo.orchestrate_deletion(self.genome)
        original_genome = copy.deepcopy(self.genome)
        (returned_genome, bed) = self.mo.tracker.collapse_list(self.genome)
        expected_name = 'deletion'
        realized_name = bed['name'][0]
        expected_bed = [['chr2', 2, 4, 'deletion', '-',  0]]
        df = pd.DataFrame(expected_bed)
        df.columns = ['chrom', 'start', 'end', 'name', 'alt', 'uid']
        self.assertTrue(bed.equals(df))
        # Assert that the insertion increased the length of the genome at chr2 by 4
        self.assertEqual(len(original_genome['chr2']) - 2 , len(returned_genome['chr2']))


    # Most complex test: validates that a single translocation creates the equivalent
    # of 2 insertions and 2 deletions
    def test_orchestrate_translocation(self):
        self.mo = Dummy_Mutation_Orchestrator()
        self.mo.orchestrate_translocation(self.genome)
        original_genome = copy.deepcopy(self.genome)
        (returned_genome, bed) = self.mo.tracker.collapse_list(self.genome)
        expected_name = 'deletion'
        realized_name = bed['name'][0]
        start = 2
        end = 4
        expected_bed = [['chr2', start, end, 'translocation(del)', '-',  0],
        ['chr1', start, end, 'translocation(del)', '-',  1],
        ['chr2', start, start, 'translocation(chr1:2-4)', 'GT',  2],
        ['chr1', start, start, 'translocation(chr2:2-4)', 'TC',  3]
        ]
        df = pd.DataFrame(expected_bed)
        df.columns = ['chrom', 'start', 'end', 'name', 'alt', 'uid']
        df = df.sort_values('uid')
        bed = bed.sort_values('uid')
        df.index = bed.index
        self.assertTrue(bed.equals(df))
    
        # Assert that the insertion increased the length of the genome at chr2 by 4
        self.assertEqual(len(original_genome['chr2']) , len(returned_genome['chr2']))
        expected_chr1 = original_genome['chr1'][:start] + original_genome['chr2'][start:end] + original_genome['chr1'][end:]
        self.assertEqual(expected_chr1, returned_genome['chr1'])
        expected_chr2 = original_genome['chr2'][:start] + original_genome['chr1'][start:end] + original_genome['chr2'][end:]
        self.assertEqual(expected_chr2, returned_genome['chr2'])


    def test_orchestrate_insertion(self):
        self.mo = Dummy_Mutation_Orchestrator()
        self.mo.orchestrate_insertion(self.genome)
        original_genome = copy.deepcopy(self.genome)
        (returned_genome, bed) = self.mo.tracker.collapse_list(self.genome)
        expected_name = 'deletion'
        realized_name = bed['name'][0]
        expected_bed = [['chr2', 2, 2, 'insertion', 'TC',  0]]
        df = pd.DataFrame(expected_bed)
        df.columns = ['chrom', 'start', 'end', 'name', 'alt', 'uid']
        print bed
        self.assertTrue(bed.equals(df))
        # Assert that the insertion increased the length of the genome at chr2 by 4
        self.assertEqual(len(original_genome['chr2']) + 2 , len(returned_genome['chr2']))

    def test_orchestrate_inversion(self):
        self.mo = Dummy_Mutation_Orchestrator()
        self.mo.orchestrate_inversion(self.genome)
        original_genome = copy.deepcopy(self.genome)
        (returned_genome, bed) = self.mo.tracker.collapse_list(self.genome)
        expected_name = 'deletion'
        realized_name = bed['name'][0]
        expected_bed = [['chr2', 2, 4, 'inversion', '-',  0]]
        df = pd.DataFrame(expected_bed)
        df.columns = ['chrom', 'start', 'end', 'name', 'alt', 'uid']
        self.assertTrue(bed.equals(df))
        # Assert that the inversion kept the length of the genome at chr2 the same
        self.assertEqual(len(original_genome['chr2']), len(returned_genome['chr2']))


