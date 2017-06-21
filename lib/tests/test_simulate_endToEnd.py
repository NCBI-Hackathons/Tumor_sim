import unittest
from .. import simulate_endToEnd
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna
import pandas as pd

class TestSimulateNormal(unittest.TestCase):

    def setUp(self):
        self.genome = {"chr1": MutableSeq("NNNNAGAGCTACGATGCTACGATGNNNNN", generic_dna),
                  "chr2": MutableSeq("NNNNNNAGAGCTACNNNGATGCGATGNN", generic_dna)}

    def test_remove_Ns(self):
        genome_out = {}
        (genome_out['chr1'], offset) = simulate_endToEnd.remove_trailing_N_characters(self.genome['chr1'])
        (genome_out['chr2'], offset) = simulate_endToEnd.remove_trailing_N_characters(self.genome['chr2'])
        self.assertEqual(genome_out, {"chr1": MutableSeq("AGAGCTACGATGCTACGATG", generic_dna),
                                      "chr2": MutableSeq("AGAGCTACNNNGATGCGATG", generic_dna)})

    def test_subtract_beds(self):
        lists = [['chr2', 6, 7, 'insertion', 'AAA', 2],['chr1', 6, 15, 'inversion', '-', 0]]
        first_bed = pd.DataFrame(lists)
        first_bed.columns = ['chrom', 'start', 'end', 'name', 'alt', 'uid']

        lists = [['chr2', 6, 7, 'insertion', 'AAA', 2]]
        second_bed = pd.DataFrame(lists)
        second_bed.columns = ['chrom', 'start', 'end', 'name', 'alt', 'uid']

        new_bed = simulate_endToEnd.subtract_beds(first_bed, second_bed)
        # Have to reset the index, or otherwise the indices will be unequal
        new_bed = new_bed.reset_index(drop=True)

        lists = [['chr1', 6, 15, 'inversion', '-', 0]]
        expected_df = pd.DataFrame(lists)
        expected_df.columns = ['chrom', 'start', 'end', 'name', 'alt', 'uid']
        self.assertTrue(expected_df.equals(new_bed))

    def test_genome_offset(self):
        genome_out = {}
        genome_offset = {}
        for chrom in self.genome:
            (genome_out[chrom], genome_offset[chrom]) = simulate_endToEnd.remove_trailing_N_characters(
                self.genome[chrom])
        self.assertEqual(genome_offset['chr1'], 4)
        self.assertEqual(genome_offset['chr2'], 6)

    def test_bed_reoffset(self):
        genome_out = {}
        genome_offset = {}
        for chrom in self.genome:
            (genome_out[chrom], genome_offset[chrom]) = simulate_endToEnd.remove_trailing_N_characters(
                self.genome[chrom])

        lists = [['chr2', 6, 7, 'insertion', 'AAA', 2],['chr1', 6, 15, 'inversion', '-', 0]]
        first_bed = pd.DataFrame(lists)
        first_bed.columns = ['chrom', 'start', 'end', 'name', 'alt', 'uid']

        corrected_bed = simulate_endToEnd.offset_bed(first_bed, genome_offset)
        lists = [['chr2', 12, 13, 'insertion', 'AAA', 2],['chr1', 10, 19, 'inversion', '-', 0]]
        expected_df = pd.DataFrame(lists)
        expected_df.columns = ['chrom', 'start', 'end', 'name', 'alt', 'uid']
        self.assertTrue(expected_df.equals(corrected_bed))



