import unittest
from .. import simulate_endToEnd
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna
import pandas as pd

class TestSimulateNormal(unittest.TestCase):

    def test_remove_Ns(self):
        genome = {"chr1": MutableSeq("NNNNAGAGCTACGATGCTACGATGNNNNN", generic_dna),
                  "chr2": MutableSeq("NNNNAGAGCTACNNNGATGCGATGNN", generic_dna)}
        genome_out = {}
        (genome_out['chr1'], offset) = simulate_endToEnd.remove_trailing_N_characters(genome['chr1'])
        (genome_out['chr2'], offset) = simulate_endToEnd.remove_trailing_N_characters(genome['chr2'])
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