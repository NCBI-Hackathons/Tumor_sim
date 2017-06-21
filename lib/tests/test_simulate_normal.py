import unittest
from .. import simulate_endToEnd
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna

class TestSimulateNormal(unittest.TestCase):

    def test_remove_Ns(self):
        genome = {"chr1": MutableSeq("NNNNAGAGCTACGATGCTACGATGNNNNN", generic_dna),
                  "chr2": MutableSeq("NNNNAGAGCTACNNNGATGCGATGNN", generic_dna)}
        genome_out = simulate_endToEnd.remove_trailing_N_characters(genome)
        self.assertEqual(genome_out, {"chr1": MutableSeq("AGAGCTACGATGCTACGATG", generic_dna),
                                      "chr2": MutableSeq("AGAGCTACNNNGATGCGATG", generic_dna)})
