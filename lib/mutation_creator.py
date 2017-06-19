from Bio import SeqIO

class Mutation_Creator:
    def create_deletion(self, mutable_seq, start, end):
        return mutable_seq

    def create_snv(self, mutable_seq, start, end, new_base):
        return mutable_seq

    def create_insertion(self, mutable_seq, start, end, new_seq):
        return mutable_seq

    def create_inversion(self, mutable_seq, start, end, new_start, new_end):
        return mutable_seq