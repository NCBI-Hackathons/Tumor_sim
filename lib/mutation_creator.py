from Bio import SeqIO

class Mutation_Creator:
    def create_deletion(self, mutable_seq, start, end):
        before_seq = mutable_seq[:start]
        after_seq = mutable_seq[end:]
        before_seq.extend(after_seq)
        return before_seq

    def create_snv(self, mutable_seq, start, new_base):
        mutable_seq[start] = new_base
        return mutable_seq

    def create_insertion(self, mutable_seq, start, new_seq):
        before_seq = mutable_seq[:start]
        after_seq = mutable_seq[start:]
        before_seq.extend(new_seq)
        before_seq.extend(after_seq)
        return before_seq

    def create_inversion(self, mutable_seq, start, end):
        inv_seq = mutable_seq[start:end]
        inv_seq.reverse()
        mutable_seq[start:end] = inv_seq
        return mutable_seq
