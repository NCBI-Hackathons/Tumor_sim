from Bio import SeqIO
from random import randint

class Mutation_Creator:
    def create_deletion(self, mutable_seq, start, end):
        return mutable_seq[:start] + mutable_seq[end:]

    ### InDels in healthy genomes havee ~400-600K InDels, ~90% sized 1-9 bp, ~10 10-20 bp, sized up to 10Kb 
    
    def create_deletion(self, mutable_seq, start, end):
        """
        For each chromosome, create a range(start_of_mappable_sequence, end_of_mappable_sequence)
        Check for base != 'N', two bases away, then create range(start, end)
        Draw 200-300 positions in range(start, end)
            For each position, draw 1-9 90% of the time, 10-20 10% of the time. Delete at end = start + draw
        random.randrange(11) ### 1, 2, ..., 10; if 10, then range(
        randint(start_mappable, end_mappable)
        """
        
      

    def create_snv(self, mutable_seq, start, new_base):
        mutable_seq[start] = new_base
        return mutable_seq

    def create_insertion(self, mutable_seq, start, new_seq):
        return mutable_seq[:start] + new_seq +  mutable_seq[start:]

    def create_inversion(self, mutable_seq, start, end):
        inv_seq = mutable_seq[start:end]
        inv_seq.reverse()
        mutable_seq[start:end] = inv_seq
        return mutable_seq

    def create_translocation(self, seq1, seq2, start1, start2, length1, length2):
        new_seq1 = self.create_deletion(seq1, start1,start1+length1)
        new_seq1 = self.create_insertion(new_seq1, start1, seq2[start2:start2+length2])
        new_seq2 = self.create_deletion(seq2, start2,start2+length2)
        new_seq2 = self.create_insertion(new_seq2, start2, seq1[start1:start1+length1])
        return (new_seq1, new_seq2)
