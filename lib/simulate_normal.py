from Bio import SeqIO
from mutation_orchestrator import Mutation_Orchestrator
import copy

input_fasta_file = '../data/subsampled_hg38.fa'
iterations = 100
output_fasta_file = 'tests/output.fasta'

def write_bam(mutated_genome):
    raise NotImplementedError

def write_fasta(mutated_genome):
    # write fasta
    output_seqs = []
    for chrom in mutated_genome:
        output_seq = SeqIO.SeqRecord(seq=mutated_genome[chrom].toseq(),
              id=chrom, description='', name='')
        output_seqs.append(output_seq)

    with open(output_fasta_file, "w") as output_handle:
        SeqIO.write(output_seqs, output_handle, "fasta")

def remove_trailing_N_characters(sequence):
    return sequence

def main():
    # read genome fasta
    original_genome = {}
    mutated_genome = {}
    for seq_record in SeqIO.parse(input_fasta_file, "fasta"):
        original_genome[seq_record.id] = seq_record
        mutated_genome[seq_record.id] = remove_trailing_N_characters(seq_record.seq.tomutable())

    orchestrator = Mutation_Orchestrator()

    # add SNVs
    for i in range(iterations):
        mutated_genome = orchestrator.snv(mutated_genome)

    write_fasta(mutated_genome)
    write_bam(mutated_genome)

if __name__ == "__main__":
    main()
