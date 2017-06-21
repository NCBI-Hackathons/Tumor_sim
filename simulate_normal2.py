from Bio import SeqIO
from mutation_orchestrator import Mutation_Orchestrator
import copy

input_fasta_file = '../data/subsampled_hg38.fa'
number_snvs = 400000
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

def count_trailing_N_characters(sequence):
    startc = 0
    lastc = 0
    for chrom in sequence:
        while sequence[chrom][0] == 'N':
            startc += 1
        while sequence[chrom][-1] == 'N':
            lastc += 1
    return [startc,lastc]

def remove_trailing_N_characters(sequence):
    for chrom in sequence:
        startc = 0
        while sequence[chrom][0] == 'N':
            startc += 1
            sequence[chrom].pop(0)
        lastc = 0
        while sequence[chrom][-1] == 'N':
            sequence[chrom].pop(-1)
            lastc += 1
    return sequence

def main():
    original_genome = {}
    trailing_N_characters = {}
    mutated_genome = {}
    for seq_record in SeqIO.parse(input_fasta_file, "fasta"):
        original_genome[seq_record.id] = seq_record
        trailing_N_characters[seq_record.id] = count_trailing_N_characters(seq_record.seq.tomutable())
        mutated_genome[seq_record.id] = remove_trailing_N_characters(seq_record.seq.tomutable())

    orchestrator = Mutation_Orchestrator()

    # add SNVs
    mutated_genome = orchestrator.snv_fast(mutated_genome, number_snvs)

    # add structural varations
    # read ge
    mutated_genome = orchestrator.generate_structural_variations(mutated_genome, 10000)

    write_fasta(mutated_genome)
    write_bam(mutated_genome)

if __name__ == "__main__":
    main()
