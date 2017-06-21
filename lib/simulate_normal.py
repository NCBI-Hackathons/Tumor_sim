from Bio import SeqIO
from mutation_orchestrator import Mutation_Orchestrator
import copy
import re

input_fasta_file = '../data/subsampled_hg38.fa'
number_snvs = 400000
output_fasta_file = 'tests/normalsim.fasta'


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
    for chrom in sequence:
        while sequence[chrom][0] == 'N':
            sequence[chrom].pop(0)
        while sequence[chrom][-1] == 'N':
            sequence[chrom].pop(-1)
    return sequence

def read_fasta():
    original_genome = {} 
    ### takes some time to load entire 3GB hg38 into memory; possible performance problem
    for seq_record in SeqIO.parse(input_fasta_file, "fasta"):
        original_genome[seq_record.id] = seq_record.upper() ## make all characters upper-case
    ## remove all 'useless' chromosomes, i.e. must match chrX, chrY or "^[a-z]{3}\d{1,2}$"
    original_genome = {k: v for k, v in original_genome.items() if re.match('^[a-z]{3}\d{1,2}$', k, re.IGNORECASE) or k in ["chrX", "chrY"]}
    


def main():
    # read genome fasta
    original_genome = {}
    mutated_genome = {}
    for seq_record in SeqIO.parse(input_fasta_file, "fasta"):
        original_genome[seq_record.id] = seq_record
        mutated_genome[seq_record.id] = remove_trailing_N_characters(seq_record.seq.tomutable())

    orchestrator = Mutation_Orchestrator()

    # add SNVs
    mutated_genome = orchestrator.snv_fast(mutated_genome, number_snvs)

    # add structural varations
    mutated_genome = orchestrator.generate_structural_variations(mutated_genome, 10000)

    write_fasta(mutated_genome)

if __name__ == "__main__":
    main()

