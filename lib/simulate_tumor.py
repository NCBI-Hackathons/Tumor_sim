from Bio import SeqIO
from mutation_orchestrator import Mutation_Orchestrator
import copy

normal_fasta_file = "tests/output.fasta" ### PATH
output_fasta_file = "tests/tumorsim.fasta"

def write_fasta(mutated_genome):  
    output_seqs = []
    for chrom in mutated_genome:
        output_seq = SeqIO.SeqRecord(seq=mutated_genome[chrom].toseq(),
              id=chrom, description='', name='')
        output_seqs.append(output_seq)

    with open(output_fasta_file, "w") as output_handle:
        SeqIO.write(output_seqs, output_handle, "fasta")


def read_fasta():
    original_genome = {} 
    ### takes some time to load entire 3GB hg38 into memory; possible performance problem
    for seq_record in SeqIO.parse(input_fasta_file, "fasta"):
        original_genome[seq_record.id] = seq_record.upper() ## make all characters upper-case
    ## remove all 'useless' chromosomes, i.e. must match chrX, chrY or "^[a-z]{3}\d{1,2}$"
    original_genome = {k: v for k, v in original_genome.items() if re.match('^[a-z]{3}\d{1,2}$', k, re.IGNORECASE) or k in ["chrX", "chrY"]}