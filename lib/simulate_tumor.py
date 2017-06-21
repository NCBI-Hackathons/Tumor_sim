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


def read_fasta_tumor(normal_genome, tumor_genome):
    ### takes some time to load entire 3GB hg38 into memory; possible performance problem
    for seq_record in SeqIO.parse(input_fasta_file, "fasta"):
        normal_genome[seq_record.id] = seq_record 
        tumor_genome[seq_record.id] = seq_record 


def main():
    # read genome fasta
    normal_genome = {}
    tumor_genome = {}
 
    read_fasta_tumor(normal_genome, tumor_genome)

    orchestrator = Mutation_Orchestrator()

    # add SNVs
    tumor_genome = orchestrator.snv_fast(tumor_genome, number_snvs)

    # add structural varations
    tumor_genome = orchestrator.generate_structural_variations(tumor_genome, 10000)

    write_fasta(tumor_genome)

if __name__ == "__main__":
    main()