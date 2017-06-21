from Bio import SeqIO
from mutation_orchestrator import Mutation_Orchestrator
import copy
import re

input_reference_fasta_file = '../data/subsampled_hg38.fa'
number_snvs = 3000000  
number_indels = 415000
number_of_tumorSVs = 10000
output_normal_fasta_file = 'tests/normalsim.fasta'
output_tumor_fasta_file = "tests/tumorsim.fasta'"
output_bedfile = "tests/output.bed"

def write_fasta(genome, output_fasta_file):  
    # write fasta
    output_seqs = []
    for chrom in genome:
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

def read_fasta_normal(input_fasta_file):
    ### takes some time to load entire 3GB hg38 into memory; possible performance problem
    genome = {}
    for seq_record in SeqIO.parse(input_fasta_file, "fasta"):
        genome[seq_record.id] = remove_trailing_N_characters(seq_record.upper().seq.tomutable())
    ## remove all 'useless' chromosomes, i.e. must match chrX, chrY or "^[a-z]{3}\d{1,2}$"
    genome = {k: v for k, v in genome.items() 
                    if re.match('^[a-z]{3}\d{1,2}$', k, re.IGNORECASE) or k in ["chrX", "chrY"]}
    return genome



def main():
    # read genome fasta
    mutated_genome = read_fasta_normal(input_reference_fasta_file)

    orchestrator = Mutation_Orchestrator()
    
    # add germilne SNVs & InDels
    mutated_genome = orchestrator.snv_fast(mutated_genome, number_snvs)
    write_fasta(mutated_genome, output_normal_fasta_file)

    orchestrator.generate_indels(mutated_genome, number_snvs)
    indeled_genome = orchestrator.generate_fasta(mutated_genome)
    write_fasta(indeled_genome, output_normal_fasta_file)
    indel_bed = orchestrator.get_pandas_dataframe
    write_bed(indel_bed)   ### write out "normalsim" bedpe

    # add structural varations
    orchestrator.generate_structural_variations(mutated_genome, number_of_tumorSVs)
    mutated_genome = orchestrator.generate_fasta(mutated_genome)
    write_fasta(mutated_genome, output_tumor_fasta_file)

    write_bed(orchestrator.get_pandas_dataframe)  ### write out "tumorsim" bedpe

if __name__ == "__main__":
    main()

