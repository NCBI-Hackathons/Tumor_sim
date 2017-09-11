from Bio import SeqIO
from mutation_orchestrator import Mutation_Orchestrator
import copy
import pandas as pd
import re
import argparse

number_snvs = 3000
number_indels = 4150
number_of_tumorSVs = 10000
output_normal_bedfile = "tests/normal.bed"
output_tumor_bedfile = "tests/tumor.bed"



def write_fasta(genome, output_fasta_file):  
    # write fasta
    output_seqs = []
    for chrom in genome:
        output_seq = SeqIO.SeqRecord(seq=genome[chrom].toseq(),
              id=chrom, description='', name='')
        output_seqs.append(output_seq)

    with open(output_fasta_file, "w") as output_handle:
        SeqIO.write(output_seqs, output_handle, "fasta")


def remove_trailing_N_characters(sequence):
    start_index = len(str(sequence)) - len(str(sequence).lstrip("N"))
    end_index = len(str(sequence).rstrip("N")) #- len(sequence) - 1
    sequence = sequence[start_index:end_index]
    offset = start_index
    return (sequence, offset)


def read_fasta_normal(input_fasta_file):
    ### takes some time to load entire 3GB hg38 into memory; possible performance problem
    genome = {}
    genome_offset = {}
    for seq_record in SeqIO.parse(input_fasta_file, "fasta"):
        genome[seq_record.id], genome_offset[seq_record.id] = remove_trailing_N_characters(
            seq_record.upper().seq.tomutable())

    ## remove all 'useless' chromosomes, i.e. must match chrX, chrY or "^[a-z]{3}\d{1,2}$"
    genome = {k: v for k, v in genome.items() 
                    if re.match('^[a-z]{3}\d{1,2}$', k, re.IGNORECASE) or k in ["chrX", "chrY"]}
    return (genome, genome_offset)

def subtract_beds(bed1, bed2):
    return bed1[~(bed1['uid'].isin(bed2['uid']))]

def write_bed(genome_offset, dframe, path):
    corrected_bed = offset_bed(dframe, genome_offset)
    corrected_bed.to_csv(path, index=False)

def offset_bed(df, genome_offset):
    for chrom in genome_offset:
        per_chrom = df[df['chrom'] == chrom]
        df.ix[per_chrom.index, 'end'] += genome_offset[chrom]
        df.ix[per_chrom.index, 'start'] += genome_offset[chrom]
    return df


def create_complementary_genome(genome):
    ### all reference FASTAs written in 5'-3'
    ### this function outputs the 3'-5' complement
    for chrom in genome:
        genome[chrom] = genome[chrom].toseq().complement()
    return genome



def main(args):
    input_reference_fasta_file = args['input_fasta']
    output_tumor_fasta_file = args['output_tumor_fasta']
    output_normal_fasta_file = args['output_normal_fasta']
    outpur_overlapping_normal_fasta_file = args['output_overlapping_normal_fasta']
    output_complement_tumor_fasta_file = args['output_complement_tumor_fasta']
    output_complement_normal_fasta_file = args['output_complement_normal_fasta']
    # read genome fasta
    (mutated_genome, genome_offset) = read_fasta_normal(input_reference_fasta_file)

    orchestrator = Mutation_Orchestrator()
    
    # add germilne SNVs & InDels
    mutated_genome = orchestrator.snv_fast(mutated_genome, number_snvs)
    write_fasta(mutated_genome, output_normal_fasta_file)

    ## output complement 3'-5' strand normal
    normal_complement = create_complementary_genome(mutated_genome)
    write_fasta(normal_complement, output_complement_normal_fasta_file)

    orchestrator.generate_indels(mutated_genome, number_snvs)
    indeled_genome = orchestrator.generate_fasta(mutated_genome)
    write_fasta(indeled_genome, output_normal_fasta_file)

    overlapping_indeled_genome = orchestrator.generate_overlapping_fasta(mutated_genome)
    write_fasta(overlapping_indeled_genome, output_normal_fasta_file)

    indel_bed = orchestrator.get_pandas_dataframe()
    write_bed(genome_offset, indel_bed, output_normal_bedfile)   ### write out "normalsim" bedpe

    # add structural varations
    orchestrator.generate_structural_variations(mutated_genome, number_of_tumorSVs)
    mutated_genome = orchestrator.generate_fasta(mutated_genome)
    write_fasta(mutated_genome, output_tumor_fasta_file)

    ## output complement 3'-5' strand tumor
    tumor_complement = create_complementary_genome(mutated_genome)
    write_fasta(tumor_complement, output_complement_tumor_fasta_file)

    tumor_bed = orchestrator.get_pandas_dataframe()
    tumor_bed = subtract_beds(tumor_bed, indel_bed)
    write_bed(genome_offset, tumor_bed, output_tumor_bedfile)  ### write out "tumorsim" bedpe


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Simulate cancer genomic structural variations')
    parser.add_argument('--input_fasta', 
                        default="../data/subsampled_hg38.fa",
                        help='file path for the input (default genome) fasta')
    parser.add_argument('--output_tumor_fasta',
                        default="tests/tumorsim.fasta",
                        help='file path for the output tumor (cancer genome) fasta')
    parser.add_argument('--output_normal_fasta',
                        default = "tests/normalsim.fasta",
                        help='file path for the output normal (SNV-added) fasta')
    parser.add_argument('--output_overlapping_normal_fasta',
                        default = "tests/normalsim_overlaps.fasta",
                        help='file path for the output normal (SNV-added) fasta; TEMPORARY')
    parser.add_argument('--output_complement_tumor_fasta',
                        default="tests/complement_tumorsim.fasta",
                        help='file path for the output complement 3-5 strand tumor (cancer genome) fasta')
    parser.add_argument('--output_complement_normal_fasta',
                        default="tests/complement_normalsim.fasta",
                        help='file path for the output complement 3-5 strand normal (SNV-added) fasta')
    args = vars(parser.parse_args())
    main(args)
