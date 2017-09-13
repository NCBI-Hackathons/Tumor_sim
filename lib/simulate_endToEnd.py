from Bio import SeqIO
from mutation_orchestrator import Mutation_Orchestrator
import copy
import pandas as pd
import re
import argparse
import copy

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
    original_len_seq = len(sequence)
    while sequence[0] == 'N':
        sequence.pop(0)
    offset = original_len_seq - len(sequence)
    while sequence[-1] == 'N':
        sequence.pop(-1)
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
    new_genome = copy.deepcopy(genome)
    for chrom in new_genome:
        new_genome[chrom].complement()
    return new_genome 

def main(args):
     # read genome fasta
    (mutated_genome, genome_offset) = read_fasta_normal(args['input_fasta'])

    orchestrator = Mutation_Orchestrator()
    # add germilne SNVs & InDels
    mutated_genome = orchestrator.snv_fast(mutated_genome, args['number_snvs'])
    orchestrator.generate_indels(mutated_genome, args['number_indels'])
    (mutated_genome, snv_and_indel_bed) = orchestrator.generate_fasta_and_bed(mutated_genome)
    ### write out "normalsim" bedpe and fasta
    write_fasta(mutated_genome, args['output_normal_fasta'])
    write_bed(genome_offset, snv_and_indel_bed, args['output_normal_bedfile'])  

    ## output complement 3'-5' strand normal
    normal_complement = create_complementary_genome(mutated_genome)
    write_fasta(normal_complement, args['output_complement_normal_fasta'])
    del normal_complement

    # add structural varations
    orchestrator.generate_structural_variations(mutated_genome, args['number_of_tumorSVs'])
    (mutated_genome, tumor_bed) = orchestrator.generate_fasta_and_bed(mutated_genome)
    write_fasta(mutated_genome, args['output_tumor_fasta'])
    write_bed(genome_offset, tumor_bed, args['output_tumor_bedfile'])  ### write out "tumorsim" bedpe

    ## output complement 3'-5' strand tumor
    tumor_complement = create_complementary_genome(mutated_genome)
    write_fasta(tumor_complement, args['output_complement_tumor_fasta'])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Simulate cancer genomic structural variations')
    parser.add_argument('--input_fasta', 
                        default="../data/subsampled_hg38.fa",
                        help='file path for the input (default genome) fasta')
    parser.add_argument('--output_tumor_fasta',
                        default="outputs/tumorsim.fasta",
                        help='file path for the output tumor (cancer genome) fasta')
    parser.add_argument('--output_normal_fasta',
                        default = "outputs/normalsim.fasta",
                        help='file path for the output normal (SNV-added) fasta')
    parser.add_argument('--output_complement_tumor_fasta',
                        default="outputs/complement_tumorsim.fasta",
                        help='file path for the output complement 3-5 strand tumor (cancer genome) fasta')
    parser.add_argument('--output_complement_normal_fasta',
                        default="outputs/complement_normalsim.fasta",
                        help='file path for the output complement 3-5 strand normal (SNV-added) fasta')
    parser.add_argument('--number_snvs',
                        default = 3000,
                        help="number of single nucleotide variants to add to the normal genome")
    parser.add_argument('--number_indels',
                        default = 4150,
                        help="number of small insertions and deletions to add to the normal genome")
    parser.add_argument('--number_of_tumorSVs',
                        default = 10000,
                        help="number of structural variations to add to the tumor genome")
    parser.add_argument('--output_normal_bedfile',
                        default = "outputs/normal.bed",
                        help='file path for the output normal (SNV-added) bedfile')
    parser.add_argument('--output_tumor_bedfile',
                        default="outputs/tumorsim.bed",
                        help='file path for the output tumor (cancer genome) bedfile')
    args = vars(parser.parse_args())
    main(args)
