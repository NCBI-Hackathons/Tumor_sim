from Bio import SeqIO
from mutation_orchestrator import Mutation_Orchestrator

fasta_file = 'test_fasta.fasta'
iterations = 100

def main():
    # read genome fasta
    genome = {}
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        genome[seq_record.id] = seq_record.seq

    orchestrator = Mutation_Orchestrator()

    # add SNVs
    mutated_genome = genome
    for i in range(iterations):
        mutated_genome = orchestrator.snv(mutated_genome)

    # write fasta
    with open("output.fasta", "w") as output_handle:
        SeqIO.write(mutated_genome, output_handle, "fasta")

if __name__ == "__main__":
    main()