from Bio import SeqIO
import numpy as np

class Mutation_Orchestrator:
    def __init__(self):
        self.structural_variations = {
        'deletion': self.orchestrate_deletion,
        'translocation': self.orchestrate_translocation,
        'duplication' : self.orchestrate_duplication,
        'inversion' : self.orchestrate_inversion,
        'insertion' : self.orchestrate_insertion
            }
    
        self.structural_variations_probabilities = {
        'deletion': 0.2,
        'translocation': 0.2,
        'duplication' : 0.2,
        'inversion' : 0.2,
        'insertion' : 0.2
            }

    def deletion(self, genome):
        return genome

    def snv(self, genome):
        return genome

    def insertion(self, genome):
        return genome

    def pick_chromosome(self, genome):
        import pdb; pdb.set_trace()
        relative_lengths = np.array([len(genome[x]) for x in genome])
        probabilities = relative_lengths / float(relative_lengths.sum())
        return np.random.choice(genome.keys(), 1, probabilities.tolist())

    def get_location_on_sequence(seq, distribution='uniform'):
        if distribution == 'uniform':
            np.randint(len(seq))

    def orchestrate_deletion(genome, distribution='uniform'):
        chrom = pick_chromosome(genome)
        get_location_on_sequence()

    def orchestrate_translocation(genome, distribution='uniform'):
        print('in orchestrate_translocation')

    def orchestrate_duplication(genome,  distribution='uniform'):
        print('in orchestrate_duplication')

    def orchestrate_inversion(genome,  distribution='uniform'):
        print ('in orchestrate_inversion')

    def orchestrate_insertion(genome,  distribution='uniform'):
        print ('in orchestrate_insertion')

    def generate_structural_variations(self, genome, number):
        variations = np.random.choice(self.structural_variations_probabilities.keys(),
         number, self.structural_variations_probabilities.values())
        mutated_genome = genome
        for variation in variations:
            mutated_genome = self.structural_variations[variation](mutated_genome)


