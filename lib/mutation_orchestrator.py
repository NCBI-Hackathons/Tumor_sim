from Bio import SeqIO
import numpy as np
from mutation_creator import Mutation_Creator

class Mutation_Orchestrator:
    def __init__(self):
        self.creator = Mutation_Creator()
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

    def snv_fast(self, genome, number):
        

    def insertion(self, genome):
        return genome

    def pick_chromosomes(self, genome, number=1, replace=True):
        relative_lengths = np.array([len(genome[x]) for x in genome])
        probabilities = relative_lengths / float(relative_lengths.sum())
        chroms = np.random.choice(genome.keys(), number, probabilities.tolist(), replace)
        return chroms

    def get_location_on_sequence(seq, distribution='uniform'):
        if distribution == 'uniform':
            np.randint(len(seq))

    def orchestrate_deletion(genome, distribution='uniform'):
        chrom = self.pick_chromosome(genome)
        start = self.get_location_on_sequence(genome[chrom])
        # To be completed

    def orchestrate_translocation(genome, distribution='uniform'):
        if len(genome) == 1:
            print('No translocations allowed: genome too small')
            return
        (chrom_source, chrom_target) = self.pick_chromosomes(genome, number = 2, replace = False)
        start_source = self.get_location_on_sequence(genome[chrom_source])
        start_target = self.get_location_on_sequence(genome[chrom_target])
        mutated_genome = creator
        print('in orchestrate_translocation')

    # Models exponential decay, discretely, within a 1-10 range. 
    def get_deletion_length(p=0.6):
        number = 1
        z = np.random.geometric(p, size=number)
        return z[0]

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


