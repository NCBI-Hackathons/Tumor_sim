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
        return

    def insertion(self, genome):
        return genome

    def pick_chromosomes(self, genome, number=1, replace=True):
        relative_lengths = np.array([len(genome[x]) for x in genome])
        probabilities = relative_lengths / float(relative_lengths.sum())
        chroms = np.random.choice(genome.keys(), number, replace=replace,p=probabilities.tolist())
        return chroms

    def get_location_on_sequence(self, seq, distribution='uniform'):
        if distribution == 'uniform':
            # Don't select a location with a N
            while True:
                location = np.random.randint(len(seq))
                if seq[location] != 'N':
                    return location
        else:
            return NotImplementedError('Only Uniform is implemented!')

    def orchestrate_deletion(self, genome, distribution='uniform'):
        chrom = self.pick_chromosomes(genome)[0]
        start = self.get_location_on_sequence(genome[chrom])
        end = start + self.get_event_length()
        if end > len(genome[chrom]):
            end = len(genome[chrom])
        genome[chrom] = self.creator.create_deletion(genome[chrom], start, end)
        return genome

    def orchestrate_translocation(self, genome, distribution='uniform'):
        if len(genome) == 1:
            print('No translocations allowed: genome too small')
            return
        (chrom_source, chrom_target) = self.pick_chromosomes(genome, number = 2, replace = False)
        start_source = self.get_location_on_sequence(genome[chrom_source])
        start_target = self.get_location_on_sequence(genome[chrom_target])
        source_event_length = self.get_event_length(p=0.001)
        target_event_length = self.get_event_length(p=0.001)
        # Make sure event 
        if start_source + source_event_length > len(genome[chrom_source]):
            source_event_length = len(genome[chrom_source]) - start_source
        if start_target + target_event_length > len(genome[chrom_target]):
            target_event_length = len(genome[chrom_target]) - start_target
        return NotImplementedError()  
        # mutated_genome = creator
        print('in orchestrate_translocation')

    # Models exponential decay, discretely, within a 1-10 range. 
    # Expected value of event is 1/p
    def get_event_length(self, p=0.6):
        number = 1
        z = np.random.geometric(p, size=number)
        return z[0]

    def orchestrate_duplication(self, genome, distribution='uniform'):
        print('in orchestrate_duplication')
        return genome

    def orchestrate_inversion(self, genome, distribution='uniform'):
        print ('in orchestrate_inversion')
        return genome

    def orchestrate_insertion(self, genome, distribution='uniform'):
        print ('in orchestrate_insertion')
        return genome

    def generate_structural_variations(self, genome, number):
        variations = np.random.choice(self.structural_variations_probabilities.keys(),
         number, self.structural_variations_probabilities.values())
        mutated_genome = genome
        for variation in variations:
            mutated_genome = self.structural_variations[variation](mutated_genome)


