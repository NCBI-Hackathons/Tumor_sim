from Bio import SeqIO
import numpy as np
from mutation_creator import Mutation_Creator
import logging

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
        self.logger = logging.basicConfig(filename='example.log',level=logging.DEBUG)

    def snv_fast(self, genome, number):
        chroms = self.pick_chromosomes(genome, number)
        new_bases = np.random.choice(['A', 'C', 'T', 'G'], number, [0.25, 0.25, 0.25, 0.25])
        for i in range(number):
            start = self.get_location_on_sequence(genome[chroms[i]])
            genome[chroms[i]] = self.creator.create_snv(genome[chroms[i]], start, new_bases[i])
            logging.info('Added base {} at position {} in chrom {}'.format(new_bases[i], str(start), chroms[i]))
        return genome

    def pick_chromosomes(self, genome, number=1, replace=True):
        relative_lengths = np.array([len(genome[x]) for x in genome])
        probabilities = relative_lengths / float(relative_lengths.sum())
        chroms = np.random.choice(list(genome.keys()), number, replace=replace,p=probabilities.tolist())
        return chroms

    def get_location_on_sequence(self, seq, distribution='uniform'):
        if distribution == 'uniform':
            # Don't select a location with a N
            while True:
                location = np.random.randint(len(seq))
                if seq[location] != 'N':
                    return location
        else:
            raise NotImplementedError('Only Uniform is implemented!')

    def orchestrate_deletion(self, genome, distribution='uniform'):
        chrom = self.pick_chromosomes(genome)[0]
        start = self.get_location_on_sequence(genome[chrom])
        end = start + self.get_event_length()
        genome[chrom] = self.creator.create_deletion(genome[chrom], start, end)
        logging.info('Orchestrated deletion from {} to {} in chrom {}'.format(start, end, chrom))
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
        (genome[chrom_source], genome[chrom_target]) = self.creator.create_translocation(
            genome[chrom_source],
         genome[chrom_target], start_source, start_target, source_event_length, target_event_length)
        logging.info('Orchestrated translocation at position {} of length {} on chrom {} to position {} of length {} on chrom {}'.format(
            start_source, source_event_length, chrom_source, start_target, target_event_length, chrom_target))
        return genome

    # Models exponential decay, discretely, within a 1-10 range. 
    # Expected value of event is 1/p
    def get_event_length(self, p=0.6):
        number = 1
        z = np.random.geometric(p, size=number)
        return z[0]

    # Duplication currently only goes one direction (forward)
    def orchestrate_duplication(self, genome, distribution='uniform'):
        chrom = self.pick_chromosomes(genome, number = 1)[0]
        start = self.get_location_on_sequence(genome[chrom])
        end = start + self.get_event_length(p=0.001)
        genome[chrom] = self.creator.create_insertion(genome[chrom], start, genome[chrom][start:end])
        logging.info('Orchestrated duplication at position {} to {} on chrom {}'.format(start, end, chrom))
        return genome

    def orchestrate_inversion(self, genome, distribution='uniform'):
        chrom = self.pick_chromosomes(genome, number = 1)[0]
        start = self.get_location_on_sequence(genome[chrom])
        end = start + self.get_event_length(p=0.001)
        genome[chrom] = self.creator.create_inversion(genome[chrom], start, end)
        logging.info('Orchestrated inversion at position {} to {} on chrom {}'.format(start, end, chrom))
        return genome

    def orchestrate_insertion(self, genome, distribution='uniform'):
        print ('in orchestrate_insertion')
        chrom = self.pick_chromosomes(genome, number = 1)[0]
        start = self.get_location_on_sequence(genome[chrom])
        event_length = self.get_event_length(p=0.6)
        new_seq_start = self.get_location_on_sequence(genome[chrom])
        new_seq_end = new_seq_start + event_length
        new_seq = genome[chrom][new_seq_start:new_seq_end]
        genome[chrom] = self.creator.create_insertion(genome[chrom], start, new_seq)
        logging.info('Orchestrated insertion at position {} on chrom {} adding bases from position {} to {}'.format(start,
            chrom, new_seq_start, new_seq_end))
        return genome

    # Keep in mind the tracker lengths have had the Ns removed, and will need to be added back to it after
    def get_tracker_object(genome):
        tracker = {}
        for chrom in genome:
            self.tracker[chrom] = [{chrom:(0, len(genome[chrom]))}]

    def generate_structural_variations(self, genome, number):
        variations = np.random.choice(list(self.structural_variations_probabilities.keys()),
         number, self.structural_variations_probabilities.values())
        self.creator.tracker = self.get_tracker_object(genome)
        mutated_genome = genome
        for variation in variations:
            mutated_genome = self.structural_variations[variation](mutated_genome)
        return mutated_genome


