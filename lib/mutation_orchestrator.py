from Bio import SeqIO
import numpy as np
from mutation_creator import Mutation_Creator
from mutation_tracker import Mutation_Tracker
import logging
from probabilities_config import structural_variations_probabilities, snv_probabilities

class Mutation_Orchestrator:
    """ Mutation_Orchestrator is a class that operates on a genome to make a mutation.
        It chooses the type, length, and start of the mutation probabilistically, using 
        the probabilities defined in probabilities_config. It stores a record of the mutations
        in the Mutation_Tracker class, and modifies the genome string using Mutation_Creator""" 
    def __init__(self):
        self.tracker = Mutation_Tracker()
        self.creator = Mutation_Creator()
        self.structural_variations = {
        'deletion': self.orchestrate_deletion,
        'translocation': self.orchestrate_translocation,
        'duplication' : self.orchestrate_duplication,
        'inversion' : self.orchestrate_inversion,
        'insertion' : self.orchestrate_insertion
        }

        self.logger = logging.basicConfig(filename='example.log',level=logging.DEBUG)

    def snv_fast(self, genome, number):
        chroms = self.pick_chromosomes(genome, number)
        ### assume for normal bases 
        new_bases = np.random.choice(list(snv_probabilities.keys()), number, snv_probabilities.values())
        for i in range(number):
            start = self.get_location_on_sequence(genome[chroms[i]])
            genome[chroms[i]] = self.creator.create_snv(genome[chroms[i]], start, new_bases[i])
            logging.info('Added base {} at position {} in chrom {}'.format(new_bases[i], str(start), chroms[i]))
        return genome

    def pick_chromosomes(self, genome, number=1, replace=True):
        relative_lengths = np.array([len(genome[x]) for x in genome])
        probabilities = relative_lengths / float(relative_lengths.sum())
        chroms = np.random.choice(list(genome.keys()), number, replace=replace ,p=probabilities.tolist())
        return chroms

    def get_location_on_sequence(self, seq, distribution='uniform'):
        if distribution == 'uniform':
            # Don't select a location with a N
            while True:
                location = np.random.randint(len(seq))
                if seq[location] != 'N':
                    return location
        else:
            raise NotImplementedError("Only Uniform is implemented!")

    # Default to being a big deletion, but p=0.6 makes it a small deletion
    def orchestrate_deletion(self, genome, distribution='uniform', p=0.001):
        chrom = self.pick_chromosomes(genome)[0]
        start = self.get_location_on_sequence(genome[chrom])
        end = self.get_end_of_event(start, genome[chrom], p)
        self.tracker.create_deletion(chrom, start, end)
        logging.info('Orchestrated deletion from {} to {} in chrom {}'.format(start, end, chrom))

    def orchestrate_translocation(self, genome, distribution='uniform'):
        if len(genome) == 1:
            return
        (chrom_source, chrom_target) = self.pick_chromosomes(genome, number = 2, replace = False)
        start_source = self.get_location_on_sequence(genome[chrom_source])
        start_target = self.get_location_on_sequence(genome[chrom_target])
        end_source = self.get_end_of_event(start_source, genome[chrom_source], p=0.001)
        end_target = self.get_end_of_event(start_target, genome[chrom_target], p=0.001)
        new_seq_source = genome[chrom_target][start_target:end_target]
        new_seq_target = genome[chrom_source][start_source:end_source]
        self.tracker.create_translocation(chrom_source, chrom_target, start_source,
                        start_target, end_source, end_target, new_seq_source, new_seq_target)
        logging.info('Orchestrated translocation at position {} to position {} on chrom {} at position {} to position {} on chrom {}'.format(
            start_source, end_source, chrom_source, start_target, end_target, chrom_target))

    # Guarantees the end of the event isn't outside the sequence
    def get_end_of_event(self, start_pos, seq, p):
        length = self.get_event_length(p)
        return np.min([start_pos + length, len(seq)])

    # Models exponential decay, discretely
    # Expected value of event is 1/p
    def get_event_length(self, p=0.6, number = 1):
        z = np.random.geometric(p, size=number)
        return z[0]

    # Duplication currently only goes one direction (forward)
    # Creates a variable amount of duplications (num_duplications, drawn from geometric dist)
    def orchestrate_duplication(self, genome, distribution='uniform'):
        chrom = self.pick_chromosomes(genome, number = 1)[0]
        start = self.get_location_on_sequence(genome[chrom])
        end = self.get_end_of_event(start, genome[chrom], p=0.001)
        num_duplications = self.get_event_length(p=0.6) # exponential ranging from 1 to 10
        new_seq = str(genome[chrom][start:end]) * num_duplications
        self.tracker.create_insertion(chrom, start, new_seq,
             name='duplication (times {})'.format(num_duplications))
        logging.info('Orchestrated duplication at po`tion {} to {} on chrom {}'.format(start, end, chrom))

    def orchestrate_inversion(self, genome, distribution='uniform'):
        chrom = self.pick_chromosomes(genome, number = 1)[0]
        start = self.get_location_on_sequence(genome[chrom])
        end = self.get_end_of_event(start, genome[chrom], p=0.001)
        self.tracker.create_inversion(chrom, start, end)
        logging.info('Orchestrated inversion at position {} to {} on chrom {}'.format(start, end, chrom))

    def orchestrate_insertion(self, genome, distribution='uniform', p=0.001):
        chrom = self.pick_chromosomes(genome, number = 1)[0]
        start = self.get_location_on_sequence(genome[chrom])
        new_seq_start = self.get_location_on_sequence(genome[chrom])
        new_seq_end = self.get_end_of_event(new_seq_start, genome[chrom], p)
        new_seq = genome[chrom][new_seq_start:new_seq_end]
        self.tracker.create_insertion(chrom, start, new_seq)
        logging.info('Orchestrated insertion at position {} on chrom {} adding bases from position {} to {}'.format(start,
         chrom, new_seq_start, new_seq_end))

    def generate_structural_variations(self, genome, number):
        variations = np.random.choice(list(structural_variations_probabilities.keys()),
                number, structural_variations_probabilities.values())
        for variation in variations:
            self.structural_variations[variation](genome)

    # Create small insertions and small deletions
    def generate_indels(self, genome, number):
        variations = np.random.choice(list(['insertion', 'deletion']), number)
        for variation in variations:
            self.structural_variations[variation](genome, p=0.6)

    # Actually collapses the list of changes    
    def generate_fasta_and_bed(self, genome):
        (genome, log_data_frame) =  self.tracker.collapse_list(genome)
        return (genome, self.bed_correct(log_data_frame))

    # We need to store insertions as same start and end. Let's correct that when outputting bed files
    def bed_correct(self, df):
        same_vals = df[df['start'] == df['end']]
        df.ix[same_vals.index, 'end'] +=1
        return df
