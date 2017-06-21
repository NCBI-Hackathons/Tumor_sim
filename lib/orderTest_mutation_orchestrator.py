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

### assume chr1----implement a dictionary for each chromosome, with value 0
length_normal_deletions_implemented = 0

    def orchestrate_normal_deletion(self, seq, distribution="uniform"):        
        """
        Keep track of function calls, and iterate by chromosome.
        For each function call, execute `len(seq)-# of deletions`
        When calling this function, do not randomly choose chromosomes! 
        Tracking the # of function calls would make no sense!
        IDEA: posisbly use `from collections import Counter` to keep track of Chromosomes
        """
         ### implement this by chromosome
         ### now assume only chr1 seq
        chrom = "chr1"
        if distribution == "uniform":
            ## avoid telomere sequences/other masks
            while True:
                start = np.random.randint(int(len(seq) - length_normal_deletions_implemented))
                if seq[start] != "N":
                    end = start + self.get_event_length()
                    length_normal_deletions_implemented += end   ## length of event
                    ## now, do this ~10K times for chromosome, and sort these events.
                    ### after sorting, when you implement the mutation, record it in a VCF
                    ###genome[chrom] = self.creator.create_deletion(genome[chrom], start, end)
                    logging.info('Orchestrated deletion from {} to {} in chrom {}'.format(start, end, chrom))
                    return genome
        else:
            raise NotImplementedError("We've only implemented a uniform distribution of Dels!")

### implement a dictinoary: keys of chromsosomes, values of [0]

length_normal_insertions_implemented = 0

    def orchestrate_normal_insertion(self, seq distribution="uniform"):        
        """
        Keep track of function calls, and iterate by chromosome.
        For each function call, execute `len(seq) + # of insertions`
        When calling this function, do not randomly choose chromosomes! 
        Tracking the # of function calls would make no sense!
        IDEA: posisbly use `from collections import Counter` to keep track of Chromosomes
        """
         ### implement this by chromosome
         ### now assume only chr1 seq
        chrom = "chr1"
        if distribution == "uniform":
            ## avoid telomere sequences/other masks
            while True:
                start = np.random.randint(int(len(seq) + length_normal_insertions_implemented)
                if seq[start] != "N":
                    end = start + self.get_event_length()
                    length_normal_insertions_implemented += end   ### length of event
                    genome[chrom] = self.creator.create_deletion(genome[chrom], start, end)
                    logging.info('Orchestrated insertion from {} to {} in chrom {}'.format(start, end, chrom))
                    return genome
        else:
            raise NotImplementedError("We've only implemented a uniform distribution of Insertions!")



### Now, fix the dictionary problem:
### It doesn't matter what the dicitonary keys are for the input FASTA
### Begin by creating two dictionaries:
### -- normal_insert_length_dict 
### -- tumor_insert_length_dict
### Create the keys from FASTA dictionary, 'original_genome' (?), and give these dicitionaries these keys
### Then initialize each value = 0 
### 
### original dictionary from FASTA is 'genome'; CHECK


## dictionary with keys from FASTA dictionary, values = 0

insertion_length_tracker = { k: 0 for k in genome.iteritems()} 

deletion_length_tracker = {k: 0 for k in genome.iteritems()}

### do SNPS first, then insertions, then deletions

### do this by chromosome
### create a list of tuples s.t. (start, end) for each
### then sort, and apply at once

    def orchestrate_normal_insertion(self, seq distribution="uniform"):        
        """
        Keep track of function calls, and iterate by chromosome.
        For each function call, execute `len(seq) + # of insertions`
        When calling this function, do not randomly choose chromosomes! 
        Tracking the # of function calls would make no sense!
        IDEA: posisbly use `from collections import Counter` to keep track of Chromosomes
        """
        ## chrom = self.pick_chromosomes(genome)[0]   ## get random chromosome
        if distribution == "uniform":
            ## avoid telomere sequences/other masks
            while True:
                start = np.random.randint(int(len(seq) + insertion_length_tracker[chrom])
                if seq[start] != "N":
                    length = self.get_event_length()
                    end = start + length
                    insertion_length_tracker[chrom] += length   ### add length of event
                    ### create list of tuples, s.t. tuple = (start, end)
                    ### sort list of tuples 
                    ### implement rearrangmenets by largest coordinate to smallest
                    genome[chrom] = self.creator.create_deletion(genome[chrom], start, end)
                    ### write this into a VCF
                    logging.info('Orchestrated insertion from {} to {} in chrom {}'.format(start, end, chrom))
                    return genome
        else:
            raise NotImplementedError("We've only implemented a uniform distribution of Insertions!")















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

    def generate_structural_variations(self, genome, number):
        variations = np.random.choice(list(self.structural_variations_probabilities.keys()),
         number, self.structural_variations_probabilities.values())
        mutated_genome = genome
        for variation in variations:
            mutated_genome = self.structural_variations[variation](mutated_genome)
        return mutated_genome


