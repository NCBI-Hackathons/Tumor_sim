from Bio import SeqIO
import numpy as np
import pandas as pd
from mutation_creator import Mutation_Creator
import logging

class Mutation_Orchestrator:
    def __init__(self):
        self.tracker = Mutation_Tracker()
        self.creator = Mutation_Creator()
        self.germline_variations = {
        'germline_deletion': self.orchestrate_germline_deletion,
        'germline_insertion' : self.orchestrate_germline_insertion
        }

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
        ### assume for normal bases 
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
            raise NotImplementedError("Only Uniform is implemented!")

    # Models exponential decay, discretely, within a 1-10 range.
    # Expected value of event is 1/p

    def get_event_length(self, p=0.6, number = 1):
        z = np.random.geometric(p, size=number)
        return z[0]

    # Default to being a big deletion, but p=0.6 makes it a small deletion


    def orchestrate_germline_deletion(self, genome, distribution='uniform', p=0.001):
        chrom = self.pick_chromosomes(genome)[0]
        start = self.get_location_on_sequence(genome[chrom])
        end = start + self.get_event_length()
        self.tracker.create_germline_deletion(chrom, start, end)
        #genome[chrom] = self.creator.create_deletion(genome[chrom], start, end)
        logging.info('Orchestrated germline deletion from {} to {} in chrom {}'.format(start, end, chrom))

    def orchestrate_deletion(self, genome, distribution='uniform'):
        chrom = self.pick_chromosomes(genome)[0]
        start = self.get_location_on_sequence(genome[chrom])
        prob = np.random.uniform(0.001, 0.0000001, 1)   ## draw prob from uniform, 0.001 to 1e-7
        end = start + self.get_event_length(p=prob[0])
        self.tracker.create_deletion(chrom, start, end)
        #genome[chrom] = self.creator.create_deletion(genome[chrom], start, end)
        logging.info('Orchestrated deletion from {} to {} in chrom {}'.format(start, end, chrom))



    def orchestrate_translocation(self, genome, distribution='uniform'):
        if len(genome) == 1:
            print('No translocations allowed: genome too small')
            return
        (chrom_source, chrom_target) = self.pick_chromosomes(genome, number = 2, replace = False)
        start_source = self.get_location_on_sequence(genome[chrom_source])
        start_target = self.get_location_on_sequence(genome[chrom_target])
        prob = np.random.uniform(0.001, 0.0000001, 1)   ## draw prob from uniform, 0.001 to 1e-7
        source_event_length = self.get_event_length(p=prob[0])
        target_event_length = self.get_event_length(p=prob[0])
        end_source = start_source+source_event_length
        end_target = start_target+target_event_length
        new_seq_source = genome[chrom_target][start_target:end_target]
        new_seq_target = genome[chrom_source][start_source:end_source]
        self.tracker.create_translocation(chrom_source, chrom_target, start_source,
                        start_target, end_source, end_target, new_seq_source, new_seq_target)
        logging.info('Orchestrated translocation at position {} of length {} on chrom {} to position {} of length {} on chrom {}'.format(
            start_source, source_event_length, chrom_source, start_target, target_event_length, chrom_target))


    # Duplication currently only goes one direction (forward)
    # Creates a variable amount of duplications (num_duplications, drawn from geometric dist)
    def orchestrate_duplication(self, genome, distribution='uniform'):
        chrom = self.pick_chromosomes(genome, number = 1)[0]
        start = self.get_location_on_sequence(genome[chrom])
        prob = np.random.uniform(0.001, 0.0000001, 1)   ## draw prob from uniform, 0.001 to 1e-7
        end = start + self.get_event_length(p=prob[0])
        duplication_prob = np.random.uniform(0.05, 0.7, 1)  ## with np.random.geometric(p, 1), these values are 14 to 2
        num_duplications = self.get_event_length(p=duplication_prob[0]) # exponential ranging from 1 to 10
        new_seq = str(genome[chrom][start:end]) * num_duplications
        self.tracker.create_insertion(chrom, start, new_seq,
             name='duplication (times {})'.format(num_duplications))
        logging.info('Orchestrated duplication at position {} to {} on chrom {}'.format(start, end, chrom))

    def orchestrate_inversion(self, genome, distribution='uniform'):
        chrom = self.pick_chromosomes(genome, number = 1)[0]
        prob = np.random.uniform(0.001, 0.0000001, 1)   ## draw prob from uniform, 0.001 to 1e-7
        start = self.get_location_on_sequence(genome[chrom])
        end = start + self.get_event_length(p=prob[0])
        self.tracker.create_inversion(chrom, start, end)
        logging.info('Orchestrated inversion at position {} to {} on chrom {}'.format(start, end, chrom))

    def orchestrate_insertion(self, genome, distribution='uniform'):
        chrom = self.pick_chromosomes(genome, number = 1)[0]
        start = self.get_location_on_sequence(genome[chrom])
        prob = np.random.uniform(0.001, 0.0000001, 1)
        event_length = self.get_event_length(p=prob[0])
        new_seq_start = self.get_location_on_sequence(genome[chrom])
        new_seq_end = new_seq_start + event_length
        new_seq = genome[chrom][new_seq_start:new_seq_end]
        self.tracker.create_insertion(chrom, start, new_seq)
        logging.info('Orchestrated insertion at position {} on chrom {} adding bases from position {} to {}'.format(start,
         chrom, new_seq_start, new_seq_end))


    def orchestrate_germline_insertion(self, genome, distribution='uniform', p=0.001):
        chrom = self.pick_chromosomes(genome, number = 1)[0]
        start = self.get_location_on_sequence(genome[chrom])
        event_length = self.get_event_length(p)
        new_seq_start = self.get_location_on_sequence(genome[chrom])
        new_seq_end = new_seq_start + event_length
        new_seq = genome[chrom][new_seq_start:new_seq_end]
        self.tracker.create_germline_insertion(chrom, start, new_seq)
        logging.info('Orchestrated insertion at position {} on chrom {} adding bases from position {} to {}'.format(start, chrom, new_seq_start, new_seq_end))



    def generate_structural_variations(self, genome, number):
        variations = np.random.choice(list(self.structural_variations_probabilities.keys()),
                number, self.structural_variations_probabilities.values())
        for variation in variations:
            self.structural_variations[variation](genome)

    # Create small insertions and small deletions, germline
    def generate_germline_indels(self, genome, number):
        variations = np.random.choice(list(['germline_insertion', 'germline_deletion']), number)
        for variation in variations:
            self.germline_variations[variation](genome, p=0.6)
    
    
    # Actually collapses the list of changes
    def generate_fasta(self, genome):
        return self.tracker.collapse_list(genome)

    def get_pandas_dataframe(self):
        return self.bed_correct(self.tracker.log_data_frame.copy())

    # We need to store insertions as same start and end. Let's correct that when outputting bed files
    def bed_correct(self, df):
        same_vals = df[df['start'] == df['end']]
        df.ix[same_vals.index, 'end'] +=1
        return df

class Mutation_Tracker:

    def __init__(self):
        self.creator = Mutation_Creator()
        self.list = []
        self.function_dict = {}
        self.log_data_frame = None

        self.mutation_functions = {
        'deletion': self.creator.create_deletion,
        'germline_deletion': self.creator.create_germline_deletion,
        'translocation': self.creator.create_translocation,
        'inversion' : self.creator.create_inversion,
        'insertion' : self.creator.create_insertion,
        'germline_insertion': self.creator.create_germline_insertion
        }


    def create_insertion(self, chrom, start, new_seq, func='insertion', name='insertion'):
        func_params = [chrom, start, new_seq]
        uid = len(self.list)
        self.list.append([chrom, start, start, name, str(new_seq), uid])
        self.function_dict[uid] = {'func':self.mutation_functions[func], 'params':func_params}
    
    def create_germline_insertion(self, chrom, start, end, func='germline_insertion', name='germline_insertion'):
        uid = len(self.list)
        func_params = [chrom, start, end]
        self.list.append([chrom, start, end, name, '-', uid])
        self.function_dict[uid] = {'func':self.mutation_functions[func], 'params':func_params}

    def create_deletion(self, chrom, start, end, func='deletion', name='deletion'):
        uid = len(self.list)
        func_params = [chrom, start, end]
        self.list.append([chrom, start, end, name, '-', uid])
        self.function_dict[uid] = {'func':self.mutation_functions[func], 'params':func_params}
  
    def create_germline_deletion(self, chrom, start, end, func='germline_deletion', name='germline_deletion'):
        uid = len(self.list)
        func_params = [chrom, start, end]
        self.list.append([chrom, start, end, name, '-', uid])
        self.function_dict[uid] = {'func':self.mutation_functions[func], 'params':func_params}


    def create_inversion(self, chrom, start, end, func='inversion', name='inversion'):
        uid = len(self.list)
        func_params = [chrom, start, end]
        self.list.append([chrom, start, end, name, '-', uid])
        self.function_dict[uid] = {'func':self.mutation_functions[func], 'params':func_params}
        

    def create_translocation(self, chrom_source, chrom_target, start_source, start_target,
                end_source, end_target, new_seq_source, new_seq_target):
        uid = len(self.list)
        # Deletions
        self.create_deletion(chrom_source, start_source, end_source, name='translocation(del)')
        self.create_deletion(chrom_target, start_target, end_target, name='translocation(del)')
        # names for insertions translocation(chrom:start-end)
        name_insertion_source = ''.join(['translocation(', chrom_target, ':', str(start_target), '-', str(end_target), ')'])
        name_insertion_target = ''.join(['translocation(', chrom_source, ':', str(start_source), '-', str(end_source), ')'])
        # insertions
        self.create_insertion(chrom_source, start_source, new_seq_source, name=name_insertion_source)
        self.create_insertion(chrom_target, start_target, new_seq_target, name=name_insertion_target)


    def collapse_list(self, genome):
        self.log_data_frame = pd.DataFrame(self.list)
        self.log_data_frame.columns = ['chrom', 'start', 'end', 'name', 'alt', 'uid']
        self.log_data_frame = self.log_data_frame.sort_values(['chrom', 'start', 'end'], ascending = [False, False, False])
        previous_starts = {}
        for chrom in genome:
            previous_starts[chrom] = float('Inf')

        mutable_genome = genome
        for idx, row in self.log_data_frame.iterrows():
            # Logic for avoiding overlaps
            uid = row.uid
            chrom = row.chrom
            if row.end <= previous_starts[chrom]:
                # Have to modify the parameters because the mutable seq needs to get passed in
                mutable_params = self.function_dict[uid]['params']
                mutable_params[0] = genome[chrom]
                mutable_genome[chrom] = self.function_dict[uid]['func'](*mutable_params)
            else:
                # Drop row from DataFrame
                self.log_data_frame = self.log_data_frame.drop(idx, axis=0)
            previous_starts[chrom] = row.start
        return(mutable_genome)
