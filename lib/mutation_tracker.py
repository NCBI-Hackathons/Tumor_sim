from mutation_creator import Mutation_Creator
import pandas as pd

class Mutation_Tracker:
    """ Mutation_Tracker tracks the history of mutations that modify a genome. 
        It shares many function names with Mutation_Creator, calling into it 
        to make changes to the mutable sequences, storing those changes in 
        bed format in self.list, and a dictionary with pointers to the function
        calls in self.function_dict. It uses self.function_dict so that the changes
        can be made in position order, maintaining correct position tracking."""
    def __init__(self):
        self.creator = Mutation_Creator()
        self.list = []
        self.function_dict = {}

        self.mutation_functions = {
        'deletion': self.creator.create_deletion,
        'translocation': self.creator.create_translocation,
        'inversion' : self.creator.create_inversion,
        'insertion' : self.creator.create_insertion
        }


    def create_insertion(self, chrom, start, new_seq, func='insertion', name='insertion'):
        func_params = [chrom, start, new_seq]
        uid = len(self.list)
        self.list.append([chrom, start, start, name, str(new_seq), uid])
        self.function_dict[uid] = {'func':self.mutation_functions[func], 'params':func_params}

    def create_deletion(self, chrom, start, end, func='deletion', name='deletion'):
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
        """ Takes all pending mutations, and orders them along the genome from end to start (descending).
        If the changes are overlapping, the one closer to the genome start whose end extends past the 
        previous mutation (by descending order's) start gets dropped"""
        log_data_frame = pd.DataFrame(self.list)
        log_data_frame.columns = ['chrom', 'start', 'end', 'name', 'alt', 'uid']
        log_data_frame = log_data_frame.sort_values(['chrom', 'start', 'end'], ascending = [False, False, False])
        previous_starts = {}
        for chrom in genome:
            previous_starts[chrom] = float('Inf')

        mutable_genome = genome
        to_drop = []
        for uid in log_data_frame.index:
        # Logic for avoiding overlaps
            chrom = log_data_frame.loc[uid, 'chrom']
            if log_data_frame.loc[uid, 'end'] <= previous_starts[chrom]:
                # Have to modify the parameters because the mutable seq needs to get passed in
                func = self.function_dict.pop(uid)
                func['params'][0] = genome[chrom]
                mutable_genome[chrom] = func['func'](*func['params'])
            else:
                # Drop row from DataFrame
                to_drop.append(uid)
            previous_starts[chrom] = log_data_frame.loc[uid, 'start']
        log_data_frame.drop(to_drop, axis=0, inplace=True)
        # Reset the tracker's list to empty, since it has been collapsed into the log_data_frame
        self.list = []
        return(mutable_genome, log_data_frame)