from Bio import SeqIO
import numpy as np
import pandas as pd
import random

class Overlap_Shuffler:
	def shuffler(self, dframe):
		return dframe.reindex(np.random.permutation(dframe.index))

	def combined_unions(self, intervals):
        ### finds overlapping intervals and takes the union, merges into one union interval
        sorted_by_lowest_interval = sorted(intervals, key=lambda tup: tup[0])
        merged_intervals = []
        for higher in sorted_by_lowest_interval:
            if not merged_intervals:
                merged_intervals.append(higher)
            else:
                lower = merged_intervals[-1]
                # check for union between 'higher' and 'lower':
                # after sorting above, by definition lower[0] <= higher[0]
                if higher[0] <= lower[1]:
                    upper_bound = max(lower[1], higher[1])
                    merged_intervals[-1] = (lower[0], upper_bound)  # replace by combined interval
                else:
                    merged_intervals.append(higher)
        return merged_intervals

    def remix_overlaps(dframe):
        ### create list of tuples
        chrom_list = []   ## append each shuffled dataframe by chromosome
        for name, group in dframe.groupby("chrom"):  ## split dataframe by group
            chrom = pd.DataFrame(group)  ### each chrom has its own dataframe
            list_of_merged_intervals = self.combined_unions(zip(chrom["start"], chrom["end"]))
            chrom["interval_ID"] = chrom["start"].apply(lambda x: next(i for i, m in enumerate(list_of_merged_intervals) if m[0] <= x <= m[1]))
            shuffled_group = chrom.groupby("interval_ID").apply(self.shuffler)
            shuffled_group = shuffled_group.drop("interval_ID", axis=1) ### drop column 'interval_ID'
            shuffled_group = shuffled_group.reset_index(drop=True)
            chrom_list.append(shuffled_group)
        total_chroms = pd.concat(chrom_list)  ### probably need to shuffle
        return total_chroms

