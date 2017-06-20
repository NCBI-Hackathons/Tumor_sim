class Mutation_Tracker():

    # Keep in mind the tracker lengths have had the Ns removed, and will need to be added back to it after
    def create_structure(genome):
        self.structure = {}
        for chrom in genome:
            self.structure[chrom] = [{chrom:(0, len(genome[chrom]))}]

    def get_interval(l):
        return l.values()[0]

    def track_deletion(self, chrom, start, end):
        # chunks_to_modify = self.find_chunks(self.structure[chrom], start, end)
        feature_list = self.structure[chrom]
        current = 0
        i = 0
        while True:
            interval = self.get_interval(feature_list[i])
            diff = np.abs(interval[0] - interval[1]) 
            if start < (current + diff):
                start_i = i
                start_offset = start - current
            if end < (current + diff):
                end_i = i
                end_offset = (current + diff) - end
                break
            i = i + 1
            current = current + diff
        chunks_to_modify = feature_list[start_i:end_i]
        interval = self.get_interval(chunks_to_modify[0])
        new_first_entry = {chunks_to_modify[0].keys()[0] : (interval[0], interval[0] + start_offset)}
        interval = self.get_interval(chunks_to_modify[-1])
        new_second_entry = {chunks_to_modify[-1].keys()[0]: (interval[1] - end_offset, interval[1])}
        # feature_list[start_i:end_i] = [new_first_entry, new_second_entry]
        return [new_first_entry, new_second_entry]

    def track_snv(self, mutable_seq, start, new_base):
        mutable_seq[start] = new_base
        return mutable_seq

    def find_chunks(feature_list, start, end)
        current = 0
        i = 0
        while True:
            diff = np.abs(feature_list[i][0] - feature_list[i][1]) 
            if start < current + diff:
                start_i = i
            if end < current + diff:
                end_i = i
            i = i + 1
            current = current + diff
        return feature_list[start_i:end_i]

    def track_insertion(self, mutable_seq, start, new_seq):
        return mutable_seq[:start] + new_seq +  mutable_seq[start:]

    def track_inversion(self, mutable_seq, start, end):
        inv_seq = mutable_seq[start:end]
        inv_seq.reverse()
        mutable_seq[start:end] = inv_seq
        return mutable_seq

    def track_translocation(self, seq1, seq2, start1, start2, length1, length2):
        new_seq1 = self.create_deletion(seq1, start1,start1+length1)
        new_seq1 = self.create_insertion(new_seq1, start1, seq2[start2:start2+length2])
        new_seq2 = self.create_deletion(seq2, start2,start2+length2)
        new_seq2 = self.create_insertion(new_seq2, start2, seq1[start1:start1+length1])
        return (new_seq1, new_seq2)