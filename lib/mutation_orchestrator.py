from Bio import SeqIO
import random
from mutation_creator import Mutation_Creator


class Mutation_Orchestrator:
    def __init__(self, prob_table_file,final_seqfile,logfile):
        self.prob_table_file = prob_table_file
        self.final_seqfile = final_seqfile
        self.logfile = logfile
        self.mc = Mutation_Creator()

    def deletion(self, genome):
        return genome

    def add_snvs(self, sequence, chrom, in_memory=True, dampen = 0.01):
        #bug: no mutation at ends
        seq = str(sequence).upper()
        final_seq = open(self.final_seqfile, "a")
        lfile = open(self.logfile, "a")
        lfile.write("chrom\tnucleotide\toriginal\tmutated")
        prob_file = open(self.prob_table_file, "r")
        prob_table = []
        key = []
        b = prob_file.readline()
        # skip the first line
        b = prob_file.readline()
        ic = 0
        # Create probability table, reading prob_file
        while b != "" and ic < 100:
            e = b.replace("\n","")
            f = e.split("\t")
            if len(f)<2:
                ic += 1
            else:
                if key.count(f[0]) == 0:
                    key.append(f[0])
                    prob_table.append([f[0], f[1],float(f[2])])
                else:
                    check = key.index(f[0])
                    prob_table[check] += [f[1],float(f[2])]
            b = prob_file.readline()
        c = 0
        ##CHANGE TO NUMPY
        random.seed()
        while c < len(seq) - 2:
            ##CHANGE TO NUMPY
            randomvalue = random.random() / dampen
            triplet = seq[c:c + 3]
            # ERROR: need to save first and last 
            original_base = triplet[1]
            if c % 1000000 == 0:
                print('In loop {} of add_snvs'.format(c))
            if 'N' not in triplet:
                check = key.index(triplet)
                # Only look at the old sequence data set 
                if randomvalue > prob_table[check][6] + prob_table[check][4] + prob_table[check][2]:
                    new_base = original_base
                elif randomvalue > prob_table[check][4] + prob_table[check][2]:
                    new_base = prob_table[check][5][1]
                elif randomvalue > prob_table[check][2]:
                    new_base = prob_table[check][3][1]
                else:
                    new_base = prob_table[check][1][1]
            else:
                new_base = original_base
            if new_base != original_base:
                lfile.write("\n" + chrom + "\t" + str(c + 2) + "\t"  + original_base + "\t" + new_base)
                if in_memory:
                    self.mc.create_snv(sequence, c, new_base)
            if not in_memory:
                final_seq.write(new_base)
            c += 1
        return sequence
    
    def insertion(self, genome):
        return genome

    def add_snvs_across_genome(self, genome):
        for chr in genome:
            genome[chr] = self.add_snvs(genome[chr], chr)
        return genome




