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

    def add_snvs(self, sequence, in_memory=True):
        #bug: no mutation at ends
        seq = str(sequence).upper()
        final_seq = open(self.final_seqfile, "w")
        lfile = open(self.logfile, "w")
        lfile.write("nucleotide\toriginal\tmutated")
        prob_file = open(self.prob_table_file, "r")
        prob_table = []
        key = []
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
                    prob_table.append(f)
                else:
                    check = key.index(f[0])
                    prob_table[check] += [f[1],float(f[2])]
            b = prob_file.readline()
        c = 0
        ##CHANGE TO NUMPY
        random.seed()
        while c < len(seq) - 2:
            ##CHANGE TO NUMPY
            randomvalue = random.random()
            triplet = seq[c:c + 3]
            # ERROR: need to save first and last 
            original_base = prob_table[check][0][1]
            if c % 1000000 == 0:
                print('In loop {} of add_snvs'.format(c))
            if 'N' not in triplet:
                # Only look at the old sequence data set 
                if randomvalue > float(prob_table[check][6]) + float(prob_table[check][4]) + float(prob_table[check][2]):
                    new_base = original_base
                elif randomvalue > float(prob_table[check][4]) + float(prob_table[check][2]):
                    new_base = prob_table[check][5][1]
                elif randomvalue > float(prob_table[check][2]):
                    new_base = prob_table[check][3][1]
                else:
                    new_base = prob_table[check][1][1]
            else:
                new_base = original_base
            if in_memory:
                if new_base != original_base:
                    self.mc.create_snv(sequence, c, new_base)
            else:
                final_seq.write(new_base)
                if new_base != original_base:
                    lfile.write("\n" + str(c + 2) + "\t" + prob_table[check][0][1] + "\t" + new_base)
            c += 1
    
        return sequence
    
    def insertion(self, genome):
        return genome

    def add_snvs_across_genome(self, genome):
        self.add_snvs(genome['chr1'])



