#!/usr/bin/python

import math
from rpy import *
import operator
     
class Sheep:
    def __init__(self,name):
        self.name = name
        self.breed = None
        self.numreads = None
        self.regions = None # List of all Regions 
        self.snps = []
        self.sequences = None # {Region:sequence}

class SNP:
    #Ref_Seq;Ref_Position;Ref_Base;Mutant_base;HQ_Depth;A;C;G;T
    def __init__(self,ref_seq,ref_position,ref_base):
        self.position = int(ref_position)
        self.region = ref_seq
        self.mutant_base = None
        self.ref_base = ref_base
        self.hq_depth = None
        self.basecount = {"A":0,"C":0,"G":0,"T":0}
        self.sheep = None
        #SNP genotype is stored as [(base,count)], len ==1 --> homozygote
        #len == 2 --> heterozygote
        #None == Genotype unknown!
        self.genotype = ()

    #Version 2 with error rate
    def calc_genotype(self):
        fac = math.factorial
        ln = math.log
        n = sorted(self.basecount.values(), reverse=True)
        all_reads = float(sum(n))
        #Part1
        part1 = n[0]*ln(n[0]/all_reads)
        #Part2
        if all_reads == n[0]:
            part2 = 0
        else:
            part2 = (all_reads-n[0])*ln((all_reads-n[0])/(3*all_reads))
        #Part3
        ln_part3 = (n[0]+n[1])/(2*all_reads)
        part3 = (n[0]+n[1])*ln(ln_part3)
        #Part4
        ln_part4 = (n[2]+n[3])/(2*all_reads)
        if ln_part4 == 0.0:
            part4 = 0
        else:
            part4 = (n[2]+n[3])*ln(ln_part4)

        lrt = 2*(part1+part2-part3-part4)
        #p_value = r.pchisq(lrt,df=1,lower_tail=False)
        return lrt
   
class Region:
    def __init__(self,name):
        self.name = name
        self.snps = []
        self.individuals = None # Refer to Sheep in List 
        self.length = None
        self.positions = None
    
    #Creates the positions for each Region in form {position:[SNP1..SNPn]}
    def store_positions(self):
        snp_pos_dic = {}
        for snp in self.snps:
            if snp_pos_dic.has_key(snp.position):
                snp_pos_dic[snp.position].append(snp)
            else:
                snp_pos_dic[snp.position] = [snp]
        self.positions = snp_pos_dic
    
    #Creates the individuals for each Region in form {sheep:{position:[SNP1..SNPn]}}
    def store_individuals(self):
        sheep_position_snp = {}
        for snp in self.snps:
            if sheep_position_snp.has_key(snp.sheep):
                sheep_position_snp[snp.sheep][snp.position] = snp
            else:
                sheep_position_snp[snp.sheep] = {snp.position:snp}
        self.individuals = sheep_position_snp
    
    # Prints easy to read tabular output
    def print_pos(self):
        pos_list = []
        print "\n#ID",
        for position, snp in self.positions.items():
            pos_list.append(int(position))
        pos_list.sort()
        for i in pos_list:
            print i,     
        for sheep,snp_dict in self.individuals.items():
            print "\n",sheep.name,
            for i in pos_list:
                if snp_dict.has_key(i):
                    #print snp_dict[i].position,
                    print snp_dict[i].ref_base+"/"+snp_dict[i].mutant_base,
                else:
                    print "-",
    # Prints input format required by fastphase
    def print_fastphase(self):
        pos_list = []
        for position, snp in self.positions.items():
            pos_list.append(int(position))
        pos_list.sort()
        print len(self.individuals)
        print len(pos_list)
        print "P",
        for i in pos_list:
            print i,
        for sheep,snp_dict in self.individuals.items():
            line1 = ""
            line2 = ""
            print "\n#"+sheep.name
            for i in pos_list:
                #ref_base = self.positions[i][0].ref_base
                if snp_dict.has_key(i):
                    #print snp_dict[i].position,
                    if snp_dict[i].genotype != None:
                        #print sheep.sequences[self.name][i-1]
                        line1 += snp_dict[i].genotype[0]
                        line2 += snp_dict[i].genotype[1]
                    elif snp_dict[i].genotype == None:
                        line1 += "?"
                        line2 += "?"
                else:
                    if sheep.sequences[self.name][i-1] == "N":
                        line1 += "?"
                        line2 += "?"
                    else:
                        
                        line1 += "?"#sheep.sequences[self.name][i-1]
                        line2 += sheep.sequences[self.name][i-1]     
            print line1
            print line2,
            
        
        
        
        

