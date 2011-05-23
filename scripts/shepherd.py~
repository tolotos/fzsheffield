#!/usr/bin/python

class Sheep:
    def __init__(self,name):
        self.name = name
        self.breed = None
        self.numreads = None
        self.regions = None # List of all Regions 
        self.snps = []

class SNP:
    #Ref_Seq;Ref_Position;Ref_Base;Mutant_base;HQ_Depth;A;C;G;T
   def __init__(self,ref_seq,ref_position,ref_base,mutant_base,hq_depth):
        self.position = int(ref_position)
        self.region = ref_seq
        self.mutant_base = mutant_base
        self.ref_base = ref_base
        self.hq_depth =int(hq_depth)
        self.basecount = {"A":0,"C":0,"G":0,"T":0}
        self.sheep = None

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
                ref_base = self.positions[i][0].ref_base
                if snp_dict.has_key(i):
                    #print snp_dict[i].position,
                    line1 += "?"#snp_dict[i].ref_base
                    line2 += snp_dict[i].mutant_base
                else:
                    line1 += "?"#ref_base
                    line2 += ref_base     
            print line1
            print line2,

