#!/usr/bin/python

import operator, pysam, math, rpy
     
class Sheep:
    def __init__(self,name):
        self.name = name
        self.breed = None
        self.numreads = None
        self.regions = [] # List of all Regions 
        self.snps = []
        self.sequences = None # {Region:sequence}
        
        
    def has_snp(self, regio_name, position):
        for region in self.regions:
            if region.name == regio_name:
                if region.snps.has_key(position):
                    return region.snp[position]
                else:
                    return False
            else:
                return False

class Position:
    #Ref_Seq;Ref_Position;Ref_Base;Mutant_base;HQ_Depth;A;C;G;T
    def __init__(self,ref_seq,ref_position,ref_base):
        self.position = int(ref_position)
        self.region = ref_seq
        self.mutant_base = None
        self.ref_base = ref_base
        self.hq_depth = None
        self.basecount = {"A":0,"C":0,"G":0,"T":0,"N":0}
        self.sheep = None
        #SNP genotype is stored as [(base,count)], len ==1 --> homozygote
        #len == 2 --> heterozygote
        #None == Genotype unknown!
        self.genotype = ()
        self.snp = False
        self.lrt = None

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
        self.lrt = lrt
        return lrt
   
class Region:
    def __init__(self,name):
        self.name = name
        self.snps = {}
        self.bases = {}
        self.sheep = None 
        self.length = None
        self.positions = None
        self.ref_seq = ""
        self.consensus = ""

    def get_bases(self,bam,site): 
        samfile = pysam.Samfile(bam, "rb")
        for pile in samfile.pileup(self.name, 0, self.length):
            pos = Position(self.name,pile.pos,self.ref_seq[pile.pos])
            pos.hq_depth = pile.n
            for pileupreads in pile.pileups:
                cur_base = pileupreads.alignment.seq[pileupreads.qpos]
                pos.basecount[cur_base] += 1
            if pile.pos == site:
                self.bases[pile.pos] = pos

           

    def get_snps(self,cutoff,bam):
        samfile = pysam.Samfile(bam, "rb")
        counter = 0
        for pile in samfile.pileup(self.name, 0, self.length):
            while counter < pile.pos:
                self.consensus+="N"
                counter += 1
            else:
                counter +=1
            pos = Position(self.name,pile.pos,self.ref_seq[pile.pos])
            pos.hq_depth = pile.n
            for pileupreads in pile.pileups:
                cur_base = pileupreads.alignment.seq[pileupreads.qpos]
                pos.basecount[cur_base] += 1
                
            sorted_bases = sorted(pos.basecount.iteritems(), 
                                  key=operator.itemgetter(1), reverse=True)
            del pos.basecount["N"]
            if sum(pos.basecount.values()) == 0:
                pos.genotype = None
                self.consensus+= "N"
            else:
                lr = pos.calc_genotype()
                if lr <= -float(cutoff):
                    pos.genotype = [sorted_bases[0][0],sorted_bases[1][0]]
                    pos.snp = True
                    self.snps[pile.pos] = pos
                    snp_iupac = sorted_bases[0][0]+sorted_bases[1][0]
                    if snp_iupac == "AG":
                        self.consensus += "R"
                    elif snp_iupac == "GA":
                        self.consensus += "R"
                    elif snp_iupac == "CT":
                        self.consensus += "Y"
                    elif snp_iupac ==  "TC":
                        self.consensus += "Y"
                    elif snp_iupac == "GT":
                        self.consensus += "K"
                    elif snp_iupac == "TG":
                        self.consensus += "K"
                    elif snp_iupac == "AC":
                        self.consensus += "M"
                    elif snp_iupac == "CA":
                        self.consensus += "M"
                    elif snp_iupac == "GC":
                        self.consensus += "S"
                    elif snp_iupac == "CG":
                        self.consensus += "S"
                    elif snp_iupac == "AT":
                        self.consensus += "W"
                    elif snp_iupac == "TA":
                        self.consensus += "W"
                elif lr >= float(cutoff):
                    pos.genotype = [sorted_bases[0][0],sorted_bases[0][0]]
                    if sorted_bases[0][0] != self.ref_seq[pile.pos]:
                        pos.snp = True
                        self.snps[pile.pos] = pos
                        self.consensus += sorted_bases[0][0]
                    else:
                        self.consensus += (sorted_bases[0][0])
                else:
                    pos.genotype = None
                    self.consensus+="N"

                    
                    
    # Prints input format required by fastphase
    def print_fastphase(self):
        print_list = []
        pos_list = []
        indiv = {}
        for position, snp in self.positions.items():
            pos_list.append(int(position))
        pos_list.sort()
        print_list.append(len(self.individuals))
        print_list.append(len(pos_list))
        print_list.append("P")
        print_list.append(pos_list)
        for sheep,snp_dict in self.individuals.items():
            line1 = ""
            line2 = ""
            indiv[sheep.name] = ["#"+sheep.name]
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
            indiv[sheep.name].append(line1)
            indiv[sheep.name].append(line2)
        print_list.append(indiv)
        return print_list
        


    #Creates the positions for each Region in form {position:[SNP1..SNPn]}
    def store_positions(self):
        snp_pos_dic = {}
        for snp in self.snps:
            if snp_pos_dic.has_key(snp.position):
                snp_pos_dic[snp.position].append(snp)
            else:
                snp_pos_dic[snp.position] = [snp]
        self.positions = snp_pos_dic


class Catalog:
    
    def __init__(self,name):
        self.name = name
        self.sites = {}
        self.consensus_seq = {} # {sheep_name: consensus ...}
        self.snpcount = None


    def init_sites(self,sheep,regio_name):
        for region in sheep.regions:
                if region.name == regio_name:
                    for position, snp in region.snps.items():
                        if self.sites.has_key(position):
                            self.sites[position][sheep.name] = snp
                        else:
                            self.sites[position] = {sheep.name:snp}
                    self.consensus_seq[sheep.name] = region.consensus


    def add_missing_genotypes(self,sheep):
        for pos,site in self.sites.items():
            if site.has_key(sheep.name):
                continue
            else:
                cur_consensus = self.consensus_seq[sheep.name]
                if len(cur_consensus) <= pos:
                    posi = Position(self.name,pos,"N")
                    posi.genotype = None
                else:
                    posi = Position(self.name,pos,cur_consensus[pos])
                    posi.genotype = (cur_consensus[pos],cur_consensus[pos])
                self.sites[pos][sheep.name] = posi
 
 
 
 
 
 
