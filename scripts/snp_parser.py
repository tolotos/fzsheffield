#!/usr/bin/python

from shepherd import *


class SNPparser:
    def __init__(self, snp_file, sheep_objects,sheep_names):
        self.snp_file = snp_file
        self.sheep_objects = sheep_objects
        self.sheep_names = sheep_names
    
    def comparison(self):
        snp_file = open(self.snp_file, "r").readlines()
        for line in snp_file:
            line = line.split()
            region_name = line[0] 
            position = line[1]
            ref_base = line[2]
            counter = 0
            for i in range(5,147,6):
                bases = line[i:i+4]
                name = self.sheep_names[counter]
                current_snp = SNP(line[0],line[1],line[2])
                current_snp.basecount["A"] = float(bases[0])
                current_snp.basecount["C"] = float(bases[1])
                current_snp.basecount["G"] = float(bases[2])
                current_snp.basecount["T"] = float(bases[3])
                current_snp.sheep = self.sheep_objects[name]
                self.sheep_objects[name].snps.append(current_snp)
                if not counter == 22:
                    counter += 1
        return self.sheep_objects
    
    def individiual(snp_file):
        #Create individual sheep, parse the name of individial from filename
        filename = snp_file.split("/")[-1]
        sheepname = filename.split("_")[0]
        sheep = Sheep(sheepname)
        #Open each file and parse content
        snp_file = open(snp_file, "r").readlines()
        for line in snp_file:    
            line = line.rstrip().split()
            #each line in format:
            #Ref_Seq;Ref_Position;Ref_Base;Mutant_base;HQ_Depth;A;C;G;T
            if not line[0] == "Ref_Seq":
                #Create a SNP object and fill it with content
                current_snp = SNP(line[0],line[1],line[2],line[3],line[4])
                current_snp.basecount["A"] = float(line[5])
                current_snp.basecount["C"] = float(line[6])
                current_snp.basecount["G"] = float(line[7])
                current_snp.basecount["T"] = float(line[8])
                current_snp.sheep = sheep
                sheep.snps.append(current_snp)
                #Assign the SNP to the individual it belongs to
        return sheep
