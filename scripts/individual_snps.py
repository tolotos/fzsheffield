#!/usr/bin/python

from optparse import OptionParser
import os
import re
import glob
from shepherd import *
from fasta import *

#=============================================================================
#Command line options=========================================================
#=============================================================================

usage = 'usage: %prog [options]'
desc="""%prog takes a folder containing txt files with SNP information for a \\
        given individual"""
        
cloptions = OptionParser(usage = usage, description=desc)
cloptions.add_option('-f', '--folder', dest = 'folder',
    help = 'Required folder', metavar='FILE',)
cloptions.add_option('-s', '--sequences', dest = 'sequences',
    help = 'Required folder', metavar='FILE',)
(options, args) = cloptions.parse_args()

#=============================================================================

#Locate all SNP files in the given directory
snp_files = []
snp_file = options.folder
#for infile in glob.glob(os.path.join(path, '*.txt')):
    #snp_files.append(infile)
#Locate all FASTA files in the given directory
fasta_files = []
fasta = options.sequences
for infile in glob.glob(os.path.join(fasta, '*.fasta')):
    fasta_files.append(infile)

#=============================================================================
#Parser for individual SNP files
def parse_snps(snp_file):
    
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

#Parser for individual fasta file, utilising Class Fasta
def parse_sequences(file):
    fasta = Fasta(file)
    fasta.load()
    return fasta
    
def extract_regions(sheep_objects):
    regions = {}
    for sheep in sheep_objects.values():
        for snp in sheep.snps:
            if regions.has_key(snp.region):
                regions[snp.region].snps.append(snp)
            else:
                regions[snp.region] = Region(snp.region)
                regions[snp.region].snps.append(snp)
    return regions

sheep_objects = {}
sheep = ["5421","5466","5467","5475","5477","5484",
         "5527","5528","5538","5554","5723","5732",
         "5734","5826","5994","6016","6048","Blue51",
         "Green417","IB008","IB017","IB021",
         "SBF1"]
for name in sheep:
    sheep_objects[name]=(Sheep(name))

def parse_snp_comparison(snp_file):
    snp_file = open(snp_file, "r").readlines()
    counter = 0
    for line in snp_file:
        line = line.split()
        region_name = line[0] 
        position = line[1]
        ref_base = line[2]
        counter = 0
        for i in range(5,147,6):
            bases = line[i:i+4]
            name = sheep[counter]
            current_snp = SNP(line[0],line[1],line[2])
            current_snp.basecount["A"] = float(bases[0])
            current_snp.basecount["C"] = float(bases[1])
            current_snp.basecount["G"] = float(bases[2])
            current_snp.basecount["T"] = float(bases[3])
            current_snp.sheep = sheep_objects[name]
            sheep_objects[name].snps.append(current_snp)
            if not counter == 22:
                counter += 1
parse_snp_comparison(snp_file)


#=============================================================================


#Iterate over all individual snp files
#for snp_file in snp_files:
    #sheep_objects.append(parse_snps(snp_file))

#Assign Sequences to sheep
for file in fasta_files:
    name = file.split("_")[-1]
    name = name.split(".")[0]
    fasta = parse_sequences(file)
    for sheep in sheep_objects.values():
        if sheep.name == name:
            sheep.sequences = fasta.sequences

regions = extract_regions(sheep_objects)


#Print file format TO DO MOVE INTO FUNCTION
#for name, region in regions.items():
    #print "#",region.name, len(region.snps)
    #for snp in region.snps:
        #print snp.sheep.name, snp.position, snp.ref_base, snp.mutant_base,\
              #snp.hq_depth, snp.basecount
#=============================================================================


for name, region in regions.items():
    print "\n"+region.name
    region.store_positions()
    region.store_individuals()
    #print region.individuals
    for snp in region.snps:
        sorted_bases = sorted(snp.basecount.iteritems(), 
                                  key=operator.itemgetter(1), reverse=True)
        if sum(snp.basecount.values()) == 0:
            snp.genotype = None
        else:
            lr = snp.calc_genotype()
            if lr <= -3.84:
                snp.genotype = [sorted_bases[0][0],sorted_bases[1][0]]
            elif lr >= 3.84:
                snp.genotype = [sorted_bases[0][0],sorted_bases[0][0]]
            else:
                snp.genotype = None
    region.print_fastphase()
