#!/usr/bin/python

from optparse import OptionParser
import os
import re
import glob
from shepherd import *

#=============================================================================
#Command line options=========================================================
#=============================================================================

usage = 'usage: %prog [options]'
desc="""%prog takes a folder containing txt files with SNP information for a \\
        given individual"""
        
cloptions = OptionParser(usage = usage, description=desc)
cloptions.add_option('-f', '--folder', dest = 'folder',
    help = 'Required folder', metavar='FILE',)
(options, args) = cloptions.parse_args()

#=============================================================================

#Locate all SNP files in the given directory
snp_files = []
path = options.folder
for infile in glob.glob(os.path.join(path, '*.txt')):
    snp_files.append(infile)
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
            current_snp.basecount["A"] = line[5]
            current_snp.basecount["C"] = line[6]
            current_snp.basecount["G"] = line[7]
            current_snp.basecount["T"] = line[8]
            current_snp.sheep = sheep
            sheep.snps.append(current_snp)
            #Assign the SNP to the individual it belongs to
    return sheep

def extract_regions(sheep_objects):
    regions = {}
    for sheep in sheep_objects:
        for snp in sheep.snps:
            if regions.has_key(snp.region):
                regions[snp.region].snps.append(snp)
            else:
                regions[snp.region] = Region(snp.region)
                regions[snp.region].snps.append(snp)
    return regions

#=============================================================================

sheep_objects = []
#Iterate over all individual snp files
for snp_file in snp_files:
    sheep_objects.append(parse_snps(snp_file))

regions = extract_regions(sheep_objects)

# Print file format TO DO MOVE INTO FUNCTION
#for name, region in regions.items():
#    print "#",region.name, len(region.snps)
#    for snp in region.snps:
#        print snp.sheep.name, snp.position, snp.ref_base, snp.mutant_base,\
#              snp.hq_depth, snp.basecounts
#=============================================================================


for name, region in regions.items():
    print "\n",region.name
    region.store_positions()
    region.store_individuals()
    #print region.individuals
    region.print_fastphase()
            
        
            
            
            
            
            
            
            
            
            
            
            
            
            
            
