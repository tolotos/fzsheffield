#!/usr/bin/python

from optparse import OptionParser
import os
import re
import glob
from shepherd import *
from fasta import *
from snp_parser import *

#=======================================================================
#Command line options===================================================
#=======================================================================

usage = 'usage: %prog [options]'
desc="""%prog takes a folder containing txt files with SNP information \\
        for a given individual"""
        
cloptions = OptionParser(usage = usage, description=desc)
cloptions.add_option('-f', '--folder', dest = 'folder',
    help = 'Required folder', metavar='FILE',)
cloptions.add_option('-s', '--sequences', dest = 'sequences',
    help = 'Required folder', metavar='FILE',)
(options, args) = cloptions.parse_args()

#=======================================================================
#Locate all FASTA files in the given directory
fasta_files = []
fasta = options.sequences
for infile in glob.glob(os.path.join(fasta, '*.fasta')):
    fasta_files.append(infile)

#Parser for individual fasta file, utilising Class Fasta
def parse_sequences(file):
    fasta = Fasta(file)
    fasta.load()
    return fasta

#=======================================================================

#Function to create regions from sheep objects
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
#=======================================================================

#Start to create Sheep objects, names hardcoded, because not provided
#in snp_comparison_file
snp_file = options.folder
sheep_objects = {}
sheep_names = ["5421","5466","5467","5475","5477","5484",
         "5527","5528","5538","5554","5723","5732",
         "5734","5826","5994","6016","6048","Blue51",
         "Green417","IB008","IB017","IB021",
         "SBF1"]
for name in sheep_names:
    sheep_objects[name]=(Sheep(name))

snp_parser = SNPparser(snp_file,sheep_objects,sheep_names)

sheep_objects = snp_parser.comparison()
regions = extract_regions(sheep_objects)

print len(regions)
#=======================================================================
#Assign Sequences to sheep
for file in fasta_files:
    name = file.split("_")[-1]
    name = name.split(".")[0]
    fasta = parse_sequences(file)
    for sheep in sheep_objects.values():
        if sheep.name == name:
            sheep.sequences = fasta.sequences
#=======================================================================

#os.mkdir("fastphase_input")
#os.chdir("fastphase_input")
#for name, region in regions.items():
    #region_file = open(region.name+".genotype","w")
    #region.store_positions()
    #region.store_individuals()
    ##print region.individuals
    #for snp in region.snps:
        #sorted_bases = sorted(snp.basecount.iteritems(), 
                                  #key=operator.itemgetter(1), reverse=True)
        #if sum(snp.basecount.values()) == 0:
            #snp.genotype = None
        #else:
            #lr = snp.calc_genotype()
            #if lr <= -3.84:
                #snp.genotype = [sorted_bases[0][0],sorted_bases[1][0]]
            #elif lr >= 3.84:
                #snp.genotype = [sorted_bases[0][0],sorted_bases[0][0]]
            #else:
                #snp.genotype = None
    #reg_print = region.print_fastphase()
    
    #region_file.write("\n"+str(reg_print[0]))
    #region_file.write("\n"+str(reg_print[1]))
    #region_file.write("\n"+str(reg_print[2])+" ")
    #for i in reg_print[3]:
        #region_file.write(str(i)+" ")
    #for name, indiv  in reg_print[4].items():
        #region_file.write("\n"+indiv[0])
        #region_file.write("\n"+indiv[1])
        #region_file.write("\n"+indiv[2])
    #region_file.close()
