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

os.mkdir("region_fasta")
os.chdir("region_fasta")
for name, region in regions.items():
    region.store_positions()
    region.store_individuals()
    region_file = open(region.name+".fa","w")
    for sheep in region.individuals.keys():
        region_file.write(">"+sheep.name)
        region_file.write("\n"+sheep.sequences[region.name]+"\n")
    region_file.close()
        
