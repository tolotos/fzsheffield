#!/usr/bin/python

from optparse import OptionParser
import shepherd, snp_parser, os, re, glob, pysam
from fasta import *
from shepherd  import *

#=======================================================================
#Command line options===================================================
#=======================================================================

usage = 'usage: %prog [options]'
desc="""%prog takes a folder containing txt files with SNP information 
        for a given individual"""
        
cloptions = OptionParser(usage = usage, description=desc)
cloptions.add_option('-f', '--folder', dest = 'folder',
    help = 'Required folder', metavar='FILE',)
cloptions.add_option('-r', '--reference', dest = 'reference',
    help = 'Required folder', metavar='FILE',)
(options, args) = cloptions.parse_args()

#=======================================================================
#Locate all sorted BAM files in the given directory
bam_files = []
bam = options.folder
for infile in glob.glob(os.path.join(bam, '*rmdump_sorted.bam')):
    bam_files.append(infile)
reference_seq = options.reference
#Parser for individual fasta file, utilising Class Fasta
def parse_fasta(file):
    fasta = Fasta(file)
    fasta.load()
    return fasta
fasta = parse_fasta(reference_seq)
regions = fasta.get_names()
sheep_objects = {}

#=======================================================================
for bam in bam_files:
    name = bam.split("/")[-1].split("_")[0]
    sheep_objects[name] = Sheep(name)
    for regio_name in regions:
        sheep_objects[name].regions.append(Region(regio_name))
        cur_regio = sheep_objects[name].regions[-1]
        cur_regio.ref_seq = fasta.sequences[cur_regio.name]
        cur_regio.length = len(cur_regio.ref_seq)
        cur_regio.get_snps(3.84,bam) # P=0.001 = 10.83; P=0.05 = 3.84
        #if len(cur_regio.snps) == 0:
            #print "Region %s has no valid snps!" % (cur_regio.name)
        #else:
            ##cur_regio.store_positions()
            #print "Region %s has been processed for sheep %s ..." % (cur_regio.name, name)
#=======================================================================


catalogs = []

for regio_name in regions:
        cur_catalog = Catalog(regio_name)
        for sheep in sheep_objects.values():
            cur_catalog.init_sites(sheep,regio_name)
        catalogs.append(cur_catalog)

for catalog in catalogs:
    for sheep in sheep_objects.values():
        catalog.add_missing_genotypes(sheep)
        
sheep_names = []
print "I", "ID",
for sheep in sheep_objects.values():
    print sheep.name, sheep.name,
    sheep_names.append(sheep.name)
print
for catalog in catalogs:
    for pos in sorted(catalog.sites.iterkeys()):
        print "M", catalog.name+"_"+str(pos),
        site = catalog.sites[pos]
        for sheep_name in sheep_names:
            if site[sheep_name].genotype == None:
                print "N", "N",
            else:
                print site[sheep_name].genotype[0],site[sheep_name].genotype[1],           
        print
