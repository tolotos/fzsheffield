#!/usr/bin/python 

from optparse import OptionParser
import pysam, subprocess, glob, os

#=======================================================================
#Command line options===================================================
#=======================================================================

usage = 'usage: %prog [options]'
desc="""%prog takes a folder containing sam files and applies different
        filters, resulting in indexed bam files"""
        
cloptions = OptionParser(usage = usage, description=desc)
cloptions.add_option('-f', '--folder', dest = 'folder',
    help = 'Required folder', metavar='FOLDER',)
cloptions.add_option('-o', '--output', dest = 'output',
    help = 'Specify output folder', metavar='FOLDER',)
(options, args) = cloptions.parse_args()

#=======================================================================

#Locate all SAM files in the given directory
sam_files = []
bam_files = []
sam = options.folder
output = options.output
for infile in glob.glob(os.path.join(sam, '*.sam')):
    sam_files.append(infile)

print sam
print "Starting to parse each single sam file and convert it to bam."
for sam in sam_files:
    name = sam.split("/")[-1]
    name = name.split(".")[0]
#=======================================================================
    cmd = "samtools view -bS -o "+output+name+".bam "+ sam
    bam_files.append(name+"_sorted.bam")
    samtools_view = subprocess.call(cmd, shell=True)
    samtools_view = 0
    if samtools_view == 0:
        print name+".sam has been exported to bam format successfully..."
    else:
        print "Export to bam format failed!"
    
#=======================================================================
    cmd = "samtools sort "+output+name+".bam "+output+name+"_sorted"
    print cmd
    print "Sorting bam file...",
    samtools_sort = subprocess.call(cmd, shell=True)
    if samtools_sort == 0:
        print "Done!"
    else:
        print "Sorting bam failed!"
#=======================================================================
    cmd = "samtools index "+output+name+"_sorted.bam"
    samtools_index = subprocess.call(cmd, shell=True)
    if samtools_index == 0:
        print "Done indexing!"
    else:
        print "Indexing bam failed!"
#=======================================================================
print "Starting to filter each sorted bam file.."
for bam in bam_files:
    name = bam.split("_")[0]
    cmd = "samtools rmdup "+output+bam+" "+output+name+"_rmdump_sorted.bam"
    samtools_filter = subprocess.call(cmd, shell=True)
    if samtools_filter == 0:
        print "Filterd out PCR duplicates!"
    else:
        print "Filtering PCR duplicates failed!"
    print "Reindexing...",
    cmd = "samtools index "+output+name+"_rmdump_sorted.bam"
    samtools_index = subprocess.call(cmd, shell=True)
    if samtools_index == 0:
        print "Done indexing!"
    else:
        print "Indexing bam failed!"

