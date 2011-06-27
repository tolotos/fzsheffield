#!/usr/bin/python

from optparse import OptionParser
from Bio import SeqIO
import os, glob,subprocess

#=======================================================================
#Command line options===================================================
#=======================================================================
usage = 'usage: %prog [options]'
desc="""%prog takes a fasta file, and return the first 200 bases of
the sequence"""
        
cloptions = OptionParser(usage = usage, description=desc)
cloptions.add_option('-f', '--folder', dest = 'folder',
    help = 'Input fasta', metavar='FILE',)
(options, args) = cloptions.parse_args()


#=======================================================================

folder = options.folder

fasta_files = []
for infile in glob.glob(os.path.join(folder, '*.fa')):
    fasta_files.append(infile)
print fasta_files
if os.path.isdir("mafft_broken_output"):
    subprocess.call("rm -r ./mafft_broken_output", shell=True)
subprocess.call("mkdir ./mafft_broken_output", shell=True)


for fasta in fasta_files:
    name = fasta.split("/")[-1]
    cmd = "mafft %s > ./mafft_broken_output/%s.aln" % (fasta, name)
    subprocess.call(cmd, shell=True)

