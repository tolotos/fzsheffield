#!/usr/bin/python

from optparse import OptionParser
from Bio import SeqIO

#=======================================================================
#Command line options===================================================
#=======================================================================
usage = 'usage: %prog [options]'
desc="""%prog takes a fasta file, and return the first 200 bases of
the sequence"""
        
cloptions = OptionParser(usage = usage, description=desc)
cloptions.add_option('-l', '--length', dest = 'length', type="int",
    help = 'Length of start fragment', metavar='FILE',)
cloptions.add_option('-f', '--fasta', dest = 'file',
    help = 'Input fasta', metavar='FILE',)
(options, args) = cloptions.parse_args()


#=======================================================================

file = options.file
stop = options.length

for seq_record in SeqIO.parse(file, "fasta"):
    print ">%s" %(seq_record.id)
    counter = 0
    sequence = seq_record.seq[:stop]
    for i in range(60,len(sequence),60):
        print sequence[counter:i]
        counter += 60
