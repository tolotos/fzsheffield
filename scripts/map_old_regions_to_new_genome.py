#!/usr/bin/python

from optparse import OptionParser
from fasta import Fasta
from Bio import SeqIO
#=======================================================================
#Command line options===================================================
#=======================================================================
usage = 'usage: %prog [options]'
desc="""%prog takes a fasta file, and return the first 200 bases of
the sequence"""
        
cloptions = OptionParser(usage = usage, description=desc)
cloptions.add_option('-f', '--fasta', dest = 'file',
    help = 'Input fasta', metavar='FILE',)
(options, args) = cloptions.parse_args()


#=======================================================================


file = options.file



for seq_record in SeqIO.parse(file, "fasta"):
    
    counter = 0
    start = 10142330
    stop = start+5000
    print ">%s_%s_%s" %(seq_record.id, start, stop)
    sequence = seq_record.seq[start:stop]
    for i in range(60,len(sequence),60):
        print sequence[counter:i]
        counter += 60
