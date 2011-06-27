#!/usr/bin/python

from optparse import OptionParser
from fasta import Fasta

#=======================================================================
#Command line options===================================================
#=======================================================================
usage = 'usage: %prog [options]'
desc="""%prog takes a fasta file, and capitalises all letters in
the sequence"""
        
cloptions = OptionParser(usage = usage, description=desc)
cloptions.add_option('-f', '--fasta', dest = 'file',
    help = 'Input fasta', metavar='FILE',)
(options, args) = cloptions.parse_args()


#=======================================================================

file = options.file

fasta = Fasta(file)
fasta.load()
fasta.capitalise()

for name, sequence in fasta.sequences.items():
    print ">"+name
    counter = 0
    for i in range(60,len(sequence),60):
        print sequence[counter:i]
        counter += 60
