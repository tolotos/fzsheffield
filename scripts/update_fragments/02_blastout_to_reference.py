#!/usr/bin/python

from optparse import OptionParser
from Bio import SeqIO
import subprocess, os

#=======================================================================
#Command line options===================================================
#=======================================================================
usage = 'usage: %prog [options]'
desc="""%prog takes a fasta file, and return the first 200 bases of
the sequence"""
        
cloptions = OptionParser(usage = usage, description=desc)
cloptions.add_option('-b', '--blastouput', dest = 'blastout', type="str",
    help = 'Blastoutput to parse', metavar='FILE',)
cloptions.add_option('-f', '--fasta', dest = 'fasta', type="str",
    help = 'Input fasta, containing all old Fragments', metavar='FILE',)
cloptions.add_option('-s', '--sheep', dest = 'sheep', type="str",
    help = 'Sheep genome to parse', metavar='FILE',)

(options, args) = cloptions.parse_args()
#=======================================================================

sheep = options.sheep
options.fasta
blastout = open(options.blastout, "r").readlines()

frag_pos = {}

for line in blastout:
    line = line.rstrip().split()
    if frag_pos.has_key(line[0]):
            if frag_pos[line[0]][1] > line[10]:
                #frag_pos[line[0]] = [int(line[8]),float(line[10]),line[1]] # line 8 or 9 depending on order in blastout
                frag_pos[line[0]] = [int(line[8]),float(line[10]),line[1]]
    else:
        frag_pos[line[0]] = [int(line[8]),float(line[10]),line[1]] #[pos,e-value,chromosome]

if os.path.isdir("mafft_input"):
    subprocess.call("rm -r ./mafft_input", shell=True)
subprocess.call("mkdir ./mafft_input", shell=True)

old_frag = {}
for seq_record in SeqIO.parse(options.fasta, "fasta"):
    seq_record.id= seq_record.id.replace("&","_") # mafft does not read &
    old_frag[seq_record.id] = seq_record.seq

for seq_record in SeqIO.parse(sheep, "fasta"):
    seq_record.id= seq_record.id.replace("&","_") # mafft does not read &
    #if os.path.dirname("./"): #!= "mafft_input":
    if os.path.exists("../mafft_input") == False:
        os.chdir("mafft_input")
    for frag_name, items in frag_pos.items():
        items[2]= items[2].replace("&","_") # mafft does not read &
        if items[2] == seq_record.id:
            frag_name = frag_name.replace("&","_") # mafft does not read &
            start= items[0]-100
            stop = items[0]+len(old_frag[frag_name])+200
            map_file = open(items[2]+frag_name+".fa","w")
            map_file.write(">%s_%s_from_%s_to_%s\n" % (seq_record.id,frag_name,start,stop))
            sequence = seq_record.seq[start:stop]
            sequence = sequence
            counter = 0
            for i in range(60,len(sequence),60):
                seq = sequence[counter:i]+"\n"
                map_file.write(str(seq))
                counter += 60
            #map_file.write(">%s\n" % (frag_name))
            #sequence = old_frag[frag_name]
            ##sequence = sequence[::-1].complement()
            #counter = 0
            #for i in range(60,len(sequence),60):                
                #seq = sequence[counter:i]+"\n"
                #map_file.write(str(seq))
                #counter += 60
            map_file.close()
