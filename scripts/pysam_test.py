#!/usr/bin/python

import pysam
samfile = pysam.Samfile("5421_sorted.bam", "rb" )
fastafile = pysam.Fastafile( "5421_ref.fasta" )
pileup_iter = samfile.pileup( stepper = "samtools", fastafile = fastafile )
snpcaller = pysam.SNPCaller(pileup_iter)
print snpcaller( "OAR10", 100 )
samfile.close()
