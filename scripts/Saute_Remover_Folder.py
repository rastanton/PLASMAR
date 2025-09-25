#!/usr/bin/env python3

from Bio import SeqIO
import sys
import glob

#Usage: Saute_Remover.py Fasta_Folder
#Removes contigs from NCBI assemblies made with Saute

def Saute_Remover(input_genome):
    Genome = list(SeqIO.parse(input_genome, 'fasta'))
    Out = open(input_genome, 'w')
    for contig in Genome:
        Description = contig.description
        if ('guided' in Description) == False:
            SeqIO.write(contig, Out, 'fasta')
    Out.close()

Folder = sys.argv[1]
Files = glob.glob(Folder + '/*.fasta')
for file in Files:
    Saute_Remover(file)
