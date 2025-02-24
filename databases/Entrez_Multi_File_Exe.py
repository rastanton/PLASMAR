import subprocess
import sys
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

###Takes in a file with entrez IDs and start positions and generates individual fastas with start positions
###Usage: python Entrez_Multi_File_Exe.py All_CP_Plasmid_Info.txt
##Written by: Rich Stanton (rstanton@cdc.gov/github.com/rastanton)

def Entrez_Fasta_File(accession, fasta):
    subprocess.call('efetch -db nucleotide -id ' + accession + ' -mode text -format fasta > ' + fasta, shell=True)

def Entrez_Multi_File(input_file, fasta_name):
    f = open(input_file, 'r')
    IDs = []
    for line in f:
        ID = line.split()[0][0:-6]
        IDs.append(ID)
    f.close()
    Current = []
    Count = 1
    for entry in IDs:
        if len(Current) == 600:
            ID = ','.join(Current)
            Entrez_Fasta_File(ID, fasta_name + str(Count) + '.fasta')
            Count += 1
            Current = [entry]
        else:
            Current.append(entry)
    if len(Current) != 0 or len(Current) != 600:
        ID = ','.join(Current)
        Entrez_Fasta_File(ID, fasta_name + str(Count) + '.fasta')
        
def Contig_Writer_All(input_fasta):
    """Writes out a fasta file of all the genes"""
    Genome = list(SeqIO.parse(input_fasta, 'fasta'))
    for gene in Genome:
        Output = open(gene.id + '.fasta', 'w')
        SeqIO.write(gene, Output, 'fasta')
        Output.close()

def Circular_Maker(gene_fasta, start):
    """Moves a circular fasta to start at the new start position"""
    gene = list(SeqIO.parse(gene_fasta, 'fasta'))
    Out = open(gene_fasta, 'w')
    New_Sequence = str(gene[0].seq[start:]) + str(gene[0].seq[0:start])
    new_gene = SeqRecord(
        Seq(New_Sequence),
        id=gene[0].id,
        description=gene[0].description,
    )
    SeqIO.write(new_gene, Out, 'fasta')

def Circular_File(input_file):
    f = open(input_file, 'r')
    for line in f:
        ID = line.split()[0]
        Start = int(line.split()[-1])
        Circular_Maker(ID, Start)
    f.close()
    
Entrez_Multi_File(sys.argv[1], 'PLASMAR_DB_')

List1 = glob.glob('PLASMAR_DB_*.fasta')
for file in List1:
    Contig_Writer_All(file)
    
Circular_File(sys.argv[1])
