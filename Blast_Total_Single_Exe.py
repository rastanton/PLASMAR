import sys
import glob
from Bio import SeqIO
from decimal import *
getcontext().prec = 5
import subprocess

##Usage:> python Blast_Total_Alignment_Exe.py Blast_Glob query_fasta Output_File
##Requires ncbi-blast+

def Max_Contig_Length(genome, contig_list):
    Max = 0
    Genome = SeqIO.to_dict(SeqIO.parse(genome, 'fasta'))
    for contig in contig_list:
        Length = len(Genome[contig].seq)
        if Length > Max:
                     Max = Length
    return Length

def Contig_Lister(input_blast):
    Out = []
    f = open(input_blast, 'r')
    for line in f:
        Contig = line.split('\t')[1]
        if (Contig in Out) == False:
            Out.append(Contig)
    return Out

def Gene_Length_Finder(input_gene):
    gene = list(SeqIO.parse(input_gene, 'fasta'))
    Length = len(gene[0].seq)
    return Length

def Overlap(a,b):
    return (a[0] >= b[0] and a[1] <= b[1]) or (b[0] >= a[0] and b[1] <= a[1]) or (a[0] <= b[1] and a[1] >= b[0]) or (b[0] <= a[1] and b[1] >= a[0])

def Overlap_Extender(a, b):
    Out = a
    if (Overlap(a, b)) == True:
        if (a[0] <= b[0]):
            Out[0] = a[0]
        else:
            Out[0] = b[0]
        if a[1] <= b[1]:
            Out[1] = b[1]
        else:
            Out[1] = a[1]
    return Out

def Multiple_Overlap_Extenders(input_list):
    Out = []
##    Add = 0
    for entry in input_list:
        Focus = entry
        for entry2 in input_list:
            if Overlap(Focus, entry2) == True:
                Focus = Overlap_Extender(Focus, entry2)
##                if Focus != entry:
##                    Add = 1
        if (Focus in Out) == False:
            Out.append(Focus)
##    return [Out, Add]
    return Out

def Recursive_Overlap(input_list):
    Current = input_list
    New = Multiple_Overlap_Extenders(Current)
    while New != Current:
        Current = New
        New = Multiple_Overlap_Extenders(Current)
    return New

def Total_Coverage(input_list, length):
    Lists = Recursive_Overlap(input_list)
    Total = 0
    for entry in Lists:
        Total += (entry[1] - (entry[0] - 1))
    Percent = Total/length
    return Percent

def Blast_Coverage(input_blast, length):
    Positions = []
    f = open(input_blast, 'r')
    for line in f:
        List1 = line.split()
        Position = [int(List1[6]), int(List1[7])]
        Positions.append(Position)
    f.close()
    Out = Total_Coverage(Positions, length)
    return Out

def Blast_Length(input_blast):
    Positions = []
    f = open(input_blast, 'r')
    for line in f:
        List1 = line.split()
        Position = [int(List1[6]), int(List1[7])]
        Positions.append(Position)
    f.close()
    Lists = Recursive_Overlap(Positions)
    Total = 0
    for entry in Lists:
        Total += (entry[1] - (entry[0] - 1))
    return Total
    
def Lean_Blast_Output5(input_blast_file, input_subject, input_query, length):
    """Reads in a blast file and calculates the % identity and % coverage"""
    Total_Coverage = Blast_Coverage(input_blast_file, length)
    Total_Percent = 0
    Total_Length = 0
##    Length = Gene_Length_Finder(input_query)
    f = open(input_blast_file, 'r')
    String1 = f.readline()
    while String1 != '':
        List1 = String1.split()
        if len(List1) > 10:
            Total_Percent = Total_Percent + (float(List1[2]) * int(List1[3]))
            Total_Length = Total_Length + int(List1[3])
        String1 = f.readline()
    f.close()
    Percent = Decimal(Total_Percent) / Decimal(Total_Length)
    Coverage_Length = Blast_Length(input_blast_file)
    Total_Coverage = Decimal(Total_Coverage) * 100
    Contig_List = Contig_Lister(input_blast_file)
    Max_Contig = Max_Contig_Length(input_subject, Contig_List)
    Out = input_subject + '\t' + input_query + '\t' + str(Percent) + '\t' + str(Total_Coverage) + '\t' + str(Coverage_Length) + '\t' + str(Max_Contig) + '\t' + str(length)
    return Out

def Blast_Single(subject, query, output_blast_file):
    """Makes a blast and returns the match data"""
    Length = Gene_Length_Finder(query)
    subprocess.call('blastn -query ' + query + ' -subject '  + subject + ' -out ' + output_blast_file + ' -outfmt 6', shell=True)
    Info = Lean_Blast_Output5(output_blast_file, subject, query, Length)
    return Info

Blast_Single(sys.argv[1], sys.argv[2], sys.argv[3])
