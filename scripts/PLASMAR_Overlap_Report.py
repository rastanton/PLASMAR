#!/usr/bin/env python3

import sys

def CP_AMR_PR_GC_List_Maker(AMR_gamma, PR_gamma):
    """Makes a list of three elements: CPs, AMRs, and PRs and the contigs each is found on"""
    Genes = [[], [], []]
    Contigs = [[], [], []]
    f = open(AMR_gamma, 'r')
    String1 = f.readline()
    for line in f:
        Gene = line.split()[0]
        Allele = Gene.split('__')[-1]
        Contig = line.split()[1]
        if ('CARBAPENEM' in Gene):
            Genes[0].append(Allele)
            Contigs[0].append(Contig)
        elif ('CARBAPENEM' in Gene) == False:
            Genes[1].append(Allele)
            Contigs[1].append(Contig)
    f.close()
    f = open(PR_gamma, 'r')
    String1 = f.readline()
    for line in f:
        Gene = line.split()[0]
        Allele = Gene.split('__')[-1]
        Contig = line.split()[1]
        Genes[2].append(Allele)
        Contigs[2].append(Contig)
    f.close()
    Out = [Genes, Contigs]
    return Out

def PR_CP_Contigs(Gene_Contig_List):
    PR_Contig_Same = []
    for entry in Gene_Contig_List[1][2]:
        if (entry in Gene_Contig_List[1][0]):
            PR_Contig_Same.append(entry)
    return PR_Contig_Same

def CP_AMR_PR_Overlaps(AMR_gamma1, PR_gamma1, AMR_gamma2, PR_gamma2):
    Info1 = CP_AMR_PR_GC_List_Maker(AMR_gamma1, PR_gamma1)
    Info2 = CP_AMR_PR_GC_List_Maker(AMR_gamma2, PR_gamma2)
    Overlap_PR_1 = PR_CP_Contigs(Info1)
    Overlap_PR_2 = PR_CP_Contigs(Info2)
    Out = [[], [], []]
    for entry in range(3):
        for gene in Info1[0][entry]:
            if (gene in Info2[0][entry]) and (gene in Out[entry]) == False:
                Out[entry].append(gene)
    PR_Overlaps = []
    for entry in Out[2]:
        Pos1 = Info1[0][2].index(entry)
        Pos2 = Info2[0][2].index(entry)
        if (Info1[1][2][Pos1] in Overlap_PR_1) and (Info2[1][2][Pos2] in Overlap_PR_2):
            entry = entry + '*'
        elif ((Info1[1][2][Pos1] in Overlap_PR_1) and (Info2[1][2][Pos2] in Overlap_PR_2) == False) or ((Info1[1][2][Pos1] in Overlap_PR_1) == False and (Info2[1][2][Pos2] in Overlap_PR_2)):
            entry = entry + '^'
        PR_Overlaps.append(entry)
    Out[2] = PR_Overlaps
    All = []
    for entry in Out:
        entry.sort()
        New_Cat = []
        for allele in entry:
            New = allele.split('_')[0]
            if allele[-1] == '*':
                New = New + '*'
            if allele[-1] == '^':
                New = New + '^'
            New_Cat.append(New)
        All.append(New_Cat)
    return All

def PLASMAR_CP_AMR_PR_Overlap_Writer(input_file, output_file, AR_Folder, PR_Folder):
    Name_Cleaner(input_file)
    f = open(input_file, 'r')
    Out = open(output_file, 'w')
    String1 = f.readline()
    Out.write('ID1\tID2\tCarbapenemase\tAR Overlaps\tPlasmid Replicon Overlaps\tPlasmids Compared\tPlasmid Overlap Score\n')
    for line in f:
        List1 = line.split('\t')
        Info = CP_AMR_PR_Overlaps(AR_Folder + List1[0] + '_AMR.gamma', PR_Folder + List1[0] + '_PR.gamma', AR_Folder + List1[1] + '_AMR.gamma', PR_Folder + List1[1] + '_PR.gamma')
        Info_Line = ''
        if '*' in List1[2]:
            CPs = Info[0]
            CPs_Measured = List1[2].split(',')
            Out_CPs = []
            for entry in CPs_Measured:
                if (CPs_Measured in CPs):
                    Out_CPs.append(entry)
                else:
                    Out_CPs.append(entry[0:-1])
            Info[0] = Out_CPs
        else:
            Info[0] = ['No shared carbapenemase']
        for entry in Info:
            New = ','.join(entry)
            Info_Line += New + '\t'
        Out.write(List1[0] + '\t' + List1[1] + '\t' + Info_Line + List1[24] + '\t' + List1[40])
    f.close()
    Out.close()
    Corrector(output_file)

def File_Expander(input_file):
    f = open(input_file, 'r')
    Lines = []
    for line in f:
        Lines.append(line)
    f.close()
    Out = open(input_file, 'w')
    Out.write(Lines[0])
    for line in Lines[1:]:
        List1 = line.split('\t')
        New_Line = [List1[1], List1[0]]
        New_Line += List1[2:]
        New_Line = ('\t').join(New_Line)
        Out.write(line)
        Out.write(New_Line)
    Out.close()

def Name_Cleaner(input_file):
    f = open(input_file, 'r')
    Lines = []
    for line in f:
        Lines.append(line)
    f.close()
    Out = open(input_file, 'w')
    Out.write(Lines[0])
    for line in Lines[1:]:
        List1 = line.split('\t')
        ID1 = List1[0].split('/')[-1][0:-17]
        ID2 = List1[1].split('/')[-1][0:-17]
        List1[0] = ID1
        List1[1] = ID2
        Line = '\t'.join(List1)
        Out.write(Line)
    Out.close()
    
def End_Rounder(input_file):
    f = open(input_file, 'r')
    Lines = []
    for line in f:
        Lines.append(line)
    f.close()
    Out = open(input_file, 'w')
    Out.write(Lines[0])
    for line in Lines[1:]:
        List1 = line.split('\t')
        Correct = str(round(float(List1[-1]), 2))
        if Correct == '0.0':
            Correct = '0.00'
        List1[-1] = Correct
        Line = '\t'.join(List1)
        Out.write(Line + '\n')
    Out.close()

def Corrector(input_file):
    File_Expander(input_file)
    End_Rounder(input_file)
    Individual_CP_Splitter(input_file)

def Individual_CP_Splitter(input_file):
    f = open(input_file, 'r')
    Lines = []
    CPs = []
    String1 = f.readline()
    Header = String1
    for line in f:
        CP = line.split('\t')[2]
        if CP == 'No shared carbapenemase':
            continue
        else:
            if (CP in CPs):
                Pos = CPs.index(CP)
                Lines[Pos].append(line)
            else:
                CPs.append(CP)
                Lines.append([line])
    f.close()
    for entry in range(len(CPs)):
        Genes = CPs[entry].split(',')
        Name = '_'.join(Genes)
        Out = open(input_file[0:-4] + '_' + Name + '.txt', 'w')
        Out.write(Header)
        for line in Lines[entry]:
            Out.write(line)
        Out.close()

PLASMAR_CP_AMR_PR_Overlap_Writer('Overlap_Comp_Probs.txt', 'Overlap_Comp_Report.txt', 'AMR/', 'PR/')
Individual_CP_Splitter('Overlap_Comp_Report.txt')
