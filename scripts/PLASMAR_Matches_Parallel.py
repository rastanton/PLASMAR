#!/usr/bin/env python3

import sys
import os
import glob
import subprocess
from statistics import mean
from multiprocessing import Pool
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import gzip
import pickle
import operator

def Overlap(a,b):
    return (a[0] >= b[0] and a[1] <= b[1]) or (b[0] >= a[0] and b[1] <= a[1]) or (a[0] <= b[1] and a[1] >= b[0]) or (b[0] <= a[1] and b[1] >= a[0]) or ((a[1] + 1) == b[0]) or ((b[1] + 1) == a[0])

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
    for entry in input_list:
        Focus = entry
        for entry2 in input_list:
            if Overlap(Focus, entry2) == True:
                Focus = Overlap_Extender(Focus, entry2)
        if (Focus in Out) == False:
            Out.append(Focus)
    return Out

def Recursive_Overlap(input_list):
    Current = input_list
    New = Multiple_Overlap_Extenders(Current)
    while New != Current:
        Current = New
        New = Multiple_Overlap_Extenders(Current)
    return New

def Match_Start_Stop_Finder(PSL_Line):
    """Finds the start and stop for the contig and gene"""
    List1 = PSL_Line.split('\t')
    Output = []
    Block_Lengths = List1[18].split(',')[0:-1]
    Blocks = []
    for lengths in Block_Lengths:
        Blocks.append(int(lengths))
    Block_Lengths = Blocks
    Blocks = []
    Gene_Starts = List1[-1].split(',')[0:-1]
    for lengths in Gene_Starts:
        Blocks.append(int(lengths))
    Gene_Starts = Blocks
    Blocks = []
    Genome_Starts = List1[-2].split(',')[0:-1]
    for lengths in Genome_Starts:
        Blocks.append(int(lengths))
    Genome_Starts = Blocks
    if List1[8] == '-':
        Genome_Start = int(List1[10]) - int(List1[12])
        Genome_End = int(List1[10]) - int(List1[11])
    else:
        Genome_Start = int(List1[11])
        Genome_End = int(List1[12])
    Genome_Length = int(List1[10])
    Gene_Start = int(List1[15])
    Gene_End = int(List1[16])
    Gene_Length = int(List1[14])
    Gene_Starts.append(Gene_End)
    Genome_Starts.append(Genome_End)       
    Genome_Start_Stop = [Genome_Start, Genome_End]
    Output.append(Genome_Start_Stop)
    Gene_Start_Stop = [Gene_Start, Gene_End]
    Output.append(Gene_Start_Stop)
    Output.append(Block_Lengths)
    Output.append(Genome_Starts)
    Output.append(Gene_Starts)
    return Output

def Gene_Match_Block_Finder(input_psl):
    """Makes a list of all the start and stops of the overlap"""
    f = open(input_psl, 'r')
    All = []
    for line in f:
        Info = Match_Start_Stop_Finder(line)
        Blocks = Info[2]
        Starts = Info[4][0:-1]
        for entry in range(len(Blocks)):
            New = [Starts[entry], Starts[entry] + Blocks[entry]]
            All.append(New)
    return All

def Gene_Match_Combinbed_Block(input_psl):
    """Makes non-overlapping start and stop block matches"""
    Info = Gene_Match_Block_Finder(input_psl)
    Combined = Recursive_Overlap(Info)
    return Combined

def Total_Length(input_list):
    """Takes in a list of non-overlapping match blocks and returns the total length"""
    Total = 0
    for entry in input_list:
        tot = entry[1] - entry[0]
        Total += tot
    return Total

def Total_Match_Length(input_psl):
    Blocks = Gene_Match_Combinbed_Block(input_psl)
    Length = Total_Length(Blocks)
    return Length

def Mutant_Counter(input_psl):
    """Counts all the mutants in a psl"""
    Count = 0
    f = open(input_psl, 'r')
    for line in f:
        Mutants = line.split()[1]
        Count += int(Mutants)
    f.close()
    return Count

def Maximum_Match(input_psl):
    """Finds the largest match block from a psl"""
    Largest = 0
    f = open(input_psl, 'r')
    for line in f:
        Large = line.split()[0]
        if int(Large) > Largest:
            Largest = int(Large)
    f.close()
    return Largest

def Match_Info(input_psl):
    f = open(input_psl, 'r')
    String1 = f.readline()
    if String1 == '':
        Out = '0\t0\t0\t0\t0\t0\t0\t0\t0'
        return Out
    Plasmid_Length = int(String1.split()[14])
    f.close()
    Max_Match = Maximum_Match(input_psl)
    Mutants = Mutant_Counter(input_psl)
    Match_Length = Total_Match_Length(input_psl)
    Percent_Overlap = (Match_Length) / (Plasmid_Length)
    Match_Total = Match_Length - Mutants
    Percent_Match = (Match_Total) / (Match_Length)
    Weighted_Match = (Percent_Overlap) * (Percent_Match)
    Contig_Info = Best_Contig(input_psl)
    Best_Contig_Percent = Contig_Info[1]
    Longest_Contig_Match = Contig_Info[0]
    Out = str(Longest_Contig_Match) + '\t' + str(Best_Contig_Percent) + '\t' + str(Percent_Match) + '\t' + str(Percent_Overlap) + '\t' + str(Weighted_Match) + '\t' + str(Max_Match) + '\t' + str(Match_Total) + '\t' + str(Match_Length) + '\t' + str(Plasmid_Length)
    return Out

def Match_Line_2(input_psl, assembly_name, plasmid_line):
    line = Match_Info(input_psl)
    Out = assembly_name + '\t' + plasmid_line + '\t' + line
    return Out

def Match_Line_Maker_3(assembly, plasmid, plasmid_line, AMR_gamma, PR_gamma, output):
    subprocess.call('blat' + ' ' + plasmid + ' '  + assembly + ' -noHead PSL/' + output + '.psl', shell=True)
    assembly_name = assembly.split('/')[-1]
    New = CP_AMR_PR_Percents_Weighted(plasmid_line, AMR_gamma, PR_gamma, 'PSL/' + output + '.psl')
    Add = str((New[0])) + '\t' + str((New[1])) + '\t' + str((New[2]))
    Line = plasmid_line + '\t' + Add
    Out = Match_Line_2('PSL/' + output + '.psl', assembly_name, Line)
    return Out

def Multi_Match_Maker_Line_List(input_list):
    ID = input_list[1].split('\t')[0]
    plasmid = Upper + '/databases/Plasmid_DB/' + ID
    assembly_name = input_list[0].split('/')[-1]
    output = assembly_name + '__' + ID
    Line = Match_Line_Maker_3(input_list[0], plasmid, input_list[1], input_list[2], input_list[3], output)
    return Line

def Multi_Match_Maker_Parallel(assembly, plasmid_line_list, AMR_gamma, PR_gamma, output_file):
    """Makes a file with lines for the match line output"""
    List1 = []
    for entry in plasmid_line_list:
        List1.append([assembly, entry, AMR_gamma, PR_gamma])
    Out = open(output_file, 'w')
    Out.write('Assembly\tPlasmid\tLength\tCP\tAMR\tPR\tCP%\tAMR%\tPR%\tCP_Contig%\tAMR_Contig%\tPR_Contig%\tMax_Contig_Match\t%Max_Contig\t%_Match\t%_Overlap\t%Weighted_Match\tMax_Match\tTotal_Match\tTotal_Match_Length\tPlasmid_Length\tMax_Overlap_%\tCP_PR_Contig_%\tAR_PR_Contig_%\n')
    with Pool() as pool:
        Lines = pool.map(Multi_Match_Maker_Line_List, List1)
    AR_List = AR_Genes(AMR_gamma)
    PR_List = PR_Genes(PR_gamma)
    Lines = Max_Matcher(Lines)
    for line in Lines:
        Line = Gene_Stars(line, AR_List, PR_List)
        CP_Contig_Overlap = CP_Overlap_Percent(AMR_gamma, PR_gamma, Line)
        AR_Contig_Overlap = AR_Overlap_Percent(AMR_gamma, PR_gamma, Line)
        Out.write(Line + '\t' + str(CP_Contig_Overlap) + '\t' + str(AR_Contig_Overlap) + '\n')
    Out.close()

def Max_Matcher(input_lines):
    Max = []
    Genes = []
    for line in input_lines:
        List1 = line.split('\t')
        CP = List1[3]
        if (CP in Genes) == False:
            Genes.append(CP)
            Max.append(int(List1[19]))
        else:
            Pos = Genes.index(CP)
            if Max[Pos] < int(List1[19]):
                Max[Pos] = int(List1[19])
    Out = []
    for line in input_lines:
        List1 = line.split('\t')
        CP = List1[3]
        Pos = Genes.index(CP)
        Percent = int(List1[19]) / Max[Pos]
        Out.append(line + '\t' + str(Percent))
    return Out
     
def AR_Genes(AR_GAMMA):
    f = open(AR_GAMMA, 'r')
    genes = []
    String1 = f.readline()
    for line in f:
        Gene = line.split()[0]
        allele = Gene.split('__')[-1]
        if (allele in genes) == False:
            genes.append(allele)
    f.close()
    return genes

def PR_Genes(PR_GAMMA):
    f = open(PR_GAMMA, 'r')
    genes = []
    String1 = f.readline()
    for line in f:
        allele = line.split()[0]
        if (allele in genes) == False:
            genes.append(allele)
    f.close()
    return genes

def Gene_Overlap(target_genes, sample_genes):
    Out = []
    for gene in target_genes:
        if (gene in sample_genes):
            Out.append(gene + '*')
        else:
            Out.append(gene)
    return Out

def Gene_Stars(input_line, AR, PR):
    List1 = input_line.split('\t')
    CP = List1[3].split(',')
    AMR = List1[4].split(',')
    Reps = List1[5].split(',')
    Name = List1[0][0:-6]
    CP_Out = Gene_Overlap(CP, AR)
    AMR_Out = Gene_Overlap(AMR, AR)
    Reps_Out = Gene_Overlap(Reps, PR)
    List1[3] = ','.join(CP_Out)
    List1[4] = ','.join(AMR_Out)
    List1[5] = ','.join(Reps_Out)
    Line = '\t'.join(List1)
    return Line

def AR_Contig_Extractor(input_GAMMA):
    Genes = []
    Contigs = []
    f = open(input_GAMMA, 'r')
    String1 = f.readline()
    for line in f:
        allele = line.split()[0]
        gene = allele.split('__')[-1]
        contig = line.split()[1]
        if (gene in Genes):
            Pos = Genes.index(gene)
            Contigs[Pos].append(contig)
        else:
            Genes.append(gene)
            Contigs.append([contig])
    f.close()
    return [Genes, Contigs]

def PR_Contig_Extractor(input_GAMMA):
    Genes = []
    Contigs = []
    f = open(input_GAMMA, 'r')
    String1 = f.readline()
    for line in f:
        gene = line.split()[0]
        contig = line.split()[1]
        if (gene in Genes):
            Pos = Genes.index(gene)
            Contigs[Pos].append(contig)
        else:
            Genes.append(gene)
            Contigs.append([contig])
    f.close()
    return [Genes, Contigs]

def Asterisk_Finder(input_string):
    List1 = input_string.split(',')
    Out = []
    if List1 == ['']:
        return Out
    else:
        for items in List1:
            if items[-1] == '*':
                Out.append(items[0:-1])
    return Out

def List_Overlap_Binary(List1, List2):
    Score = 0
    for entry in List1:
        if (entry in List2):
            Score = 1
    return Score

  
def CP_Overlap_Percent(AR_GAMMA, PR_GAMMA, PLASMAR_Line):
    List1 = PLASMAR_Line.split('\t')
    CP = Asterisk_Finder(List1[3])
    PR = Asterisk_Finder(List1[5])
    AR_Contigs = AR_Contig_Extractor(AR_GAMMA)
    PR_Contigs = PR_Contig_Extractor(PR_GAMMA)
    CP_Contigs = []
    for entry in CP:
        Pos = AR_Contigs[0].index(entry)
        CP_Contigs.append(AR_Contigs[1][Pos])
    Local_PR_Contigs = []
    for entry in PR:
        Pos = PR_Contigs[0].index(entry)
        Local_PR = PR_Contigs[1][Pos]
        for entry2 in Local_PR:
            Local_PR_Contigs.append(entry2)
    Count = 0
    for entry in CP_Contigs:
        Score = List_Overlap_Binary(entry, Local_PR_Contigs)
        Count += Score
    Percent = Count / len(CP)
    return Percent

def AR_Overlap_Percent(AR_GAMMA, PR_GAMMA, PLASMAR_Line):
    List1 = PLASMAR_Line.split('\t')
    CP = Asterisk_Finder(List1[4])
    PR = Asterisk_Finder(List1[5])
    AR_Contigs = AR_Contig_Extractor(AR_GAMMA)
    PR_Contigs = PR_Contig_Extractor(PR_GAMMA)
    CP_Contigs = []
    for entry in CP:
        Pos = AR_Contigs[0].index(entry)
        CP_Contigs.append(AR_Contigs[1][Pos])
    Local_PR_Contigs = []
    for entry in PR:
        Pos = PR_Contigs[0].index(entry)
        Local_PR = PR_Contigs[1][Pos]
        for entry2 in Local_PR:
            Local_PR_Contigs.append(entry2)
    Count = 0
    for entry in CP_Contigs:
        Score = List_Overlap_Binary(entry, Local_PR_Contigs)
        Count += Score
    if len(CP) == 0:
        return 0
    else:
        Percent = Count / len(CP)
    return Percent

def Repeat_Remover(any_list):
    """Removes repeats for any list"""
    new_list = []
    for items in any_list:
        if (items in new_list) == False:
            new_list.append(items)
    return new_list

def List_Matcher(List1, List2):
    Correct1 = Repeat_Remover(List1)
    Correct2 = Repeat_Remover(List2)
    if Correct1 == ['']:
        Correct1 = []
    if Correct2 == ['']:
        Correct2 = []
    Total = len(Correct1)
    Count = 0
    for entry in Correct1:
        if (entry in Correct2):
            Count += 1
    Out = [Count, Total]
    return Out

def List_Match_Percent(List1, List2):
    Info = List_Matcher(List1, List2)
    if Info[1] == 0:
        Percent = 1
    else:
        Percent = (Info[0]) / (Info[1])
    return Percent

def Contig_Match_Block_Finder(input_psl, contig):
    """Makes a list of all the start and stops of the overlap for a specific contig"""
    f = open(input_psl, 'r')
    All = []
    for line in f:
        ID = line.split('\t')[9]
        if ID == contig:
            Info = Match_Start_Stop_Finder(line)
            Blocks = Info[2]
            Starts = Info[3][0:-1]
            for entry in range(len(Blocks)):
                New = [Starts[entry], Starts[entry] + Blocks[entry]]
                All.append(New)
    return All

def Contig_Match_Combined_Block(input_psl, contig):
    """Makes non-overlapping start and stop block matches"""
    Info = Contig_Match_Block_Finder(input_psl, contig)
    Combined = Recursive_Overlap(Info)
    return Combined

def Contig_Finder_gamma(input_gamma, gene_list):
    All = []
    f = open(input_gamma, 'r')
    String1 = f.readline()
    for line in f:
        gene = line.split('\t')[0]
        gene = gene.split('__')[-1]
        if (gene in gene_list):
            Contig = line.split('\t')[1]
            All.append(Contig)
    f.close()
    return All

def Contig_Length_Finder(input_psl, contig):
    Length = 0
    f = open(input_psl, 'r')
    for line in f:
        ID = line.split('\t')[9]
        if ID == contig:
            Length = line.split('\t')[10]
            Length = int(Length)
            break
    return Length

def Contig_Coverage_List(input_psl, contig):
    """Makes a list with the coverage and the contig length"""
    Coverage = [0, 0]
    List1 = Contig_Match_Combined_Block(input_psl, contig)
    Contig_Length = Contig_Length_Finder(input_psl, contig)
    Cov_Length = Total_Length(List1)
    if Contig_Length != 0:
        Coverage = [Cov_Length, Contig_Length]
    return Coverage

def All_Contigs(input_psl):
    f = open(input_psl, 'r')
    Contigs = []
    for line in f:
        ID = line.split('\t')[9]
        if (ID in Contigs) == False:
            Contigs.append(ID)
    f.close()
    return Contigs

def All_Contig_Coverage(input_psl):
    Contigs = All_Contigs(input_psl)
    Data = []
    for contig in Contigs:
        Info = Contig_Coverage_List(input_psl, contig)
        Data.append([contig, Info[0], Info[1]])
    Data.sort(key=lambda x: x[1], reverse=True)
    return Data

def Best_Contig(input_psl):
    Data = All_Contig_Coverage(input_psl)
    Percent = Data[0][1] / Data[0][2]
    return [Data[0][1], Percent]

def Genes_Contig_Coverage(gene_list, gamma, input_psl):
    if gene_list == ['']:
        return 1
    All = []
    Contigs = Contig_Finder_gamma(gamma, gene_list)
    Contigs = Repeat_Remover(Contigs)
    Coverage_Length = 0
    Contig_Length = 0
    for contig in Contigs:
        Cov = Contig_Coverage_List(input_psl, contig)
        Coverage_Length += Cov[0]
        Contig_Length += Cov[1]
    if Contig_Length == 0:
        All_Weighted = 0
    else:
        All_Weighted = Coverage_Length / Contig_Length
    return All_Weighted

def CP_AMR_PR_Percents(plasmid_line, AMR_gamma, PR_gamma, plasmid_psl):
    """Takes in a plasmid line and calculates the % CP, %AMR, and % plasmid markers for a plasmid match"""
    if plasmid_line[-1] == '\n':
        plasmid_line = plasmid_line[0:-1]
    List1 = plasmid_line.split('\t')
    CP = List1[2].split(',')
    CP = Repeat_Remover(CP)
    AMR = List1[3].split(',')
    AMR = Repeat_Remover(AMR)
    PR = List1[4].split(',')
    PR = Repeat_Remover(PR)
    CP_Cov = Genes_Contig_Coverage(CP, AMR_gamma, plasmid_psl)
    AMR_Cov = Genes_Contig_Coverage(AMR, AMR_gamma, plasmid_psl)
    PR_Cov = Genes_Contig_Coverage(PR, PR_gamma, plasmid_psl)
    Out = [CP_Cov, AMR_Cov, PR_Cov]
    return Out

def CP_AMR_PR_Percents_Weighted(plasmid_line, AMR_gamma, PR_gamma, plasmid_psl):
    """Takes in a plasmid line and calculates the % CP, %AMR, and % plasmid markers for a plasmid match weighted by the fraction of overlapping genes"""
    Raw_Contig_Matches = CP_AMR_PR_Percents(plasmid_line, AMR_gamma, PR_gamma, plasmid_psl)
    Raw_Fraction = Plasmid_Match_Finder_4_Line(AMR_gamma, PR_gamma, plasmid_line)
    Out = []
    for entry in range(3):
        Weighted = Raw_Contig_Matches[entry] * Raw_Fraction[entry]
        Out.append(Weighted)
    return Out

def Plasmid_Match_Finder_4_Line(AMR_gamma, PR_gamma, plasmid_line):
    """Eliminates the AMR alignment requirement"""
    Info = CP_AMR_PR_List_Maker(AMR_gamma, PR_gamma)
    List1 = plasmid_line[0:-1].split('\t')
    CP = List1[2].split(',')
    AMR = List1[3].split(',')
    PR = List1[4].split(',')
    CP_Match = List_Match_Percent(CP, Info[0])
    AMR_Match = List_Match_Percent(AMR, Info[1])
    PR_Match = List_Match_Percent(PR, Info[2])
    Matches = [CP_Match, AMR_Match, PR_Match]
    return Matches

def Plasmid_Match_Finder_5(AMR_gamma, PR_gamma, plasmid_file, run_matches):
    """Eliminates the AMR alignment requirement"""
    Info = CP_AMR_PR_List_Maker(AMR_gamma, PR_gamma)
    Matches = []
    f = open(plasmid_file, 'r')
    String1 = f.readline()
    for line in f:
        List1 = line[0:-1].split('\t')[0:-1]
        line = '\t'.join(List1)
        CP = List1[2].split(',')
        AMR = List1[3].split(',')
        PR = List1[4].split(',')
        CP_Match = List_Match_Percent(CP, Info[0])
        AMR_Match = List_Match_Percent(AMR, Info[1])
        PR_Match = List_Match_Percent(PR, Info[2])
        if run_matches == 0:
            if CP_Match == 1 and PR_Match > 0 and (AMR_Match > 0 or AMR == ['']):
                line = line + '\t' + str(CP_Match) + '\t' + str(AMR_Match) + '\t' + str(PR_Match)
                Matches.append(line)
        elif run_matches == 1:
            if CP_Match == 1:
                line = line + '\t' + str(CP_Match) + '\t' + str(AMR_Match) + '\t' + str(PR_Match)
                Matches.append(line)
    return Matches

def Plasmid_Match_Finder_Test(AMR_gamma, PR_gamma, plasmid_file):
    """Eliminates the AMR alignment requirement"""
    Info = CP_AMR_PR_List_Maker(AMR_gamma, PR_gamma)
    Matches = []
    f = open(plasmid_file, 'r')
    String1 = f.readline()
    for line in f:
        List1 = line[0:-1].split('\t')
        CP = List1[6].split(',')
        AMR = List1[7].split(',')
        PR = List1[8].split(',')
        CP_Match = List_Match_Percent(CP, Info[0])
        AMR_Match = List_Match_Percent(AMR, Info[1])
        PR_Match = List_Match_Percent(PR, Info[2])
        if CP_Match > 0:
            line = line[0:-1] + '\t' + str(CP_Match) + '\t' + str(AMR_Match) + '\t' + str(PR_Match)
            Matches.append(line)
    return Matches

def CP_AMR_PR_List_Maker(AMR_gamma, PR_gamma):
    """Makes a list of three elements: CPs, AMRs, and PRs"""
    Out = [[], [], []]
    f = open(AMR_gamma, 'r')
    String1 = f.readline()
    for line in f:
        Gene = line.split()[0]
        Allele = Gene.split('__')[-1]
        if ('CARBAPENEM' in Gene):
            Out[0].append(Allele)
        else:
            Out[1].append(Allele)
    f.close()
    f = open(PR_gamma, 'r')
    String1 = f.readline()
    for line in f:
        Gene = line.split()[0]
        Allele = Gene.split('__')[-1]
        Out[2].append(Allele)
    f.close()
    return Out

def Plasmid_Report_Maker_4(assembly, AMR_gamma, PR_gamma, plasmid_file, output_file, run_matches):
    Plasmid_List = Plasmid_Match_Finder_5(AMR_gamma, PR_gamma, plasmid_file, run_matches)
    Multi_Match_Maker_Parallel(assembly, Plasmid_List, AMR_gamma, PR_gamma, output_file)

def Array_Maker(input_list):
    Out = np.array(input_list)
    return Out

def Array_Maker_Plasmar_Float(input_file, start, stop):
    """Makes an array from an input Plasmar for positions given"""
    List1 = []
    f = open(input_file, 'r')
    String1 = f.readline()
    for line in f:
        Line_List = line[0:-1].split('\t')
        Line = []
        for entry in Line_List[start:stop]:
            Line.append(float(entry))
        List1.append(Line)
    f.close()
    Out = Array_Maker(List1)
    return Out

def PLASMAR_Single_Line_List(input_plasmar):
    Out = []
    f = open(input_plasmar, 'r')
    String1 = f.readline()
    for line in f:
        List1 = line[0:-1].split('\t')
        Out.append(List1)
    f.close()
    return Out

def GZ_Pickle(input_pickle):
    with gzip.open(input_pickle, 'rb') as file:
        Model = pickle.load(file)
    return Model

def RF_Strat_Kfold_Predictions_single_Pickle_Scored(PLASMAR, model):
    """Takes in an X and y array and generates a file with the predicted values based on an input training set and unknown data"""
    f = open(PLASMAR, 'r')
    String1 = f.readline()
    String1 = f.readline()
    f.close()
    if String1 == '':
        Out = open(PLASMAR[0:-4] + '_Prob.txt', 'w')
        Out.write('Assembly\tPlasmid\tLength\tCP\tAMR\tPR\tCP%\tAMR%\tPR%\tCP_Contig%\tAMR_Contig%\tPR_Contig%\tMax_Contig_Match\t%Max_Contig\t%_Match\t%_Overlap\t%Weighted_Match\tMax_Match\tTotal_Match\tTotal_Match_Length\tPlasmid_Length\tMax_Overlap_%\tCP_PR_Contig_%\tAR_PR_Contig_%\tProbability\n')
        Out.close()
    else:
        Array = Array_Maker_Plasmar_Float(PLASMAR, 6, 22)
        Lines = PLASMAR_Single_Line_List(PLASMAR)
        rf_model = GZ_Pickle(model)
        Predictions = []
        Data = rf_model.predict_proba(Array)
        Predictions.append(Data)
        for entry in range(len(Lines)):
            Lines[entry].append(Predictions[0][entry][1])
        Lines.sort(key=operator.itemgetter(24), reverse=True)
        Out = open(PLASMAR[0:-4] + '_Prob.txt', 'w')
        Out.write('Assembly\tPlasmid\tLength\tCP\tAMR\tPR\tCP%\tAMR%\tPR%\tCP_Contig%\tAMR_Contig%\tPR_Contig%\tMax_Contig_Match\t%Max_Contig\t%_Match\t%_Overlap\t%Weighted_Match\tMax_Match\tTotal_Match\tTotal_Match_Length\tPlasmid_Length\tMax_Overlap_%\tCP_PR_Contig_%\tAR_PR_Contig_%\tProbability\n')
        for line in Lines:
            Line = '\t'.join(line[0:-1])
            Out.write(Line + '\t' + str(line[-1]) + '\n')
        Out.close()
    subprocess.call('rm ' + PLASMAR, shell=True)

script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
Upper = os.path.dirname(script_directory)
subprocess.call('mkdir PSL', shell=True)
subprocess.call('mkdir PLASMAR_Matches', shell=True)
run_matches = sys.argv[1]

List1 = glob.glob('*.fasta')
for entry in List1:
    AMR = 'AMR/' + entry[0:-6] + '_AMR.gamma'
    Reps = 'PR/' + entry[0:-6] + '_PR.gamma'
    Plasmar = 'PLASMAR_Matches/' + entry[0:-6] + '_Matches.txt'
    Plasmid_Report_Maker_4(entry, AMR, Reps, Upper + '/databases/All_CP_Plasmid_Info.txt', Plasmar, int(run_matches))
    RF_Strat_Kfold_Predictions_single_Pickle_Scored(Plasmar, Upper + '/models/PLASMAR_Match_Val_SKL_1.6.1.pkl.gz')
 
