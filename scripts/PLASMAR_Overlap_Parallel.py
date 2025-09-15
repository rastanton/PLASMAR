#!/usr/bin/env python3

import sys
import os
import glob
import subprocess
from statistics import mean
from operator import itemgetter
from multiprocessing import Pool
from math import sqrt
import gzip
import pickle
from sklearn.ensemble import ExtraTreesClassifier
import numpy as np
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

def Overlap_Length(List1, List2):
    """Finds the overlap length of two lists"""
    Out = 0
    if Overlap(List1, List2):
        New = []
        if List1[0] >= List2[0]:
            New.append(List1[0])
        else:
            New.append(List2[0])
        if List1[1] <= List2[1]:
            New.append(List1[1])
        else:
            New.append(List2[1])
        Out = New[1] - New[0]
    return Out

def Total_Overlap(List1, List2):
    Out = 0
    for item in List1:
        for item2 in List2:
            Distance = Overlap_Length(item, item2)
            Out += Distance
    return Out

def PSL_Overlaps5(psl1, psl2):
    """Finds the overlap and match lengths of two psls"""
    List1 = Gene_Match_Combinbed_Block(psl1)
    List2 = Gene_Match_Combinbed_Block(psl2)
    Length1 = Total_Length(List1)
    Length2 = Total_Length(List2)
    Overlap = Total_Overlap(List1, List2)
    List3 = List1 + List2
    Combined = Recursive_Overlap(List3)
    Short = min([Length1, Length2])
    Long = max([Length1, Length2])
    Length_All = Total_Length(Combined)
    return [Overlap, Length_All, Short, Long]

def List_Overlaps(List1, List2):
    """Finds if there are overlaps in two lists"""
    Out = []
    for items in List1:
        if items in List2:
            Out.append(items)
    return Out

def Plasmid_CP_Extractor(PLASMAR):
    Out = []
    CPs = []
    f = open(PLASMAR, 'r')
    String1 = f.readline()
    for line in f:
        CP = line.split('\t')[3]
        plasmid = line.split()[1]
        if (CP in CPs):
            Pos = CPs.index(CP)
            Out[Pos].append(plasmid)
        else:
            CPs.append(CP)
            Out.append([plasmid])
    f.close()
    return [CPs, Out]

def Plasmid_Overlap_Data(PLASMAR1, PLASMAR2):
    Plasmids1 = Plasmid_CP_Extractor(PLASMAR1)
    Plasmids2 = Plasmid_CP_Extractor(PLASMAR2)
    Out = []
    for entry in range(len(Plasmids1[0])):
        for entry2 in range(len(Plasmids2[0])):
            if Plasmids1[0][entry] == Plasmids2[0][entry2]:
                Overlaps = List_Overlaps(Plasmids1[1][entry], Plasmids2[1][entry2])
                Out.append([Plasmids1[0][entry], Overlaps])
    return Out

def PLASMAR_ID(PLASMAR):
    f = open(PLASMAR, 'r')
    String1 = f.readline()
    String1 = f.readline()
    ID = String1.split()[0]
    return ID

def Fraction_Finder5(PLASMAR1, PLASMAR2, PSL_Folder):
    Overlap_Plasmids = Plasmid_Overlap_Data(PLASMAR1, PLASMAR2)
    ID1 = PLASMAR_ID(PLASMAR1)
    ID2 = PLASMAR_ID(PLASMAR2)
    Genes = []
    Overlaps = []
    for entry in Overlap_Plasmids:
        Genes.append(entry[0])
        Overlaps.append(entry[1])
    Data = []
    for entry in Overlaps:
        if entry == []:
            All = [[0, 0, 0, 0]]
        else:
            All = []
            for plasmid in entry:
                Info = PSL_Overlaps5(PSL_Folder + ID1 + '__' + plasmid + '.psl', PSL_Folder + ID2 + '__' + plasmid + '.psl')
                All.append(Info)
        Data.append(All)
    if len(Data) == 0:
        Data.append([0, 0, 0, 0])
    return Genes, Data

def Fraction_5_Averages(input_fraction_list):
    Out = []
    for entry in input_fraction_list:
        if entry == [0, 0, 0, 0]:
            New = [0, 0, 0]
        else:
            All = entry[0] / entry[1]
            Small = entry[0] / entry[2]
            Large = entry[0] / entry[3]
            New = [All, Small, Large]
        Out.append(New)
    return Out

def List_Tab(input_list):
    Out = ''
    for entry in input_list:
        Out = Out + str(entry) + '\t'
    return Out[0:-1]


def Fraction_List5(PLASMAR1, PLASMAR2, PSL_Folder):
    Plasmids = Fraction_Finder5(PLASMAR1, PLASMAR2, PSL_Folder)
    Lines = []
    for Info in Plasmids[1]:
        Percents = Fraction_5_Averages(Info)
        Max = str(max(list(map(itemgetter(0), Info))))
        Total = sum(list(map(itemgetter(0), Info)))
        All = sum(list(map(itemgetter(1), Info)))
        Short = sum(list(map(itemgetter(2), Info)))
        Long = sum(list(map(itemgetter(3), Info)))
        Scores = []
        for entry in [All, Short, Long]:
            if entry == 0:
                Scores.append(0)
            else:
                score = Total / entry
                Scores.append(score)
        Means = []
        for entry in range(3):
            Data = mean(list(map(itemgetter(entry), Percents)))
            Means.append(Data)
        Out = Max + '\t' + List_Tab(Scores) + '\t' + List_Tab(Means) + '\n'
        Lines.append(Out)
    return [Plasmids[0], Lines]

def Similarity_Coefficient(number1, number2):
    """Uses the square root similarity coefficient"""
    Sum1 = number1 * number1
    Sum2 = number2 * number2
    if Sum1 == 0 and Sum2 == 0:
        Coefficient = 1
    elif Sum1 == 0 or Sum2 == 0:
        Coefficient = 0
    else:
        Coefficient = min([Sum1, Sum2]) / sqrt(Sum1 * Sum2)
    return Coefficient

def List_Means(List1):
    Out = []
    for entry in range(len(List1[0])):
        Mean = mean(list(map(itemgetter(entry), List1)))
        Out.append(Mean)
    return Out

def List_Overlaps(List1, List2):
    """Finds if there are overlaps in two lists"""
    Out = []
    for items in List1:
        if items in List2:
            Out.append(items)
    return Out

def Plasmid_Extractor(PLASMAR):
    Out = []
    f = open(PLASMAR, 'r')
    String1 = f.readline()
    for line in f:
        plasmid = line.split()[1]
        Out.append(plasmid)
    f.close()
    return Out


def Tab_File_Maker(input_file):
    """Makes a list of file lines skipping line 1"""
    f = open(input_file, 'r')
    String1 = f.readline()
    Out = []
    for line in f:
        Info = line[0:-1].split('\t')
        Out.append(Info)
    f.close()
    return Out

def Float_Maker(input_list):
    Out = []
    for item in input_list:
        Out.append(float(item))
    return Out

def Repeat_Remover(any_list):
    """Removes repeats for any list"""
    new_list = []
    for items in any_list:
        if (items in new_list) == False:
            new_list.append(items)
    return new_list

def List_Overlap_Score(List1, List2):
    Score = 0
    for entry in List1:
        if (entry in List2):
            Score += 1
    Out = Score / len(List2)
    return Out
   
def Gene_Overlaps(genes1, genes2):
    G1 = []
    if (genes1 == ['']) and (genes2 == ['']):
        return 1
    elif (genes1 == ['']) or (genes2 == ['']):
        return 0
    else:
        for genes in genes1:
            if genes[-1] == '*' and (genes in G1) == False:
                G1.append(genes)
        G2 = []
        for genes in genes2:
            if genes[-1] == '*' and (genes in G2) == False:
                G2.append(genes)
        All = G1 + G2
        All = Repeat_Remover(All)
        Overlap = []
        for genes in G1:
            if (genes in Overlap) == False and (genes in G2):
                Overlap.append(genes)
        if len(All) == 0:
            Score = 1
        else:
            Score = List_Overlap_Score(Overlap, All)
    return Score

def Gene_Overlap_Scorer(Line1, Line2):
    Out = []
    List1 = Line1[3:6]
    List2 = Line2[3:6]
    for entry in range(3):
        New1 = List1[entry].split(',')
        New2 = List2[entry].split(',')
        Score = Gene_Overlaps(New1, New2)
        Out.append(Score)
    return Out


def Weighted_Mean(input_list):
    All = 0
    Total = 0
    for entry in input_list:
        Value = Similarity_Coefficient(entry[0], entry[1])
        Max = max(entry)
        Out = Max * Value
        All += Out
        Total += Max
    if Total == 0:
        Data = 0
    else:
        Data = All / Total
    return Data

def Weight_List_Means(input_list):
    All = []
    for entry in range(len(input_list[0][0])):
        All.append([])
    for entry in input_list:
        for Pos in range(len(All)):
            All[Pos].append([entry[0][Pos], entry[1][Pos]])
    Out = []
    for entry in All:
        Value = Weighted_Mean(entry)
        Out.append(Value)
    return Out


def Plasmid_Average_Similarity(PLASMAR1, PLASMAR2):
    Info = Plasmid_Overlap_Data(PLASMAR1, PLASMAR2)
    Lines = []
    CPs = []
    for entry in range(len(Info)):
        if Info[entry][1] == []:
            Line = Zero_Maker(28)
        else:
            Line = Plasmid_Average_Similarity_Single(PLASMAR1, PLASMAR2, Info[entry][1])
        CPs.append(Info[entry][0])
        Lines.append(Line)
    return CPs, Lines
    

def Plasmid_Average_Similarity_Single(PLASMAR1, PLASMAR2, Plasmids):
    List1 = Tab_File_Maker(PLASMAR1)
    List2 = Tab_File_Maker(PLASMAR2)
    Probability_Ranges = [0, 0, 0, 0]
    Probs_Sims = []
    Max_Sim = 0
    Gene_Comps = []
    Comps = []
    for plasmid in Plasmids:
        for line1 in List1:
            for line2 in List2:
                if line1[1] == plasmid and line2[1] == plasmid:
                    Genes = Gene_Overlap_Scorer(line1, line2)
                    vals1 = Float_Maker(line1[6:24])
                    vals2 = Float_Maker(line2[6:24])
                    Comps.append([vals1, vals2])
                    vals1 = float(line1[-1])
                    vals2 = float(line2[-1])
                    Gene_Comps.append(Genes)
                    Min = min([vals1, vals2])
                    if Min > Max_Sim:
                        Max_Sim = Min
                    if vals1 > 0 or vals2 > 0:
                        Probs_Sims.append([vals1, vals2])
                    break
    if len(Comps) > 0:
        Out = Weight_List_Means(Comps)
    else:
        Out = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    Gene_Means = List_Means(Gene_Comps)
    Out = Gene_Means + Out
    String1 = ''
    for entry in Out:
        String1 = String1 + str(entry) + '\t'
    if len(Probs_Sims) == 0:
        Mean_Sim = -1
    else:
        Mean_Sim = Weighted_Mean(Probs_Sims)
    Both = []
    for entry in Probs_Sims:
        if (entry[0] > 0) and (entry[1] > 0):
            Both.append(entry)
    if len(Probs_Sims) != 0:
        Fraction = len(Both) / len(Probs_Sims)
    else:
        Fraction = 0
    Overlap_Sim_Mean = Weighted_Mean(Both)
    String1 = String1 + str(len(Plasmids)) + '\t' + str(Fraction) + '\t' + str(Overlap_Sim_Mean) + '\t' + str(Mean_Sim) +'\t' + str(len(Both)) + '\t' + str(len(Probs_Sims)) + '\t' + str(Max_Sim)
    return String1

def Zero_Maker(count):
    Out = ''
    for entry in range(count):
        Out += '0\t'
    return Out[0:-1]

def PLASMAR_Overall_Similarity_Scores(PLASMAR1, PLASMAR2, PSL_Folder):
    Lines1 = Plasmid_Average_Similarity(PLASMAR1, PLASMAR2)
    Lines = []
    Lines2 = Fraction_List5(PLASMAR1, PLASMAR2, PSL_Folder)
    for entry in range(len(Lines1[0])):
        List1 = Lines1[1][entry].split()
        New1 = '\t'.join(List1[0:-1])
        Out = Lines1[0][entry] + '\t' + New1 + '\t' + Lines2[1][entry][0:-1] + '\t' + List1[-1] + '\n'
        Lines.append(Out)
    return Lines

def PLASMAR_Overall_Similarity_Scores_List(input_list):
    try:
        Info = PLASMAR_Overall_Similarity_Scores(input_list[0], input_list[1], input_list[2])
        Out = ''
        for entry in Info:
            Out = Out + input_list[0] + '\t' + input_list[1] + '\t' + entry
    except:
        Out = input_list[0] + '\t' + input_list[1] + '\tError\n'
    return Out

##def Comp_Maker(PLASMAR_Glob, PSL_Folder, output_file):
##    Pairs = Glob_Pairs_Maker(PLASMAR_Glob)
##    Out = open(output_file, 'w')
##    Out.write('CP%\tAMR%\tPR%\tCP%\tAMR%\tPR%\tCP_Contig%\tAMR_Contig%\tPR_Contig%\t%_Match\t%_Overlap\t%Weighted_Match\tMax_Match\tTotal_Match\tTotal_Match_Length\tProbability_Overlap\tNumber_of_overlaps\tLongest_Length\tAverage_Overlap\tPR CP Contig\tPR AR Contig\tPlasmids\tFraction_>0\tOverlap_>1_Sim_Score\tOverall_Sim_Score\tLen_Both_Score>0\tLen_Any_Score>0\tAll Overlap Fraction\tOverlap Fraction Min\tOverlap Fraction Max\tMean All Overlap\tMean Short Overlap\tMean Long Overlap\tMax Overlap\n')
##    for pair in Pairs:
##        Line = PLASMAR_Overall_Similarity_Scores(pair[0], pair[1], PSL_Folder)
##        Out.write(pair[0] + '\t' + pair[1] + '\t' + Line)
##    Out.close()

def List_Splitter(input_list, count):
    Out = []
    New = []
    for entry in input_list:
        if len(New) == count:
            Out.append(New)
            New = [entry]
        else:
            New.append(entry)
    Out.append(New)
    return Out

def Pair_Maker(input_list):
    Pairs = []
    for entry in range(0, len(input_list)):
        for entry2 in range(entry + 1, len(input_list)):
            Pairs.append([input_list[entry], input_list[entry2]])
    return Pairs

def Comp_Maker_Pool_Folder(PLASMAR_Matches_Folder, PSL_Folder, output_file):
    Pair_List = []
    List1 = glob.glob(PLASMAR_Matches_Folder + '/*_Matches_Prob.txt')
    Inputs = Pair_Maker(List1)
    for entry in Inputs:
        Pair_List.append([entry[0], entry[1], PSL_Folder])
    Pairs = List_Splitter(Pair_List, 1000)
    Out = open(output_file, 'w')
    Out.write('ID1\tID2\tCarbapenemase\tCP%\tAMR%\tPR%\tCP_Fraction_%\tAMR_Fratcion_%\tPR_Fraction_%\tCP_Contig%\tAMR_Contig%\tPR_Contig%\t%_Match\t%_Overlap\t%Weighted_Match\tMax_Match\tTotal_Match\tTotal_Match_Length\tProbability_Overlap\tNumber_of_overlaps\tLongest_Length\tAverage_Overlap\tPR CP Contig\tPR AR Contig\tPlasmids\tFraction_>0\tOverlap_>1_Sim_Score\tOverall_Sim_Score\tLen_Both_Score>0\tLen_Any_Score>0\tMax Overlap\tAll Overlap Fraction\tOverlap Fraction Min\tOverlap Fraction Max\tMean All Overlap\tMean Short Overlap\tMean Long Overlap\tMax Match Prob\n')
    for entry in Pairs:
        with Pool() as pool:
            Lines = pool.map(PLASMAR_Overall_Similarity_Scores_List, entry)
        for line in Lines:
            Out.write(line)
    Out.close()

def Array_Maker(input_list):
    Out = np.array(input_list)
    return Out

def Array_Maker_Plasmar_Float(input_file, start, stop):
    """Makes an array from an input Plasmar for positions given"""
    List1 = []
    Y_list = []
    f = open(input_file, 'r')
    String1 = f.readline()
    for line in f:
        Line_List = line[0:-1].split('\t')
        Line = []
        Y = 0
        if len(Line_List) > stop:
            for entry in Line_List[start:stop]:
                Line.append(float(entry))
            List1.append(Line)
            Y = 1
        Y_list.append(Y)
    f.close()
    Out = Array_Maker(List1)
    return Out, Y_list

def Array_Maker_Positions(input_PLASMAR_scores, pos_list):
    X_list = []
    Y_list = []
    f = open(input_PLASMAR_scores, 'r')
    String1 = f.readline()
    for line in f:
        List1 = line.split()
        X_line = List1
        X_correct = []
        Y = 0
        if len(X_line) > max(pos_list):
            for entry in range(len(X_line)):
                if (entry in pos_list):
                    X_correct.append(float(X_line[entry]))
            X_list.append(X_correct)
            Y = 1
        Y_list.append(Y)
    f.close()
    Out = []
    X_data = Array_Maker(X_list)
    Out = X_data
    return Out, Y_list

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

def Overlap_Pickle_Scored_Combined(PLASMAR, all_model, float_model):
    """Takes in an X and y array and generates a file with the predicted values based on an input training set and unknown data"""
    Pos = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 25, 26, 27, 31, 32, 33, 34, 35, 36]
    Float_Array, Y_Float = Array_Maker_Positions(PLASMAR, Pos)
    Array, Y_Array = Array_Maker_Plasmar_Float(PLASMAR, 3, 37)
    Lines = []
    if len(Array) > 0:
        Lines_All = PLASMAR_Single_Line_List(PLASMAR)
        Lines = []
        for entry in range(len(Y_Float)):
            if Y_Float[entry] == 1:
                Lines.append(Lines_All[entry])    
        All_model = GZ_Pickle(all_model)
        Float_model = GZ_Pickle(float_model)
        All_Predictions = []
        Float_Predictions = []
        All_Data = All_model.predict_proba(Array)
        Float_Data = Float_model.predict_proba(Float_Array)
        All_Predictions.append(All_Data)
        Float_Predictions.append(Float_Data)
        for entry in range(len(Lines)):
            Lines[entry][-1] = float(Lines[entry][-1])
            Lines[entry].append(All_Predictions[0][entry][1])
            Lines[entry].append(Float_Predictions[0][entry][1])
            Best = max(Lines[entry][-3:])
            Lines[entry].append(Best)
        Lines.sort(key=operator.itemgetter(-1), reverse=True)
    f = open(PLASMAR, 'r')
    String1 = f.readline()
    Cats = String1[0:-1]
    f.close()
    Out = open(PLASMAR[0:-4] + '_Probs.txt', 'w')
    Out.write(Cats + '\tAll Prob\tFloat Prob\tMax Prob\n')
    for line in Lines:
        Line = '\t'.join(line[0:-4])
        Out.write(Line + '\t' + str(line[-4]) + '\t' + str(line[-3]) + '\t' + str(line[-2]) + '\t' + str(line[-1]) + '\n')
    Out.close()
    subprocess.call('rm ' + PLASMAR, shell=True)

script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
Upper = os.path.dirname(script_directory)

Comp_Maker_Pool_Folder('PLASMAR_Matches/', 'PSL/', 'Overlap_Comp.txt')
Overlap_Pickle_Scored_Combined('Overlap_Comp.txt', Upper + '/models/PLASMAR_Overlap_All_ET_SKL_1.6.1.pkl.gz', Upper + '/models/PLASMAR_Overlap_Float_ET_SKL_1.6.1.pkl.gz')
