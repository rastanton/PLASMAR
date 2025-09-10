import numpy as np
from scipy.stats import iqr
import glob
import pickle
from sklearn.ensemble import HistGradientBoostingClassifier
import gzip
import sys
import os

def Array_Maker(input_list):
    Out = np.array(input_list)
    return Out

def Gene_Puller(input_file):
    f = open(input_file, 'r')
    String1 = f.readline()
    Genes = []
    for line in f:
        List1 = line.split('\t')
        gene = List1[3]
        if (gene in Genes) == False:
            Genes.append(gene)
    f.close()
    return Genes
    
def Data_Puller_Gene(input_file, gene, start, stop):
    f = open(input_file, 'r')
    String1 = f.readline()
    Data = []
    for entry in range(stop - start):
        Data.append([])
    for line in f:
        List1 = line.split('\t')
        if (List1[3] == gene):
            for position in list(range(start,stop)):
                Data[position - start].append(float(List1[position]))
    f.close()
    Out = []
    for entry in Data:
        Info = Array_Maker(entry)
        Out.append(Info)
    return Out

def IQR_NP(input_array):
    iqr = np.subtract(*np.percentile(input_array, [75, 25]))
    return iqr

def Data_Lister(input_array):
    List = [min(input_array), max(input_array), np.median(input_array), np.mean(input_array), np.std(input_array), IQR_NP(input_array)]
    return List

def PLASMAR_Condenser(input_plasmar):
    Genes = Gene_Puller(input_plasmar)
    All = []
    for entry in Genes:
        Data = Data_Puller_Gene(input_plasmar, entry, 6, 24)
        All.append(Data)
    Lines = []
    for Pos in range(len(All)):
        Data = All[Pos]
        Line = []
        for entry in Data:
            line = Data_Lister(entry)
            Line = Line + line
        Line.append(len(Data[0]))
        Line = [Genes[Pos], Line]
        Lines.append(Line)
    return Lines

def GZ_Pickle(input_pickle):
    with gzip.open(input_pickle, 'rb') as file:
        Model = pickle.load(file)
    return Model

def PLASMAR_Scorer(condensed_data, Model):
    Scores = []
    for entry in condensed_data:
        Array = Array_Maker(entry[1])
        score = Model.predict_proba([Array])
        score = score[0][1]
        Scores.append([entry[0], score])
    return Scores

def PLASMAR_Condenser_glob(Model_Pickle, output_file):
    Model = GZ_Pickle(Model_Pickle)
    List1 = glob.glob('P10/*Prob.txt')
    Out = open(output_file, 'w')
    Out.write('ID\tCarbapenemase\tPlasmid_Prob\n')
    for entry in List1:
        Info = PLASMAR_Condenser(entry)
        Scores = PLASMAR_Scorer(Info, Model)
        Name = entry[4:-17]
        for line in Scores:        
            Out.write(Name + '\t' + str(line[0]) + '\t' + str(line[1]) + '\n')
    Out.close()

script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
Upper = os.path.dirname(script_directory)
        
PLASMAR_Condenser_glob(Upper + '/models/PLASMAR_Presence_HGB_SKL_1.6.1.pkl.gz', 'Plasmid_Presence.txt')
