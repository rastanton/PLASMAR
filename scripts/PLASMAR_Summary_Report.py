import sys
import glob
import subprocess

def Asterisk_Remover(input_string):
    Out = ''
    for entry in input_string:
        if entry != '*':
            Out += entry
    return Out

def Top_Plasmids(input_PLASMAR):
    f = open(input_PLASMAR, 'r')
    String1 = f.readline()
    CPs = []
    Plasmids = []
    Scores = []
    Counts = []
    Lengths = []
    ARs = []
    PFs = []
    for line in f:
        List1 = line.split('\t')
        CP = Asterisk_Remover(List1[3])
        Plasmid = List1[1][0:-6]
        Score = round(float(List1[-1]), 2)
        Length = "{:,}".format(int(List1[2]))
        AR = List1[4]
        PF = List1[5]
        if (CP in CPs) == False:
            Plasmids.append(Plasmid)
            CPs.append(CP)
            Scores.append(Score)
            ARs.append(AR)
            PFs.append(PF)
            Lengths.append(Length)
            Counts.append(1)
        else:
            Pos = CPs.index(CP)
            Counts[Pos] += 1
    f.close()
    Out = []
    for entry in range(len(Plasmids)):
        Data = [CPs[entry], str(Counts[entry]), Plasmids[entry], ARs[entry], PFs[entry], Lengths[entry], str(Scores[entry])]
        Out.append(Data)
    return Out

def CP_List_Maker(AMR_gamma):
    """Makes a list of carbapenemases found in a GAMMA file"""
    Genes = []
    f = open(AMR_gamma, 'r')
    String1 = f.readline()
    for line in f:
        Gene = line.split()[0]
        Allele = Gene.split('__')[-1]
        Contig = line.split()[1]
        if ('CARBAPENEM' in Gene) and (Allele in Genes) == False:
            Genes.append(Allele)
    f.close()
    return Genes 
        
def Overlaps_Over_20(input_overlap_file, ID):
    CPs = []
    Matches = []
    f = open(input_overlap_file, 'r')
    String1 = f.readline()
    for line in f:
        List1 = line.split('\t')
        CP = List1[2]
        Match = List1[1]
        Score = float(List1[-1])
        if (List1[0] == ID):
            if Score >= 0.2:
                if (CP in CPs):
                    Pos = CPs.index(CP)
                    Matches[Pos].append(Match)
                else:
                    CPs.append(CP)
                    Matches.append([Match])
    f.close()
    Out = []
    for entry in range(len(CPs)):
        Out.append([CPs[entry], Matches[entry]])
    return Out

def Plasmid_Presence_Probs(input_presence_file, ID):
    CPs = []
    Scores = []
    f = open(input_presence_file, 'r')
    String1 = f.readline()
    for line in f:
        List1 = line.split()
        if ID == List1[0]:
            CP = Asterisk_Remover(List1[1])
            Score = str(round(float(List1[2]), 2))
            CPs.append(CP)
            Scores.append(Score)
    f.close()
    Out = []
    for entry in range(len(CPs)):
        Out.append([CPs[entry], Scores[entry]])
    return Out

def Report_Lines(AMR_GAMMA, PLASMAR, overlap_file, presence_file, ID):
    """Makes a report line for an ID"""
    Lines = []
    CPs = CP_List_Maker(AMR_GAMMA)
    Plasmids = Top_Plasmids(PLASMAR)
    Overlaps = Overlaps_Over_20(overlap_file, ID)
    Presence_Probs = Plasmid_Presence_Probs(presence_file, ID)
    Scored_Matches = []
    for entry in Plasmids:
        List1 = entry
        CP = List1[0]
        Scored_Matches.append(CP)
        True_Overlaps = 'None'
        for entry2 in Presence_Probs:
            if (entry2[0] == CP):
                List1.append(entry2[1])
        for entry3 in Overlaps:
            if (entry[0] == CP):
                True_Overlaps = ','.join(entry3[1])
        List1.append(True_Overlaps)
        Line = '\t'.join(List1)
        Line = ID + '\t' + Line
        Lines.append(Line)
    for entry in CPs:
        if (entry in Scored_Matches) == False:
            Line = ID + '\t' + entry + '\t0\tNo matching plasmids\tN/A\tN/A\tN/A\tNA/A\tN/A\tN\A'
            Lines.append(Line)
    return Lines
            
def Report_Glob(input_folder):
    List1 = glob.glob(input_folder + '/*.fasta')
    Out = open(input_folder + '/PLASMAR_Summary.txt', 'w')
    Out.write('ID\tCarbapenemase\tDatabase Matches\tClosest Plasmid\tClosest Plasmid AR Genes\tClosest Plasmid Replicons\tClosest Plasmid Length (bp)\tClosest Plasmid Score\tCarbapenemase Plasmid Presence Score\tOverlaps > 0.2\n')
    for entry in List1:
        Name = entry.split('/')[-1][0:-6]
        Lines = Report_Lines(input_folder + '/AMR/' + Name + '_AMR.gamma', input_folder + '/P10/' + Name + '_Matches_Prob.txt', input_folder + '/Overlap_Comp_Report.txt',
                             input_folder + '/Plasmid_Presence.txt', Name)
        for line in Lines:
            Out.write(line + '\n')
    Out.close()

Report_Glob('./')
subprocess.call('mkdir Reports', shell=True)
subprocess.call('mv *.txt Reports/', shell=True)
