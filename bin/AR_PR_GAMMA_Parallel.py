import sys
import os
import glob
import subprocess
from multiprocessing import Pool

subprocess.call('mkdir AMR', shell=True)
subprocess.call('mkdir PR', shell=True)

def GAMMA_List(input_list):
    Version = input_list[0]
    input_fasta = input_list[1]
    input_DB = input_list[2]
    Folder = input_list[3]
    Name = input_fasta.split('/')[-1]
    subprocess.call(Version + ' ' + Name + ' ' + input_DB + ' ' + Folder + Name[0:-6] + '_' + Folder[0:-1], shell=True)

def GAMMA_Parallel(AR_DB, PF_DB):
    List1 = glob.glob('*.fasta')
    Run_List = []
    for files in List1:
        Run_List.append(['GAMMA.py', files, AR_DB, 'AMR/'])
        Run_List.append(['GAMMA-S.py', files, PF_DB, 'PR/'])
    with Pool() as pool:
        pool.map(GAMMA_List, Run_List)

script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
Upper = os.path.dirname(script_directory)
AR_DB = Upper + '/databases/AR_Database.fasta'
PF_DB = Upper + '/databases/Replicon_Database.fasta'

GAMMA_Parallel(AR_DB, PF_DB)
