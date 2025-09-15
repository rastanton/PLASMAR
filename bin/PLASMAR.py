#!/usr/bin/env python3

import sys
import os
import subprocess
import argparse
import glob

def PLASMAR_Folder(input_folder, matches_only, run_matches):
    script_path = os.path.abspath(__file__)
    script_directory = os.path.dirname(script_path)
    os.chdir(input_folder)
    List1 = glob.glob('*.fasta')
    if len(List1) == 0:
        print('No .fasta files detected')
    elif len(List1) == 1 or matches_only == 1:
        subprocess.call(script_directory + '/AR_PR_GAMMA_Parallel.py', shell=True)
        subprocess.call(script_directory + '/PLASMAR_Matches_Parallel.py ' + str(run_matches), shell=True)
        subprocess.call(script_directory + '/PLASMAR_Presence.py', shell=True)
    else:
        subprocess.call(script_directory + '/AR_PR_GAMMA_Parallel.py', shell=True)
        subprocess.call(script_directory + '/PLASMAR_Matches_Parallel.py ' + str(run_matches), shell=True)
        subprocess.call(script_directory + '/PLASMAR_Overlap_Parallel.py', shell=True)
        subprocess.call(script_directory + '/PLASMAR_Overlap_Report.py', shell=True)
        subprocess.call(script_directory + '/PLASMAR_Presence.py', shell=True)
        subprocess.call(script_directory + '/PLASMAR_Summary_Report.py', shell=True)
        subprocess.call(script_directory + '/PLASMAR_Heatmap.py', shell=True)
    

def main():
    parser = argparse.ArgumentParser(description="""PLASMAR estimates the probability that multiple short read assemblies carry the same carbapenemase plasmid\nUSAGE: PLASMAR path/to/fasta/folder""")
    parser.add_argument("-f", "--folder", required=True,
                        help="location of the folder with the *.fasta files (required)")
    parser.add_argument("-m", "--matches", action="store_true",
                        help="only determines match probabilities to known plasmids and not potential overlaps between multiple sequences")
    parser.add_argument("-a", "--all", action="store_true",
                        help="considers all plasmids with the same carbapenemase allele as potential matches, only suggested for rare carbapenemases since the models were not trained on this data")
    args = parser.parse_args()
    run_matches = 0
    matches_only = 0
    if args.all:
        run_matches = 1
    if args.matches:
        matches_only = 1
    PLASMAR_Folder(args.folder, matches_only, run_matches)
    subprocess.call('rm -rf PSL/', shell=True)

if __name__ == "__main__":
    sys.exit(main())

