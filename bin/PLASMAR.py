#!/usr/bin/env python3

import sys
import os
import subprocess
import argparse

def PLASMAR_Folder(input_folder):
    script_path = os.path.abspath(__file__)
    script_directory = os.path.dirname(script_path)
    os.chdir(input_folder)
    subprocess.call('python ' + script_directory + '/AR_PR_GAMMA_Parallel.py', shell=True)
    subprocess.call('python ' + script_directory + '/PLASMAR_Matches_Parallel.py', shell=True)
    subprocess.call('python ' + script_directory + '/PLASMAR_Overlap_Parallel.py', shell=True)
    subprocess.call('python ' + script_directory + '/PLASMAR_Overlap_Report.py', shell=True)
    subprocess.call('python ' + script_directory + '/PLASMAR_Presence.py', shell=True)
    subprocess.call('python ' + script_directory + '/PLASMAR_Summary_Report.py', shell=True)
    

def main():
    parser = argparse.ArgumentParser(description="""PLASMAR estimates the probability that multiple short read assemblies carry the same carbapenemase plasmid\nUSAGE: PLASMAR path/to/fasta/folder""")
    parser.add_argument("-f", "--folder", required=True,
                        help="location of the folder with the *.fasta files")
    args = parser.parse_args()
    PLASMAR_Folder(args.folder)

if __name__ == "__main__":
    sys.exit(main())

