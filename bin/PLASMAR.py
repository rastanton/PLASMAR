import os
import subprocess

script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)


subprocess.call('python ' + script_directory + '/AR_PR_GAMMA_Parallel.py', shell=True)
subprocess.call('python ' + script_directory + '/PLASMAR_Matches_Parallel.py', shell=True)
subprocess.call('python ' + script_directory + '/PLASMAR_Overlap_Parallel.py', shell=True)
subprocess.call('python ' + script_directory + '/PLASMAR_Overlap_Report.py', shell=True)
subprocess.call('python ' + script_directory + '/PLASMAR_Presence.py', shell=True)
subprocess.call('python ' + script_directory + '/PLASMAR_Summary_Report.py', shell=True)
