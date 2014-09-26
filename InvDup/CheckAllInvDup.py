#!/usr/bin/env python



import os
import sys
import subprocess

dirs = os.listdir(".")
command = "/bin/rm summary.txt"
subprocess.call(command.split())
for d in dirs:
    assemblyName = d + "/region/region.ctg.consensus.fasta"
    if (os.path.isdir(d) and os.path.exists(assemblyName)):
        command = "/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/InvDup/CheckInvDup.py {} --regionFile {} --summary summary.txt".format(assemblyName ,d + "/region.txt")
        subprocess.call(command.split())
