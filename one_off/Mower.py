#!/usr/bin/env python

import sys
import subprocess
import re
import time

cpuRe = re.compile("cpu=(\d+)\:(\d+)\:(\d+)")


killedFile = open("ended_jobids.txt", "a+")

while (True):
    command = "/net/eichler/vol5/home/mchaisso/scripts/Runtimes.sh cmd"
    proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    (out, err) = proc.communicate()
    for line in out.splitlines():
        vals = line.split()
        if (len(vals) < 2):
            continue
        m = cpuRe.match(vals[1])
        jobid = vals[0]
        if (m is not None and len(m.groups()) == 3):
            g = m.groups()
            h = int(g[0])
            m = int(g[1])
            s = int(g[2])
            if (h >= 2):
                command = "qdel " + jobid
                print command
                subprocess.call(command.split())
                killedFile.write(jobid + "\n")

    time.sleep(20)

