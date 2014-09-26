#!/usr/bin/env python

class Alignment:
    def __init__(self):
        self.qname = ""
        self.tname = ""
        self.qstat = 0
        self.qend  = 0
        self.qstrand = 0
        self.qlen  = 0
        self.tstart = 0
        self.tend   = 0
        self.tstrand = 0
        self.tlen  = 0
        self.score = 0
        self.number = 0
        self.identity = 0

    def ToString():
        members = [qname, tname, qstat, qend, qstrand, qlen, tstart, tend, tstrand, tlen, score, number,        identity]
        return [

class M4Reader:
    def __init__(filename):
        self.fh = open(filename)
        self.prev = None

    def GetNext():
        line = self.fh.readline()
        if (line == ""):
            return None
        vals = line.split()
        a = Alignment()
        a.qname = line[0]
        a.tname = line[1]
        a.tstrand = int(line[2])
        a.qstrand = int(line[3])
        a.identity = float(line[4])
        a.score = int(line[5])
        a.tstart = int(line[6])
        a.tend = int(line[7])
        a.tlen = int(line[8])
        a.qstart = int(line[9])
        a.qend = int(line[10])
        a.qlen = int(line[11])
        if (self.prev is not None and self.prev.qname == a.qname):
            a.number = slelf.prev.number + 1
        self.prev = a
        return a

