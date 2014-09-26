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

    def ToString(self):
        members = ["qname", "tname", "qstart", "qend", "qstrand", "qlen", "tstart", "tend", "tstrand", "tlen", "score", "number", "identity"]

        #return str(self.__dict__.values())
        return " ".join([str(getattr(self,members[i])) for i in range(len(members))])


class M4Reader:
    def __init__(self, filename):
        self.fh = open(filename)
        self.prev = None

    def GetNext(self):
        line = self.fh.readline()
        if (line == ""):
            return None
        vals = line.split()
        a = Alignment()
        a.qname = vals[0]
        a.tname = vals[1]
        a.tstrand = int(vals[2])
        a.qstrand = int(vals[3])
        a.score = int(vals[4])
        a.identity = float(vals[5])
        a.tstart = int(vals[6])
        a.tend = int(vals[7])
        a.tlen = int(vals[8])
        a.qstart = int(vals[9])
        a.qend = int(vals[10])
        a.qlen = int(vals[11])
        if (self.prev is not None and self.prev.qname == a.qname):
            a.number = self.prev.number + 1
        self.prev = a
        return a

