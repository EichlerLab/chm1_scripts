#!/usr/bin/env python


import argparse
import sys
import subprocess
import tempfile

ap = argparse.ArgumentParser(description="Look for gaps in quiver consensi.")
ap.add_argument("fasta", help="Quiver consensus fasta")
ap.add_argument("genome", help="Reference genome.")
ap.add_argument("tmpdir", help="Work in this dorectory.")
ap.add_argument("--chunk", help="Split consensus into this chunk size", type=int, default=50000)
ap.add_argument("--exe", help="Path to blasr executables.", default="/net/eichler/vol5/home/mchaisso/software/blasr_1/cpp")
args = ap.parse_args()

splitFileName = tempfile.mktemp(suffix=".split.fasta", dir=args.tmpdir)
splitCmd = "{}/sequtils/bin/tileGenome {} {} {} -stride {}".format(args.exe, args.fasta, splitFileName, args.chunk, args.chunk)
subprocess.call(splitCmd.split())

samFileName = tempfile.mktemp(suffix=".sam", dir=args.tmpdir)
alnCmd = "{}/alignment/bin/blasr {} {} -sa {}.sa -sam -sdpTupleSize 14 -minMatch 30 -maxMatch 50 -out {} -nproc 12 -affineOpen 20 -affineAlign -affineOpen 10 -affineExtend 0 -insertion 10 -deletion 10  -bestn 2  ".format(args.exe, splitFileName, args.genome, args.genome, samFileName)

subprocess.call(alnCmd.split())


gapCommand = "/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintGaps.py {} {} --outFile {}.gaps".format(args.genome, samFileName, args.fasta)


subprocess.call(gapCommand.split())
