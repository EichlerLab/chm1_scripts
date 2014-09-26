#!/usr/bin/env python

import argparse
import sys

ap = argparse.ArgumentParser(description="using left and right coordinates of split alignments, try and find adjacent alignments")
ap.add_argument("left", help="left bed.")
ap.add_argument("right", help="right bed")
ap.add_argument("--delta", help="How close the sides should be.", default=500)
args = ap.parse_args()

leftFile = open(args.left)
rightFile = open(args.right)

left = leftFile.readlines()
right = rightFile.readlines()
if (len(left) != len(right)):
    print args.left + " must have a 1-1 correspondence with " + args.right
    print "too few lines."
    sys.exit(1)
    
for i in range(len(left)):
    leftv = left[i].split()
    rightv = right[i].split()
    if (leftv[3][0:-2] != rightv[3][0:-2]):
        print args.left + " must have a 1-1 correspondence with " + args.right
        print leftv[3]
        print rightv[3]
        sys.exit(1)
    if (leftv[0] == rightv[0] and abs(int(leftv[2]) - int(rightv[1])) < args.delta):
        print left[i].strip() + " " + right[i].strip()
    
        

    
