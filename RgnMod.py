#!/usr/bin/env python

import Tools
import argparse
ap = argparse.ArgumentParser(description="Modify a region.")
ap.add_argument("region", help="UCSC style region")
ap.add_argument("--left", help="Move left by this much.", type=int, default=0)
ap.add_argument("--right", help="Move right by this much", type=int, default=0)

args = ap.parse_args()

region = Tools.ParseRegionStr(args.region)

rgn2 = region[0] + ":" + str(region[1] + args.left) + "-" + str(region[2] + args.right)

print rgn2
