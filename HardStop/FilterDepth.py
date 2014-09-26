#!/usr/bin/env python
#chr13:30215500-30216000 output/rgn_1    1       26.84
#chr13:37798000-37798500 output/rgn_10   2       19.87

import argparse
ap = argparse.ArgumentParser()
ap.add_argument("depth", help="depth table, 4 fields, region dir nContigs targetDepth")
ap.add_argument("--ncontigs", help="number of contigs.", default=None,type=int)
ap.add_argument("--depth", help="min depth")


