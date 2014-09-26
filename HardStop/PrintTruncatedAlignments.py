#!/usr/bin/env python

import sys
import argparse

ap = argparse.ArgumentParser(description="Print alignments that are not full length.")
ap.add_argument("input", help="Input m0 file.")
ap.add_argument("--minTrunc", help="
