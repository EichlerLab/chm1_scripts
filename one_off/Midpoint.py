#!/usr/bin/env python

import Tools
import sys
region = Tools.ParseRegionStr(sys.argv[1])
print region[0] + ":" + str((region[2] + region[1])/2) + "-" + str((region[2] + region[1])/2 + 1)
