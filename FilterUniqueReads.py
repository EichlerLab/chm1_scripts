#!/usr/bin/env python
import sys

prevRead = ""
for line in sys.stdin.readlines():
    vals = line.split()
    title = vals[0]
    read = '/'.join(title.split('/')[0:2])

    if (read != prevRead):
        sys.stdout.write(line)
    prevRead = read

