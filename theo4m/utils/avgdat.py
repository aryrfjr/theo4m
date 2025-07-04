#!/usr/bin/python
#

#
# This is a script that averages a XY .dat file based on a dX.
#

# Libraries
import sys

dx = int(sys.argv[1])

fileobj = open(sys.argv[2])
datl = fileobj.readlines()
fileobj.close()

i = 0
while i < len(datl):
    if i + int(dx / 2) > len(datl): break
    x = float(datl[i + int(dx / 2)].split()[0])
    y = 0.0
    for j in range(dx):
        y += float(datl[i + j].split()[1])
    y = y / dx
    i += dx
    print x, y

