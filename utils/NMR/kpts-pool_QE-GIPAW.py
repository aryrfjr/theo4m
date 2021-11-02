#!/usr/bin/python
#

import sys

# reading the input file
fileobj = open(sys.argv[1])
lines = fileobj.readlines()
fileobj.close()

max_cores = int(sys.argv[2])

for line in lines:
    print line
    kpts = float(line.split()[1])
    for i in range(1, 1001):
        div = kpts / i
        if (div).is_integer():
            npool = int(kpts / div)
            if npool * 16 <= max_cores and (float(npool) * 16.0) % 24.0 == 0.0:
                print "div: %d, npool: %d, cores: %d, nodes: %f" % (int(div), npool, npool * 16, float(npool) * 16.0 / 24.0)
    print ""

