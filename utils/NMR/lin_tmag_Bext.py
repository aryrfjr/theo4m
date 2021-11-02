#!/usr/bin/python

# python ~/dev-python/examples/lin_tmag_Bext.py 4 0.001 0.003 0.005 0.007
# http://scipy-cookbook.readthedocs.io/items/LinearRegression.html

#
# Subroutines
#

def frange(start, end, jump):
  while start <= end:
    yield start
    start += jump

#
# The script
#

import sys
sys.path.append("/home/aryjr/dev-python")
from numpy import corrcoef
from scipy import polyfit
from aseext.nmr.gipaw.io import read_qe
from aseext.nmr.gipaw.gipawrun import GIPAWRun
from aseext.nmr.gipaw.gipawatom import GIPAWAtom
from aseext.nmr.nuclide import Nuclide

ntmag = int(sys.argv[1])

x = [0.0]*ntmag
y = [0.0]*ntmag

for i in range(ntmag):
    folder = "ss-20-90-720-fd-0.008-false-sp-" + sys.argv[2+i]
    gas = read_qe(folder, "CuAl2", ["hyperfine_art", "hyperfine_nart"], None, [1])
    gas.compute_sigma_s(False, True, [1])
    x[i] = gas.get_total_magnetization()
    y[i] = gas.get_Bext()
    print x[i], y[i]

cc = corrcoef(x, y)[0,1]
(ar,br) = polyfit(x, y, 1)

print ""
print ar, br, cc
print ""

for tmag in frange(0.0000,0.003,0.0005):
    print br + ar * tmag, tmag

print ""

print (float(sys.argv[2+i+1]) - br) / ar, sys.argv[2+i+1]

