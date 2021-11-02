#!/usr/bin/python

# use at CENAPAD-SP like this:
#
# python ~/dev-python/theo4m/utils/dump2extxyz.py ~/UFSCar/MG-NMR/ML/big-data-full/Zr49Cu49Al2/c/md/lammps/100/7 3 Zr Cu Al

from numpy import dot
import sys
sys.path.append("/home/aryjr/dev-python")
from theo4m.io import read_lammps

# the LAMMPS $RUN_DIR
outdir = sys.argv[1]
# chemical elements in the same order of LAMMPS input file
ref_symbs = []
for i in range(3, int(sys.argv[2])+3):
    ref_symbs.append(sys.argv[i])

# 
rsteps, ratoms, rccell = read_lammps(lmpoutput = outdir+"/zca-th300.dump", 
                                     spc_symbs = ref_symbs, 
                                     frac = True, items = [-2])

atoms = ratoms[0]
ccell = rccell[0]

# writting 
ftw = open(outdir+"/2000_extxyzf.xyz","w")
for iats in range(len(ratoms)):
    atoms = ratoms[iats]
    ccell = rccell[iats]
    a = ccell[0][0]
    b = ccell[1][1]
    c = ccell[2][2]
    ftw.write("%d\n" % len(atoms))
    ftw.write("Lattice=\"%10f 0.0 0.0 0.0 %10f 0.0 0.0 0.0 %10f\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n" % (a, b, c))
    symbs = atoms.get_chemical_symbols()
    geom = atoms.get_positions()
    geom = dot(geom, [[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]])
    for i in range(len(atoms)):
        ftw.write("%s %10f %10f %10f\n" % (symbs[i], geom[i][0], geom[i][1], geom[i][2]))
ftw.close()

"""
ftw = open(outdir+"/extxyzf_eos_1999.xyz","w")
for pc in range(-12,13):
    a = ccell[0][0] #+ (ccell[0][0]*float(pc)/float(100))
    b = ccell[1][1] #+ (ccell[1][1]*float(pc)/float(100))
    c = ccell[2][2] + (ccell[2][2]*float(pc)/float(100))
    ftw.write("%d\n" % len(atoms))
    ftw.write("Lattice=\"%10f 0.0 0.0 0.0 %10f 0.0 0.0 0.0 %10f\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n" % (a, b, c))
    symbs = atoms.get_chemical_symbols()
    geom = atoms.get_positions()
    geom = dot(geom, [[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]])
    for i in range(len(atoms)):
        ftw.write("%s %10f %10f %10f\n" % (symbs[i], geom[i][0], geom[i][1], geom[i][2]))
ftw.close()    
"""

