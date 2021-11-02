#!/usr/bin/python

#
# use at CENAPAD-SP like this:
#
# python ~/dev-python/theo4m/utils/dump2data.py <DUMP_FILE> 3 Zr Cu Al <PX> <PY> <PZ>
#

from numpy import dot
import sys
sys.path.append("/home/aryjr/dev-python")
from theo4m.io import read_lammps

# the LAMMPS $RUN_DIR
dfile = sys.argv[1]
# chemical elements in the same order of LAMMPS input file
ref_symbs = []
for i in range(3, int(sys.argv[2])+3):
    ref_symbs.append(sys.argv[i])

px = int(sys.argv[i+1])
py = int(sys.argv[i+2])
pz = int(sys.argv[i+3])

# 
rsteps, ratoms, rccell = read_lammps(lmpoutput = dfile, 
                                     spc_symbs = ref_symbs, 
                                     frac = True, items = [-2]	) # only the last item
atoms = ratoms[0]
ccell = rccell[0]

print "# Coordinates\n"
print str(len(atoms)*px*py*pz)+" atoms"
print sys.argv[2]+" atom types\n"

a = ccell[0][0]
b = ccell[1][1]
c = ccell[2][2]

print "  0.0      "+str(a*px)+" xlo xhi"
print "  0.0      "+str(b*py)+" ylo yhi"
print "  0.0      "+str(c*pz)+" zlo zhi\n"
#print "-10.0      "+str((b*py)+10)+" ylo yhi" # vacuum of 10 angst. around the rod
#print "-10.0      "+str((c*pz)+10)+" zlo zhi\n" # vacuum of 10 angst. around the rod
print "Masses\n"
print "1 91.2240   # Zr"
print "2 63.5463   # Cu"
print "3 26.9815   # Al\n"
print "Atoms\n"

symbs = atoms.get_chemical_symbols()
geom = atoms.get_positions()
geom = dot(geom, [[a, 0.0, 0.0], [0.0, b, 0.0], [0.0, 0.0, c]])
iat = 1
for ix in range(px):
    for iy in range(py):
        for iz in range(pz):
            for i in range(len(atoms)):
                print "%d %d %10f %10f %10f" % (iat, ref_symbs.index(symbs[i])+1, geom[i][0]+(a*ix), geom[i][1]+(b*iy), geom[i][2]+(c*iz))
                iat += 1

