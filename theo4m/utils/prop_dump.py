#!/usr/bin/python

#
# python ~/dev-python/theo4m/utils/prop_dump.py <DUMP_FILE> <PX> <PY> <PZ> <VX> <VY> <VZ> <OUTPUT_TYPE>
#
# TODO: generalize to other compositions
#


from numpy import dot
import sys
sys.path.append("/home/aryjr/dev-python")
from theo4m.io import read_lammps

# the LAMMPS $RUN_DIR
dfile = sys.argv[1]
# chemical elements in the same order of LAMMPS input file
ref_symbs = ["Zr", "Cu", "Al"]
# propagation
px = int(sys.argv[2])
py = int(sys.argv[3])
pz = int(sys.argv[4])
# vacuum 
vx = int(sys.argv[5])
vy = int(sys.argv[6])
vz = int(sys.argv[7])

output_type = sys.argv[8]

# Reading the last step of the LAMMPS .dump file
rsteps, ratoms, rccell = read_lammps(lmpoutput = dfile, 
                                     spc_symbs = ref_symbs, 
                                     frac = True, item = -2)
atoms = ratoms[0]
ccell = rccell[0]
a = ccell[0][0]
b = ccell[1][1]
c = ccell[2][2]

if output_type == "ext_xyz": # this is just for visualization 
    print "%d" % len(atoms))
    print "Lattice=\"%10f 0.0 0.0 0.0 %10f 0.0 0.0 0.0 %10f\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"" % ((a*px)+vx, (b*py)+vy, (c*pz)+vz))
elif output_type == "data":
    print "# Coordinates\n"
    print str(len(atoms)*px*py*pz)+" atoms"
    print sys.argv[2]+" atom types\n"
    print "   0.0      "+str(a*px)+" xlo xhi"
    print "-100.0      "+str((b*py)+100)+" ylo yhi" # vacuum of 100 angst. around the rod
    print "-100.0      "+str((c*pz)+100)+" zlo zhi\n" # vacuum of 100 angst. around the rod
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
                #print "%d %d %10f %10f %10f" % (iat, ref_symbs.index(symbs[i])+1, geom[i][0]+(a*ix), geom[i][1]+(b*iy), geom[i][2]+(c*iz))
                #iat += 1
                print "%s %10f %10f %10f" % (symbs[i], geom[i][0]+(a*ix), geom[i][1]+(b*iy), geom[i][2]+(c*iz))

