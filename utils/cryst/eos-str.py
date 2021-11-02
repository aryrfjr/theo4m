#!/usr/bin/python

import sys

from numpy import array, dot
import sys

ref_symbs = ["Zr","Cu","Al"] # in the order of LAMMPS

fileobj = open(sys.argv[1])
lines = fileobj.readlines()
fileobj.close()

# reading the cell
cell = [x.split()[0:3] for x in lines[0:3]]
cell = array([[float(col) for col in row] for row in cell])
natoms = int(lines[3])
species = [line.split()[0] for line in lines[4:4+natoms]]
geom = array([[float(col) for col in line.split()[1:4]]
    for line in lines[4:4+natoms]])
ic = 1
for i in range(-6, 7, 1):
    m = i / 100.0
    ncell = cell * ((1.0+m)**(1.0/3.0))
    ngeom = dot(geom, ncell)
    # the general files .ncell and .ncoord
    fileobj = open(str(ic)+".ncell","w")
    fileobj.write("%16.8f %16.8f %16.8f\n" % (ncell[0,0],ncell[0,1],ncell[0,2]))
    fileobj.write("%16.8f %16.8f %16.8f\n" % (ncell[1,0],ncell[1,1],ncell[1,2]))
    fileobj.write("%16.8f %16.8f %16.8f\n" % (ncell[2,0],ncell[2,1],ncell[2,2]))
    fileobj.close()
    fileobj = open(str(ic)+".ncoord","w")
    for ia in range(natoms):
        fileobj.write("%s %16.8f %16.8f %16.8f\n" % (species[ia], ngeom[ia][0], ngeom[ia][1], ngeom[ia][2]))
    fileobj.close()
    # the LAMMPS data files
    fileobj = open(str(ic)+"-LAMMPS.ndata","w")
    fileobj.write("# Coordinates\n\n")
    fileobj.write("%d atoms\n" % natoms)
    fileobj.write("3 atom types\n\n")
    fileobj.write("0.0 %16.8f xlo xhi\n" % ncell[0,0])
    fileobj.write("0.0 %16.8f ylo yhi\n" % ncell[1,1])
    fileobj.write("0.0 %16.8f zlo zhi\n\n" % ncell[2,2])
    fileobj.write("Masses\n\n")
    fileobj.write("1 91.2240   # Zr\n")
    fileobj.write("2 63.5463   # Cu\n")
    fileobj.write("3 26.9815   # Al\n\n")
    fileobj.write("Atoms\n\n")
    for ia in range(natoms):
        fileobj.write("%d %d %16.8f %16.8f %16.8f\n" % (ia+1, 
                                                     ref_symbs.index(species[ia])+1, 
                                                     ngeom[ia][0], ngeom[ia][1], ngeom[ia][2]))
    fileobj.close()
    ic+=1

