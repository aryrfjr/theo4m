#!/usr/bin/python

import sys
#sys.path.append("/home/aryjr/opt/ase-3.11.0")
sys.path.append("/home/aryjr/dev-python")

from ase import Atoms
from ase import Atom
from ase.visualize import view
from ase import io

inpfile = 'cell_coords.xyz'
fileobj = open(inpfile)
lines = fileobj.readlines()
fileobj.close()

cell = []
symbols = []
positions = []

ntypat = int(lines[0].split()[0])

ref_symbs = ['None'] * ntypat
masses = [0.0] * ntypat

for i in range(1, ntypat+1):
    ref_symbs[i-1] = lines[i].split()[0]
    masses[i-1] = float(lines[i].split()[1])

for i in range(ntypat+1,ntypat+4):
    line = lines[i].split()
    cell.append([float(line[0]), float(line[1]), float(line[2])])

nat = int(lines[ntypat+4].split()[0])

for i in range(ntypat+5,ntypat+5+nat):
    line = lines[i].split()
    symbols.append(line[0])
    positions.append([float(line[1]), float(line[2]), float(line[3])])

cryst = Atoms(cell=cell, positions=positions)
cryst.set_chemical_symbols(symbols)

px = int(sys.argv[1])
py = int(sys.argv[2])
pz = int(sys.argv[3])

pcell = cryst.repeat((px, py, pz))

print '# Coordinates'
print ''
print '%d atoms' % (nat * px * py * pz)
print '%d atom types' % ntypat
print ''

print '0.0 %16.8f xlo xhi' % pcell.get_cell()[0][0]
print '0.0 %16.8f ylo yhi' % pcell.get_cell()[1][1]
print '0.0 %16.8f zlo zhi' % pcell.get_cell()[2][2]
print ''
print 'Masses'
print ''
for i in range(ntypat):
    print '%d %16.8f   # %s' % (i+1, masses[i], ref_symbs[i])
print ''
print 'Atoms'
print ''

for i in range(len(pcell.get_positions())):
    print '%d %d %16.8f %16.8f %16.8f' % (i+1, ref_symbs.index(pcell.get_chemical_symbols()[i])+1, pcell.get_positions()[i][0], pcell.get_positions()[i][1], pcell.get_positions()[i][2])

