# http://nbviewer.ipython.org/github/jochym/qe-doc/blob/master/Primitive_unit_cells.ipynb

#print '###########################################################################'
#print '# do not forget export PYTHONPATH=/opt/spglib-1.7.0/python/ase/lib/python #'
#print '#               export PYTHONPATH=$PYTHONPATH:/home/aryjr/opt/ase-3.11.0  #'
#print '###########################################################################'
#print ''

import sys
#sys.path.append("/home/aryjr/opt/ase-3.11.0")

# ASE system
import ase
from ase import Atom, Atoms
from ase import io
from ase.lattice.spacegroup import crystal
from ase.visualize import view
from numpy import array
from ase.io import write

# Spacegroup/symmetry library
from spglib import spglib

a = float(sys.argv[1])
b = float(sys.argv[2])
c = float(sys.argv[3])
alpha = float(sys.argv[4])
beta = float(sys.argv[5])
gamma = float(sys.argv[6])

spg = int(sys.argv[7])

ntypat = int(sys.argv[8])

symbs = ['None'] * ntypat
coords = [(0.0,0.0,0.0)] * ntypat

nt = 9
for i in range(ntypat):
    symbs[i] = sys.argv[nt]
    coords[i] = (float(sys.argv[nt + 1]),float(sys.argv[nt + 2]),float(sys.argv[nt + 3]))
    nt += 4

cryst = crystal(symbs, # Atoms in the crystal
coords,                # Atomic positions (fractional coordinates)
spacegroup=spg,        # International number of the spacegroup of the crystal
cellpar=[a, b, c, alpha, beta, gamma])  # Unit cell (a, b, c, alpha, beta, gamma) in Angstrom, Degrees

# Write the image to disk file
#ase.io.write('crystal.png',       # The file where the picture get stored
#             cryst,               # The object holding the crystal definition
#             format='png',        # Format of the file
#             show_unit_cell=2,    # Draw the unit cell boundaries
#             rotation='115y,15x', # Rotate the scene by 115deg around Y axis and 15deg around X axis
#             scale=35)            # Scale of the picture

view(cryst)

# http://www.python-course.eu/python3_formatted_output.php

print '%16.8f %16.8f %16.8f' % tuple(cryst.get_cell()[0])
print '%16.8f %16.8f %16.8f' % tuple(cryst.get_cell()[1])
print '%16.8f %16.8f %16.8f\n' % tuple(cryst.get_cell()[2])

for i in range(len(cryst.get_positions())):
    print '%s %16.8f %16.8f %16.8f' % (cryst.get_chemical_symbols()[i], cryst.get_positions()[i][0], cryst.get_positions()[i][1], cryst.get_positions()[i][2])

print 'Space group:', spglib.get_spacegroup(cryst)

