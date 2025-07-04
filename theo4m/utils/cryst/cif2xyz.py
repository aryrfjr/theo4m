# https://wiki.fysik.dtu.dk/ase/ase/spacegroup/spacegroup.html

import sys
# ASE new
from ase.spacegroup import crystal
from ase.visualize import view
import ase.spacegroup as aspg
# CIF library
import CifFile

#
# Subroutines
#

def cifstr2float(cif):
    """Convert CIF string to float, discarding precisions."""
    try:  # scalar
        return float(cif.split("(", 1)[0])
    except AttributeError:  # list
        return [float(n.split("(", 1)[0]) for n in cif]

#
# The script
#

# Reading the .cif file
cf = CifFile.ReadCif(sys.argv[1])
cif_data = cf[cf.keys()[0]]

primitive = sys.argv[2] == "True"

wtd = sys.argv[3] # what to do?

a = cifstr2float(cif_data["_cell_length_a"])
b = cifstr2float(cif_data["_cell_length_b"])
c = cifstr2float(cif_data["_cell_length_c"])

alpha = cifstr2float(cif_data["_cell_angle_alpha"])
beta = cifstr2float(cif_data["_cell_angle_beta"])
gamma = cifstr2float(cif_data["_cell_angle_gamma"])

spg = int(cif_data["_symmetry_Int_Tables_number"])

atoms = cif_data.GetLoop("_atom_site_label")

# Building the cell
symbs = []
coords = []
for atom in atoms:
    if len(atom[0]) == 3 and "B" in atom[0]:
        symbs.append(str(atom[0][0:len(atom[0])-2]))
    elif len(atom[0]) == 3 and "O" in atom[0]:
        symbs.append(str(atom[0][0:len(atom[0])-2]))
    else:
        symbs.append(str(atom[0][0:len(atom[0])-1]))
    coords.append((cifstr2float(atom[4]), cifstr2float(atom[5]), cifstr2float(atom[6])))

cryst = crystal(symbs, # Atoms in the crystal
coords,                # Atomic positions (fractional coordinates)
spacegroup=spg,        # International number of the spacegroup of the crystal
cellpar=[a, b, c, alpha, beta, gamma],  # Unit cell (a, b, c, alpha, beta, gamma) in Angstrom, Degrees
primitive_cell=primitive)

if wtd == "show-all":
    view(cryst)
    if primitive:
        print "primitive cell"
    else:
        print "conventional cell"
    print ""

if wtd == "show-all" or wtd == "show-cell":
    print '%16.8f %16.8f %16.8f' % tuple(cryst.get_cell()[0])
    print '%16.8f %16.8f %16.8f' % tuple(cryst.get_cell()[1])
    print '%16.8f %16.8f %16.8f' % tuple(cryst.get_cell()[2])

if wtd == "show-all" or wtd == "show-nat": print len(cryst.get_positions())

if wtd == "show-all" or wtd == "show-xyz":
    for i in range(len(cryst.get_positions())):
        print '%s %16.8f %16.8f %16.8f' % (cryst.get_chemical_symbols()[i], cryst.get_positions()[i][0], cryst.get_positions()[i][1], cryst.get_positions()[i][2])

if len(sys.argv) > 4:
    px = int(sys.argv[4])
    py = int(sys.argv[5])
    pz = int(sys.argv[6])
    crystp = cryst.repeat((px, py, pz))

    view(crystp)

    # http://www.python-course.eu/python3_formatted_output.php

    print "conventional (propagated)"
    print ""

    print '%16.8f %16.8f %16.8f' % tuple(crystp.get_cell()[0])
    print '%16.8f %16.8f %16.8f' % tuple(crystp.get_cell()[1])
    print '%16.8f %16.8f %16.8f\n' % tuple(crystp.get_cell()[2])
    print len(crystp.get_positions())
    for i in range(len(crystp.get_positions())):
        print '%s %16.8f %16.8f %16.8f' % (crystp.get_chemical_symbols()[i], crystp.get_positions()[i][0], crystp.get_positions()[i][1], crystp.get_positions()[i][2])

