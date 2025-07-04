#!/usr/bin/python

import sys
sys.path.append("/home/aryjr/dev-python")

from aseext.nmr.gipaw.io import read_qe
from aseext.nmr.gipaw.gipawrun import GIPAWRun
from aseext.nmr.gipaw.gipawatom import GIPAWAtom
from aseext.nmr.nuclide import Nuclide

nuclides = [None]*2
nuclides[0] = Nuclide()
nuclides[0].set_symbol("Al")
nuclides[0].set_quadrupolar_moment(14.661) # obs, in QE-GIPAW this value must be multiplied by 100
nuclides[0].set_nuclear_spin(5/2)
nuclides[1] = Nuclide()

ifile = sys.argv[1]
fileobj = open(ifile)
prefixes = []
symbols = []
folders = []
lines = fileobj.readlines()
for i in range(len(lines)):
    prefixes.append(lines[i].split()[0])
    symbols.append(lines[i].split()[1])
    folders.append(lines[i].split()[2])
fileobj.close()

tasks = ["hyperfine_art", "hyperfine_nart", "efg"]

for isys in range(len(prefixes)):
    # loading data
    nuclides[1].set_symbol(symbols[isys])
    nuclides[1].set_quadrupolar_moment(0.0)
    nuclides[1].set_nuclear_spin(0.0)
    run = read_qe(folders[isys], prefixes[isys], tasks, nuclides, [1])
    run.compute_efg_parameters()
    # writing data
    print prefixes[isys]
    atoms = run.get_atoms_by_symbol("Al")
    for ia in range(len(atoms)):
        print "%d %12.8f %12.8f %12.8f %6.2f" % (ia, atoms[ia].get_rho_s_bare(0), 
                      atoms[ia].get_rho_s_GIPAW_art(0), atoms[ia].get_rho_s_core_relax_art(0), 
                      abs(atoms[ia].get_Cq()))
    print ""

