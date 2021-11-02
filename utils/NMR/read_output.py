#!/usr/bin/python

# python ~/dev-python/theo4m/utils/NMR/read_QE-GIPAW.py <SYS_NAME> <NUCLEUS> <MASS_NUMBER>

import sys
sys.path.append("/home/aryjr/dev-python")

from theo4m.nmr.gipaw.io import read_qe
from theo4m.nmr.gipaw.gipawrun import GIPAWRun
from theo4m.nmr.gipaw.gipawatom import GIPAWAtom
#from theo4m.nmr.lapw.io import read_elk
#from theo4m.nmr.lapw.lapwrun import LAPWRun
#from theo4m.nmr.lapw.lapwatom import LAPWAtom
from theo4m.nmr.nuclide import Nuclide

# Quadrupolar moments can be found at (must be multiplied by 100)
# https://www.psi.ch/low-energy-muons/DocumentsEN/nuclear-moments.pdf
# and at
# https://onlinelibrary.wiley.com/doi/pdf/10.1002/cmr.10035

# reading the nuclices.info file
fileobj = open("/home/aryjr/UFSCar/UCLA-KS/nuclides.info")
lines = fileobj.readlines()
fileobj.close()
nuclides = [None]*int(lines[0].split()[0])
i = 0
il = 1
while i < len(nuclides):
    line = lines[il].split()
    if not '#' in line:
        nuclides[i] = Nuclide()
        nuclides[i].set_symbol(line[0])
        nuclides[i].set_mass_number(line[1])
        nuclides[i].set_quadrupolar_moment(float(line[2]))
        nuclides[i].set_nuclear_spin(line[3])
        nuclides[i].set_gyromangnetic_ratio(float(line[4]))
        i += 1
    il += 1

# now the arguments ...
# which code (QE-GIPAW or Elk)
code = sys.argv[1]

# type of input:
# 0 - a single path
# 1 - a file with a set of paths
itype = int(sys.argv[2])
if itype == 0:
    folders = [sys.argv[3]]
else:
    ifile = sys.argv[3]
    fileobj = open(ifile)
    folders = []
    lines = fileobj.readlines()
    for i in range(len(lines)):
        folders.append(lines[i].split()[0])
    fileobj.close()

# prefix of the calculation
prefix = sys.argv[4]

# what to write:
# - 0: run attributes: nkpts, total_energy, Bext
# - 1: per atom efg parameters: Cq and eta
# - 2: per atom hyperfine sigma_s
w2w = int(sys.argv[5])
if w2w in [1]: # tasks that require the atom id
    wa = int(sys.argv[6]) # which atom

# setting the tasks
if w2w == 1:
    tasks = ["efg"]
elif w2w == 2:
    tasks = ["hyperfine", "hyperfine_art"]
else:
    tasks = []

for ifold in range(len(folders)):
    # loading data
    if code == "QE-GIPAW":
        run = read_qe(folders[ifold], prefix, tasks, nuclides, [1])
    elif code == "Elk":
        run = read_elk(folders[ifold], prefix, tasks, nuclides)

    # processing data
    if "efg" in tasks:
        run.compute_efg_parameters()
    elif "hyperfine" in tasks and code == "Elk":
        run.compute_sigma_s()
    elif code == "QE-GIPAW":
        art = "hyperfine_art" in tasks
        nart = "hyperfine" in tasks or "hyperfine_nart" in tasks
        if art or nart: run.compute_sigma_s(art, nart, [1])

    # writing data
    if w2w == 0:
        print "%d %.8f" % (run.get_k_points_number(), run.get_total_energy())#, run.get_Bext())
    elif w2w == 1:
        atom = run.get_atom(wa)
        print "%s %.2f %.2f" % (atom.symbol, atom.get_Cq(), atom.get_eta())
    elif w2w == 2:
        atoms = run.get_atoms_by_symbol(sys.argv[6])
        if code == "QE-GIPAW":
            tp = ""
            for ia in range(len(atoms)):
                ########################
                tp += "%.2f " % (atoms[ia].get_sigma_s_art(0))
                #tp += "%.2f %.2f %.2f %.2f " % (atoms[ia].get_sigma_s_art(0), atoms[ia].get_sigma_s_bare(0), atoms[ia].get_sigma_s_GIPAW_art(0), atoms[ia].get_sigma_s_core_relax_art(0))
                #tp += "%.4f %.2f %.3f %.3f %.3f " % (atoms[ia].get_Bhf()*1000, atoms[ia].get_sigma_s_nart(0), atoms[ia].get_rho_s_bare(0)*1000, atoms[ia].get_rho_s_GIPAW_nart(0)*1000, atoms[ia].get_rho_s_core_relax_nart(0)*1000)
                #tp += "%d %.8f %.2f " % (atoms[ia].tag, atoms[ia].get_rho_s_total_nart(0), atoms[ia].get_sigma_s_nart(0))
                #tp += "%d %.8f %.2f " % (atoms[ia].tag, atoms[ia].get_rho_s_total_nart(0), atoms[ia].get_sigma_s_art(0))
                #tp += "%.3f %.3f %.3f " % (atoms[ia].get_rho_s_bare(0)*1000, atoms[ia].get_rho_s_GIPAW_nart(0)*1000, atoms[ia].get_rho_s_core_relax_nart(0)*1000)
                #tp += "%.8f %.2f " % (atoms[ia].get_rho_s_total_nart(), atoms[ia].get_sigma_s_nart(0))
                #tp += "%.3f %.3f %.3f " % (atoms[ia].get_rho_s_bare(0)*1000000, atoms[ia].get_rho_s_GIPAW_art(0)*1000000, atoms[ia].get_rho_s_core_relax_art(0)*1000000)
            tp += "%.2f %.8f" % (run.get_Bext(), run.get_total_energy())
            print tp
        elif code == "Elk":
            tp = ""
            for ia in range(len(atoms)):
                tp += "%.8f %.2f " % (atoms[ia].get_Bhf(), atoms[ia].get_sigma_s())
            tp += "%.2f" % (run.get_Bext())
            print tp

