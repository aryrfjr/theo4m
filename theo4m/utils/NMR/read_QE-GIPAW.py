#!/usr/bin/python
#
# Quadrupolar moments can be found at (must be multiplied by 100)
# https://www.psi.ch/low-energy-muons/DocumentsEN/nuclear-moments.pdf
# and at
# https://onlinelibrary.wiley.com/doi/pdf/10.1002/cmr.10035
# 
# python ~/dev-python/theo4m/utils/NMR/read_QE-GIPAW.py 0 <INPUT_FILE> <WHAT_TO_SHOW> <NUCLIDE> <N_SITES> <SITES>
# with <N_SITES> <SITES> the number of different sites of the nuclide at issue and their respective ids
#
# or
# 
# python ~/dev-python/theo4m/utils/NMR/read_QE-GIPAW.py 1 <PREFIX> <DIR> <ELEMENT> <WHAT_TO_SHOW> [OPTIONAL: --nuc=<NUCLIDE> --crx=<CORE_RELAX>]
#

import sys
sys.path.append("/home/aryjr/dev-python")

from theo4m.nmr.gipaw.io import read_qe
from theo4m.nmr.gipaw.gipawrun import GIPAWRun
from theo4m.nmr.gipaw.gipawatom import GIPAWAtom
from theo4m.nmr.nuclide import Nuclide

# reading the nuclides.info file
fileobj = open("/home/aryjr/dev-python/theo4m/utils/NMR/nuclides.info")
lines = fileobj.readlines()
fileobj.close()
nuclides = {}
i = 0
il = 1
while i < int(lines[0].split()[0]):
    line = lines[il].split()
    if not '#' in line:
        nuclide = Nuclide()
        nuclide.set_symbol(line[0])
        nuclide.set_mass_number(line[1])
        nuclide.set_quadrupolar_moment(float(line[2]))
        nuclide.set_nuclear_spin(float(line[3]))
        nuclide.set_gyromangnetic_ratio(float(line[4]))
        nuclides[str(line[1])+line[0]] = nuclide
        i += 1
    il += 1

wtd = int(sys.argv[1])
# read a .nmr input file
# python ~/dev-python/theo4m/utils/NMR/read_QE-GIPAW.py 0 <INPUT_FILE> <WHAT_TO_SHOW> <NUCLIDE> <N_SITES> <SITES>
# with <N_SITES> <SITES> the number of different sites of the nuclide at issue and their respective ids
if wtd == 0:
    # Reading command line arguments
    wts = sys.argv[3]
    nuclide_label = sys.argv[4]
    n_sites = int(sys.argv[5]) # the set of sites to show information of
    sites = [-1]*n_sites
    for i in range(n_sites):
        sites[i] = int(sys.argv[i + 6])
    # Reading the input file
    fileobj = open(sys.argv[2])
    lines = fileobj.readlines()
    fileobj.close()
    # the number of simulations
    nsim = int(lines[0].split()[0])
    # information about the computational reference
    ref_prefix = lines[1].split()[0]
    print "Reference: " + ref_prefix
    ref_site = int(lines[1].split()[1])
    ref_exp_shift = float(lines[1].split()[2])
    ref_sims = [None]*nsim
    isim = 0
    for il in range(2, 2+nsim):
        print "Loading " + lines[il].split()[0]
        ref_sims[isim] = read_qe(lines[il].split()[0], ref_prefix, ["nmr"])
        isim += 1
    ll = il + 1
    # now loading the information about the system
    sys_prefix = lines[ll].split()[0]
    sims_labels = [""]*nsim
    print "System: " + sys_prefix
    sys_sims_ss = [None]*nsim
    sys_sims_so = [None]*nsim
    isim = 0
    for il in range(ll + 1, ll + nsim * 3, 3):
        sims_labels[isim] = lines[il].split()[0]
        print "System label: " + sims_labels[isim]
        print "Loading: " + lines[il+1].split()[0]
        sys_sims_ss[isim] = read_qe(lines[il+1].split()[0], sys_prefix, ["hyperfine_art","hyperfine_nart"], nuclide_label=nuclide_label)
        sys_sims_ss[isim].compute_sigma_s(True, True)
        print "Loading: " + lines[il+2].split()[0]
        sys_sims_so[isim] = read_qe(lines[il+2].split()[0], sys_prefix, ["nmr", "efg"], nuclide_label=nuclide_label)
        sys_sims_so[isim].compute_efg_parameters(nuclides[nuclide_label])
        isim += 1
    print ""
    
    # Now the output
    if wts == "kiso":
        for i in range(n_sites):
            ss_site = sys_sims_ss[0].get_atom(sites[i])
            print "K_iso:", ss_site.symbol+"("+str(sites[i])+")"
            print "%15s %15s %15s" % ("simulation", "nart", "art")
            for ism in range(nsim):
                sref = ref_sims[ism].get_atom(ref_site).get_sigma_o()
                ssnart = sys_sims_ss[ism].get_atom(sites[i]).get_sigma_s_nart()
                ssart = sys_sims_ss[ism].get_atom(sites[i]).get_sigma_s_art()
                so = sys_sims_so[ism].get_atom(sites[i]).get_sigma_o()
                Kison = (sref + ref_exp_shift) - (ssnart + so)
                Kisoa = (sref + ref_exp_shift) - (ssart + so)
                print "%15s %15.2f %15.2f" % (sims_labels[ism], Kison, Kisoa)
            print ""
    elif wts == "efg":
        for i in range(n_sites):
            ss_site = sys_sims_ss[0].get_atom(sites[i])
            print "EFG:", ss_site.symbol+"("+str(sites[i])+")"
            print "%15s %15s %15s %15s" % ("simulation", "Cq", "eta", "freq.")
            for ism in range(nsim):
                at = sys_sims_so[ism].get_atom(sites[i])
                print "%15s %15.2f %15.2f %15.2f" % (sims_labels[ism], at.get_Cq()*1000, at.get_eta(), at.get_nuQ())
            print ""
# read a directory containing a HF task simulation
# python ~/dev-python/theo4m/utils/NMR/read_QE-GIPAW.py 1 <PREFIX> <DIR> <ELEMENT> [OPTIONAL: --nuc=<NUCLIDE> --crx=<CORE_RELAX>]
elif wtd == 1:
    prefix = sys.argv[2]
    ssdir = sys.argv[3]
    element = sys.argv[4]
    iarg = 5
    nuclide_label = ""
    core_relax = ""
    if len(sys.argv) > iarg and sys.argv[iarg][:6] == "--nuc=":
        nuclide_label = sys.argv[iarg][6:]
        iarg += 1
    if len(sys.argv) > iarg and sys.argv[iarg][:6] == "--crx=":
        core_relax = sys.argv[iarg][6:]
    hfi = read_qe(ssdir, prefix, ["hyperfine_art","hyperfine_nart"], nuclide_label=nuclide_label, core_relax=core_relax)
    hfi.compute_sigma_s(True, True)
    atoms = hfi.get_atoms_by_symbol(element)
    #tp = ""
    tp = "%.2f " % (hfi.get_Bext())
    for ia in range(len(atoms)):
        tp += "%.2f " % (atoms[ia].get_Bhf())
        #tp += "%.2f " % (atoms[ia].get_sigma_s_art())
        #tp += "%.2f " % (atoms[ia].get_sigma_s_bare())
        #tp += "%.2f " % (atoms[ia].get_sigma_s_GIPAW_art())
        #tp += "%.2f " % (atoms[ia].get_sigma_s_core_relax_art())
        #tp += "%.2f " % (atoms[ia].get_sigma_s_nart())
        #tp += "%.8f " % (atoms[ia].get_rho_s_total_art())
        #tp += "%.8f \n" % (atoms[ia].get_rho_s_total_nart())
    #tp += "%.2f %.8f" % (hfi.get_Bext(), hfi.get_total_energy())
    print tp
    #print hfi.get_Bext()

