"""

Miscellaneous database utils for ML models.

"""

import gc
import random
import os
import numpy as np
#import quippy
import copy
from ase import Atom
from ase import Atoms
from ase.cell import Cell
from theo4m.atom.desc.soap import soaplist
from theo4m.utils.supercell import get_clusters_from_cell
from theo4m.utils.supercell import get_item_ext_xyz

##########################################################################
#
# Subroutines
#
##########################################################################

#############################################################################################
#
# Originally devised for the EXAFS-guided MC geometry 
# optimization. See more on 21/07/2020.
#
#############################################################################################

def get_dist(atomA, atomB, cell):
    # the cell is cubic, so very simple PBC implementation
    dist_v = [0.0]*3
    for i in range(0,3):
        dist_v[i] = abs(atomA.position[i] - atomB.position[i])
        if dist_v[i] > cell[i][i]/2: dist_v[i] -= cell[i][i]
    d = np.linalg.norm(dist_v)
    return d

#############################################################################################
#
# Originally devised for the EXAFS-guided MC geometry 
# optimization. See more on 24/04/2020-(1).
#
#############################################################################################

def move_atom(iat, atoms, dd):
    mx = random.uniform(-dd, dd)
    my = random.uniform(-dd, dd)
    mz = random.uniform(-dd, dd)
    g = atoms.get_positions()
    g[iat][0] += mx
    g[iat][1] += my
    g[iat][2] += mz
    atoms.set_positions(g)
    #print iat, mx, my, mz
    return mx, my, mz

def unmove_atom(iat, atoms, ux, uy, uz):
    g = atoms.get_positions()
    g[iat][0] -= ux
    g[iat][1] -= uy
    g[iat][2] -= uz
    atoms.set_positions(g)

def swap_atoms(iata, iatb, atoms):
    #print "swapping", iata, iatb
    s = atoms.get_chemical_symbols()
    tmp = s[iata]
    s[iata] = s[iatb]
    s[iatb] = tmp
    atoms.set_chemical_symbols(s)
    return atoms

def __DEPRECATED_topological_index__(atoms, tbonds, tcoords, spc_symbs, short_bonds):
    ind_bonds = 0.0
    ind_coords = 0.0
    symbs = atoms.get_chemical_symbols()
    short_bonded = []
    for i in range(len(atoms)):
        tcc = {}
        for ispc in range(len(spc_symbs)):
            tcc[symbs[i]+"-"+spc_symbs[ispc]] = 0
        cl = atoms.get_1st_coord(3.75, i)
        #print i, symbs[i], cl.get_number_of_neighbors()
        dist_diffs = 0.0
        for inb in range(cl.get_number_of_neighbors()):
            bt = symbs[i]+"-"+cl.get_neighbor_symbol(inb)
            diff = abs(tbonds[bt]-cl.get_distance(inb))
            if diff > short_bonds:
                short_bonded.append(i)
                #print "********* HERE!!! *********"
            #print bt, tbonds[bt], cl.get_distance(inb), diff
            dist_diffs += diff
            tcc[bt] += 1
        #ind_bonds += dist_diffs / float(cl.get_number_of_neighbors()) # using average previously
        ind_bonds += dist_diffs
        dist_diffs = 0.0
        for key, value in tcc.items():
            diff = abs(tcoords[key] - value)
            #print key, value, tcoords[key], diff
            dist_diffs += diff
        #ind_coords += dist_diffs / float(len(tcc)) # using average previously
        ind_coords += dist_diffs
        #print i, ind_bonds, ind_coords
        #print ""
        #print "----------------------------------------------------"
        #cl.to_xyz()
        #print "----------------------------------------------------"
        #print ""
    return ind_bonds, ind_coords, short_bonded

# new version implemented on 18/06/2020
def topological_index(atoms, tbonds, tcoords, spc_symbs, short_bond_thr):
    ind_bonds = 0.0 # index for bonds
    ind_coords = 0.0 # index for coordination numbers
    short_bonded = [] # too short bonds
    # building the local environments in each cell
    # all the clusters will be set at the center of a 
    # 20.0 Angst. cubic cell with no PBC (like a molecule)
    fcell = np.array([[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]])
    clusters = get_clusters_from_cell(atoms, fcell, centered = True)
    ###############################################################################
    # TEST_TAG: tests_to_visualize_local_environments
    #ftest = open("test.xyz","w")
    #ftest.write(get_item_ext_xyz(clusters[44], {}, False))
    #ftest.close()
    ###############################################################################
    for i in range(len(clusters)):
        cluster = clusters[i]
        tcc = {} # types of bonds to compute coordination numbers
        for ispc in range(len(spc_symbs)): # setting up a dictionary
            tcc[cluster[0].symbol+"-"+spc_symbs[ispc]] = 0
        dist_diffs = 0.0 # absolute differences in bond distances
        # looping over the neighbors of the central atoms
        for inb in range(1, len(cluster)):
            # checking the absolute differences in bond distances
            bt = cluster[0].symbol+"-"+cluster[inb].symbol # the bond type
            d = get_dist(cluster[0], cluster[inb], atoms.get_cell())
            # checking if the distance is less than 
            # the short bond threshold
            if d < short_bond_thr:
                short_bonded.append(i)
            # increasing the absolute differences in bond distances
            diff = abs(tbonds[bt]-d)
            dist_diffs += diff
            # increasing the number of bonds of the current type
            tcc[bt] += 1
        # increasing the total absolute differences in bond distances
        ind_bonds += dist_diffs
        # the absolute differences in coordination numbers
        coord_diffs = 0.0
        for key, value in tcc.items():
            diff = abs(tcoords[key] - value)
            coord_diffs += diff
        ind_coords += coord_diffs
    return ind_bonds, ind_coords, short_bonded

def thermodynamic_index(atoms):
    thind = 0.0
    return thind

#############################################################################################
#
# Originally devised for the comparison between the 
# databases generated with EXAFS-guided MC geometry 
# optimization and snapshots from QfM CMD simulations. 
# See more on 11/06/2020.
#
#############################################################################################

# Writing the extended .xyz files to be used by ASAP
# https://github.com/BingqingCheng/ASAP
def write_ASAP_ext_xyz(dir2w, prefix, atoms_set, params_set, lmax, nmax, pcs_SOAPs = True):
    #
    # First of writing an extended .xyz file with all cells 
    # (all ase.Atoms objects) with their respective per-cell 
    # average SOAP vectors.
    soapsl = soaplist.SOAPList(atoms_set, verb = False)
    cutOff = 3.75
    print ("Calculating average SOAPs for %d cells ..." % len(atoms_set))
    soapsl.compute_average(cutoff = str(cutOff), 
                  l_max = lmax, 
                  n_max = nmax, 
                  n_Z = "3", 
                  Z = "{13 29 40}", 
                  n_species = "3", 
                  species_Z = "{13 29 40}")
    print ("Done!")
    ###############################################################################
    # TEST_TAG: tests_to_compare_per-atom_SOAPs
    #soapsl.compute_per_atom(cutoff = str(cutOff), 
    #              l_max = lmax, 
    #              n_max = nmax, 
    #              n_Z = "3", 
    #              Z = "{13 29 40}", 
    #              n_species = "3", 
    #              species_Z = "{13 29 40}")
    ###############################################################################
    file_name = dir2w+"/"+prefix+"_per-cell_avg_SOAPs.xyz"
    ftw = open(file_name,"a")
    ###############################################################################
    # TEST_TAG: tests_to_evaluate_ln_max_convergence
    #tmpsoaps = []
    ###############################################################################
    for i in range(len(atoms_set)):
        asoap = soapsl.get_average_soap(0,i)
        # adding SOAP-related extra params to each cell
        params_set[i]["SOAP_args"] = soapsl.get_asoaps_args(0)
        params_set[i]["SOAPs"] = " ".join(map(str, asoap))
        ###############################################################################
        # TEST_TAG: tests_to_evaluate_ln_max_convergence
        #tmpsoaps.append(asoap)
        #if i > 0:
        #    print (np.dot(tmpsoaps[i-1], tmpsoaps[i]))
        ###############################################################################
        ftw.write(get_item_ext_xyz(atoms_set[i], params_set[i], False))
        ###############################################################################
        # TEST_TAG: tests_to_compare_per-atom_SOAPs
        #if i == 1: # which frame
        #    pasoaps, cats = soapsl.get_per_atom_soaps(0,i)
        #    for isp in range(len(pasoaps)): # per-atom SOAPs of the first frame
        #        if isp == 10: # a chosen per-atom SOAP
        #            ftest = open(dir2w+"/check_per_atom_SOAPs","w")
        #            ftest.write(str(i)+" "+str(cats[isp]-1)+"\n")
        #            for ii in range(len(pasoaps[isp])):
        #                ftest.write(str(pasoaps[isp][ii])+"\n")
        #            ftest.close()
        ###############################################################################
    ftw.close()
    # releasing memory
    del soapsl
    del ftw
    print ("File \"%s\" written." % file_name)
    if pcs_SOAPs:
        #
        #
        # Now, writing the per-specie extended .xyz files with 
        # all corresponding chemical environments from all cells 
        # separated in each frame; with their respective 
        # per-central-atom SOAP vectors.
        # create a set of clusters for each cell (or frame)
        for i in range(len(atoms_set)): # looping over all frames (cells)
            print ("Creating clusters for cell %d ..." % i)
            # all the clusters will be set at the center of a 
            # 20.0 Angst. cubic cell with no PBC (like a molecule)
            fcell = np.array([[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0]])
            clusters = get_clusters_from_cell(atoms_set[i], fcell, centered = True)
            clusters_params = []
            # removing the average SOAPs key from the parameters
            del params_set[i]["SOAPs"]
            # now looping over all atoms in the current frame (atoms_set[i])
            for j in range(len(atoms_set[i])):
                cluster = clusters[j]
                # setting the per-cluster parameter for the extended .xyz file
                parcl = copy.deepcopy(params_set[i]) # TODO: redundant information
                parcl["CentralAtom"] = j
                # getting the per-central-atom SOAP vectors of this cluster
                soapsl = soaplist.SOAPList([cluster], verb = False)
                soapsl.compute_per_atom(cutoff = str(cutOff), 
                          l_max = lmax, 
                          n_max = nmax, 
                          n_Z = "3", 
                          Z = "{13 29 40}", 
                          n_species = "3", 
                          species_Z = "{13 29 40}")
                pasoaps, cats = soapsl.get_per_atom_soaps(0,0)
                parcl["CentralAtom_SOAP"] = " ".join(map(str, pasoaps[0]))
                clusters_params.append(parcl)
                ###############################################################################
                # TEST_TAG: tests_to_visualize_local_environments
                #if i == 0 and j == 55:
                #    ftest = open(dir2w+"/test.xyz","w")
                #    ftest.write(get_item_ext_xyz(cluster, parcl, False))
                #    ftest.close()
                ###############################################################################
                # releasing memory
                del soapsl
            # appending/writing the extended .xyz file
            # with the per-central-atom SOAP vectors.
            file_nameZr = dir2w+"/"+prefix+"_per-central_Zr_atom_SOAPs.xyz"
            file_nameCu = dir2w+"/"+prefix+"_per-central_Cu_atom_SOAPs.xyz"
            file_nameAl = dir2w+"/"+prefix+"_per-central_Al_atom_SOAPs.xyz"
            print ("Writing %d clusters of cell %d in the per-species files:\n - \"%s\"\n - \"%s\"\n - \"%s\" ..." 
                                           % (len(clusters), i, file_nameZr, file_nameCu, file_nameAl))
            ftwZr = open(file_nameZr,"a")
            ftwCu = open(file_nameCu,"a")
            ftwAl = open(file_nameAl,"a")
            for j in range(len(clusters)):
                specie = clusters[j].get_chemical_symbols()[0]
                if specie == "Zr": 
                    ftwZr.write(get_item_ext_xyz(clusters[j], clusters_params[j], False))
                elif specie == "Cu":
                    ftwCu.write(get_item_ext_xyz(clusters[j], clusters_params[j], False))
                elif specie == "Al":
                    ftwAl.write(get_item_ext_xyz(clusters[j], clusters_params[j], False))
                ###############################################################################
                # TEST_TAG: tests_to_compare_per-atom_SOAPs
                #if clusters_params[j]["DumpedItem"] == "2" and clusters_params[j]["CentralAtom"] == 0:
                #    ftest = open(dir2w+"/check_per_atom_SOAPs_xyz","w")
                #    ftest.write(str(int(clusters_params[j]["DumpedItem"])-1)+" "+str(clusters_params[j]["CentralAtom"])+"\n")
                #    for ii in range(len(pasoaps[0])):
                #        ftest.write(str(pasoaps[0][ii])+"\n")
                #    ftest.close()
                ###############################################################################
            print ("Done!")
            ftwZr.close()
            ftwCu.close()
            ftwAl.close()
            # releasing memory
            del clusters
            del clusters_params
            del ftwZr
            del ftwCu
            del ftwAl
            # https://pymotw.com/2/gc/
            #gc.collect()

#############################################################################################
#
# Originally devised for the test planned on 20/06/2020-(2).
#
#############################################################################################

# Writing the extended .xyz files to be used by ASAP
# https://github.com/BingqingCheng/ASAP
def write_QE_pwi(dir2w, optid, atoms, kpts, ecut, calculation):
    nat = len(atoms)
    cell = atoms.get_cell()
    symbs = atoms.get_chemical_symbols()
    file_name = dir2w+"/opt-"+optid+".pw.inp"
    ftw = open(file_name, "w")
    ftw.write(" &CONTROL\n")
    ftw.write("                       title = \'opt-"+optid+"\',\n")
    ftw.write("                 calculation = \'"+calculation+"\',\n")
    ftw.write("                restart_mode = \'from_scratch\',\n")
    ftw.write("                      outdir = \'/scratch/rmnvm2/ary.junior2/work/DBFS/OPT/opt-"+optid+"\',\n")
    ftw.write("                      wfcdir = \'/scratch/rmnvm2/ary.junior2/work/DBFS/OPT/opt-"+optid+"\',\n")
    ftw.write("                  pseudo_dir = \'/scratch/rmnvm2/ary.junior2/UPFs/\',\n")
    ftw.write("                      prefix = \'opt-"+optid+"\',\n")
    ftw.write("                   verbosity = \'high\',\n")
    ftw.write("                     tprnfor = .true.,\n")
    ftw.write("                     tstress = .true.,\n")
    ftw.write("                       nstep = 10,\n")
    ftw.write("               etot_conv_thr = 0.1,\n")
    ftw.write("               forc_conv_thr = 1.0\n")
    ftw.write(" /\n")
    ftw.write(" &SYSTEM\n")
    ftw.write("                       ibrav = 0,\n")
    ftw.write("                   celldm(1) = 1.8897265,\n")
    ftw.write("                         nat = "+str(nat)+",\n")
    ftw.write("                        ntyp = 3,\n")
    ftw.write("                     ecutwfc = "+str(ecut)+",\n")
    ftw.write("                     ecutrho = "+str(ecut*8)+",\n")
    ftw.write("                 occupations = \'smearing\',\n")
    ftw.write("                    smearing = \'fd\',\n")
    ftw.write("                     degauss = 0.008,\n")
    ftw.write("                       nosym = .true.,\n")
    ftw.write("                       noinv = .true.\n")
    ftw.write(" /\n")
    ftw.write(" &ELECTRONS\n")
    ftw.write("            electron_maxstep = 100,\n")
    ftw.write("                    conv_thr = 1.0D-6,\n")
    ftw.write("                 mixing_beta = 0.4\n")
    ftw.write(" /\n")
    ftw.write(" &IONS\n")
    ftw.write(" /\n")
    ftw.write("CELL_PARAMETERS alat\n")
    ftw.write("%10f %10f %10f\n" % (cell[0,0], cell[0,1], cell[0,2]))
    ftw.write("%10f %10f %10f\n" % (cell[1,0], cell[1,1], cell[1,2]))
    ftw.write("%10f %10f %10f\n" % (cell[2,0], cell[2,1], cell[2,2]))
    ftw.write("ATOMIC_SPECIES\n")
    ftw.write("   Al  26.9815   Al.aryjr-1.0.0-ld1.UPF\n")
    ftw.write("   Cu  63.5463   Cu.aryjr-1.0.0-ld1.UPF\n")
    ftw.write("   Zr  91.2240   Zr.aryjr-1.0.0-ld1.UPF\n")
    ftw.write("ATOMIC_POSITIONS\n")
    geom = atoms.get_positions()
    for i in range(nat):
        ftw.write("%s %10f %10f %10f\n" % (symbs[i], geom[i][0], geom[i][1], geom[i][2]))
    ftw.write("K_POINTS automatic\n")
    ftw.write("  "+kpts+" "+kpts+" "+kpts+"   1 1 1\n")
    ftw.close()

