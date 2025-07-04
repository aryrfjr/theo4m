"""

Miscellaneous database utils for ML models.

"""

import numpy as np
from ase import neighborlist
from ase import Atoms
from ase import Atom
import copy

def adjust_pbc_sa(iat, atoms):
    c = atoms.get_cell()
    g = atoms.get_positions()
    for ic in range(0, 3):
        if g[iat][ic] > c[ic][ic]: g[iat][ic] -= c[ic][ic]
        if g[iat][ic] < 0.0: g[iat][ic] += c[ic][ic]
    atoms.set_positions(g)

def adjust_pbc(atoms):
    c = atoms.get_cell()
    g = atoms.get_positions()
    for i in range(len(atoms)):
        for ic in range(0, 3):
            if g[i][ic] > c[ic][ic]: g[i][ic] -= c[ic][ic]
            if g[i][ic] < 0.0: g[i][ic] += c[ic][ic]
    atoms.set_positions(g)

def cryst2angs(atoms):
    c = atoms.get_cell()
    g = atoms.get_positions()
    g = np.dot(g, c)
    atoms.set_positions(g)

# TODO: working for a ortogonal cell
def angs2cryst(atoms):
    c = atoms.get_cell()
    g = atoms.get_positions()
    for i in range(len(atoms)):
        g[i][0] /= c[0][0]
        g[i][1] /= c[1][1]
        g[i][2] /= c[2][2]
    atoms.set_positions(g)

# Dump a cell in a ase.Atoms object (adump) in the file fdump
def dump_supercell(fdump, adump, title="title", spc_symbs=[], output="xyz", coord_type="angs"):
    ocell = adump.get_cell()
    symbs = adump.get_chemical_symbols()
    geom = adump.get_positions()
    if output == "xyz":
        fdump.write("%d\n" % len(adump))
        # The order of the lattice vector was wrong in a previous version 
        # but not now it is correct. See more at:
        # https://libatoms.github.io/QUIP/io.html#module-ase.io.extxyz
        fdump.write("Lattice=\"%10f %10f %10f  %10f %10f %10f  %10f %10f %10f\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n" 
                  % (ocell[0,0], ocell[0,1], ocell[0,2], ocell[1,0], ocell[1,1], ocell[1,2], ocell[2,0], ocell[2,1], ocell[2,2]))
        for i in range(len(adump)):
            fdump.write("%s %10f %10f %10f\n" % (symbs[i], geom[i][0], geom[i][1], geom[i][2]))
    else: # dumping a .pw.inp file
        fdump.write(" &CONTROL\n")
        fdump.write("                       title = '%s',\n" % title)
        fdump.write("                 calculation = 'scf',\n")
        fdump.write("                restart_mode = 'from_scratch',\n")
        fdump.write("                      outdir = '.',\n")
        fdump.write("                      wfcdir = '.',\n")
        fdump.write("                  pseudo_dir = '.',\n")
        fdump.write("                      prefix = '%s',\n" % title)
        fdump.write("                   verbosity = 'high',\n")
        fdump.write("                     tprnfor = .true.,\n")
        fdump.write("                  wf_collect = .true.\n")
        fdump.write(" /\n")
        fdump.write(" &SYSTEM\n")
        fdump.write("                       ibrav = 0,\n")
        fdump.write("                   celldm(1) = 1.8897265,\n")
        fdump.write("                         nat = %d,\n" % len(adump))
        fdump.write("                        ntyp = %d,\n" % len(spc_symbs))
        fdump.write("                     ecutwfc = 0,\n")
        fdump.write("                     ecutrho = 0,\n")
        fdump.write(" /\n")
        fdump.write("CELL_PARAMETERS alat\n")
        fdump.write("%10f %10f %10f\n" % (ocell[0,0], ocell[0,1], ocell[0,2]))
        fdump.write("%10f %10f %10f\n" % (ocell[1,0], ocell[1,1], ocell[1,2]))
        fdump.write("%10f %10f %10f\n" % (ocell[2,0], ocell[2,1], ocell[2,2]))
        fdump.write("ATOMIC_SPECIES\n")
        for i in range(len(spc_symbs)):
            fdump.write("   %s  0.0   psp.UPF\n" % spc_symbs[i])

        if coord_type == "angs":
            fdump.write("ATOMIC_POSITIONS\n")
        else:
            fdump.write("ATOMIC_POSITIONS crystal\n")
        for i in range(len(adump)):
            fdump.write("%s %10f %10f %10f\n" % (symbs[i], geom[i][0], geom[i][1], geom[i][2]))
        fdump.write("K_POINTS automatic\n")
        fdump.write("  1 1 1   1 1 1 \n")

#############################################################################################
#
# Originally devised for the comparison between the 
# databases generated with EXAFS-guided MC geometry 
# optimization and snapshots from QfM CMD simulations. 
# See more on 18/06/2020.
#
#############################################################################################
def get_clusters_from_cell(atoms, out_cell, centered = True):
    clusters = []
    nl = neighborlist.NeighborList([3.75/2]*len(atoms), skin=0.0, bothways=True, self_interaction=False)
    nl.update(atoms)
    csymbs = atoms.get_chemical_symbols()
    atomsp = atoms.get_positions()
    atomst = atoms.get_tags()
    for icl in range(len(atoms)):
        indices, offsets = nl.get_neighbors(icl)
        cluster = Atoms(cell=out_cell, pbc=not centered)
        # the central atom
        cluster.append(Atom(csymbs[icl], copy.deepcopy(atomsp[icl]), tag=atomst[icl]))
        # and its first-neighbors
        for k, offset in zip(indices, offsets):
            #cluster.append(Atom(csymbs[k], copy.deepcopy(atomsp[k]) + np.dot(offset, atoms.get_cell()), tag=atomst[k]))
            cluster.append(Atom(csymbs[k], copy.deepcopy(atomsp[k]), tag=atomst[k]))
        # translation of the cluster to the center of the cell
        # expected to be cubic if pbc == False
        if centered:
            cluster.translate(-atomsp[icl])
            center = out_cell[0][0]/2.0
            cluster.translate([center, center, center])
        clusters.append(cluster)
    return clusters

#############################################################################################
#
# Originally devised for the comparison between the 
# databases generated with EXAFS-guided MC geometry 
# optimization and snapshots from QfM CMD simulations. 
# See more on 11/06/2020.
#
#############################################################################################

# Returns an item of an ordinary extended .xyz with extra parameters
# https://web.archive.org/web/20190811094343/https://libatoms.github.io/QUIP/io.html#extendedxyz
def get_item_ext_xyz(atoms, params, frac, 
                     extra_propstr = "", extra_prop = []): # TODO generalize extra properties
    item_str = ""
    item_str += ("%d\n" % len(atoms))
    # The order of the lattice vector was wrong in a previous version 
    # but not now it is correct. See more at:
    # https://libatoms.github.io/QUIP/io.html#module-ase.io.extxyz
    cell = atoms.get_cell()
    parstr = ("Lattice=\"%10f %10f %10f  %10f %10f %10f  %10f %10f %10f\" " 
              % (cell[0,0], cell[0,1], cell[0,2], cell[1,0], cell[1,1], cell[1,2], cell[2,0], cell[2,1], cell[2,2]))
    parstr += "Properties=species:S:1:pos:R:3"
    if extra_propstr != "": parstr += ":" + extra_propstr
    parstr += " "
    for kp in params.keys():
        if isinstance(params[kp], str): # if it is a string
            parstr += kp+"=\""+params[kp]+"\" "
        else:
            parstr += kp+"="+str(params[kp])+" "
    if atoms.get_pbc()[0] and atoms.get_pbc()[1] and atoms.get_pbc()[2]:
        parstr += "pbc=\"T T T\" "
    else:
        parstr += "pbc=\"F F F\" "
    item_str += (parstr+"\n")
    symbs = atoms.get_chemical_symbols()
    geom = atoms.get_positions()
    if frac: geom = np.dot(geom, cell)
    for i in range(len(atoms)):
        item_str += ("%s %10f %10f %10f" % (symbs[i], geom[i][0], geom[i][1], geom[i][2]))
        # TODO generalize extra properties
        item_str += (" %15f %15f %15f\n" % (extra_prop[i][0], extra_prop[i][1], extra_prop[i][2]))
    return item_str

# Writing an ordinary extended .xyz with extra parameters
# https://web.archive.org/web/20190811094343/https://libatoms.github.io/QUIP/io.html#extendedxyz
def write_ext_xyz(dir2w, file_name, atoms, params, append = True):
    if append:
        wtd = "a"
    else:
        wtd = "w"
    ftw = open(dir2w+"/"+file_name+".xyz", wtd)
    ftw.write(get_item_ext_xyz(atoms, params, False))
    ftw.close()

