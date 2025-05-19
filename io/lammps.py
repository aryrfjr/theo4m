"""

Loads a supercell from LAMMPS output files.

According to the documentation, the commands thermodynamics, dump, and 
fix do write output files. Moreover, various output-related commands work 
with three different styles of data:
 - a global datum is one or more system-wide values, e.g. the temperature 
   of the system;
 - a per-atom datum is one or more values per atom, e.g. the kinetic energy 
   of each atom;
 - local datums are calculated by each processor based on the atoms it owns, 
   but there may be zero or more per atom, e.g. a list of bond distances.

"""
from numpy import array
from numpy import dot
from ase import Atom
from ase import Atoms

__all__ = [
        "read_lammps", 
        "read_lammps_dump", 
        "read_lammps_out"
        "read_lammps_data"
        ]

"""

items can be: [-1] (all), [-2] (last), or [#, ...] (a specific set of items)

Deprecated! Use read_lammps_dump instead.

"""
def read_lammps(lmpoutput = "", items = [-1], spc_symbs = [], frac = False):
    return __dump_atom__all__(lmpoutput, items, spc_symbs, frac)

"""

items can be: [-1] (all), [-2] (last), or [#, ...] (a specific set of items)

"""
def read_lammps_dump(lmpdump_file = "", items = [-1], spc_symbs = [], frac = False):
    # Note: frac = True means that the LAMMPS dump file has fractional 
    # coordinates and the returned positions must be absolute.
    return __dump_atom__all__(lmpdump_file, items, spc_symbs, frac)

# Reading the .lmp.out file
def read_lammps_out(lmpout_file = ""):
    heads = []
    blocks = []
    fileobj = open(lmpout_file)
    lines = fileobj.readlines()
    fileobj.close()
    inblock = False
    for il, line in enumerate(lines):
        if "Step" in line: # creating a block
            inblock = True
            heads.append(line.split())
            block = []
        elif "Loop time" in line and inblock: # finishing a block
            blocks.append(block)
            inblock = False
        elif inblock: 
            block.append(line.split())
    return heads, blocks

# Reading the .data file
# https://lammps.sandia.gov/doc/read_data.html
# TODO: make it compatible with the LAMMPS read_data command above
def read_lammps_data(lmpdata_file = "", spc_symbs = [], frac = False):
    fileobj = open(lmpout_file)
    lines = fileobj.readlines()
    fileobj.close()


# Processing an individual item in a .dump file
def __process_item__(i, lines, spc_symbs, frac):
    # the LAMMPS step
    step = lines[i+1].split()[0]
    # number of atoms
    nat = int(lines[i+3].split()[0])
    # the cell vectors of the cell centered at (0,0,0)
    a = [float(lines[i+5].split()[0]), float(lines[i+5].split()[1])]
    b = [float(lines[i+6].split()[0]), float(lines[i+6].split()[1])]
    c = [float(lines[i+7].split()[0]), float(lines[i+7].split()[1])]
    cell = array([a, b, c])
    # the cubic cell with origin at (0,0,0)
    ac = cell[0][1] - cell[0][0]
    bc = cell[1][1] - cell[1][0]
    cc = cell[2][1] - cell[2][0]
    ocell = array([[ac,0.0,0.0],[0.0,bc,0.0],[0.0,0.0,cc]])
    # atomic coordinates
    lstart = i+9
    lstop = lstart+nat
    ca_ids = [int(line.split()[0]) for line in lines[lstart:lstop]]
    ca_symbs = [spc_symbs[int(line.split()[1])-1] for line in lines[lstart:lstop]]
    if frac:
        geom = dot(array([[float(col) for col in line.split()[2:5]]
            for line in lines[lstart:lstop]]), ocell)
    else:
        geom = array([[float(col) for col in line.split()[2:5]]
            for line in lines[lstart:lstop]])
    atoms = Atoms(cell=ocell,pbc=True)
    for i in range(nat):
        atoms.append(Atom(ca_symbs[i], geom[i], tag=ca_ids[i]))
    return step, atoms, ocell

# Reading the .dump file
def __dump_atom__all__(lmpdump_file, items, spc_symbs, frac):
    fileobj = open(lmpdump_file)
    lines = fileobj.readlines()
    fileobj.close()
    rsteps = []
    ratoms = []
    rcell = []
    iitem = 1
    if items[0] == -2: # only the last item
        nati = int(lines[3].split()[0]) # the number of atoms of the first item
        il = len(lines) - (nati + 9)
        istep, iatoms, icell = __process_item__(il, lines, spc_symbs, frac)
        rsteps.append(istep)
        ratoms.append(iatoms)
        rcell.append(icell)
    else:
        for il, line in enumerate(lines):
            if "ITEM: TIMESTEP" in line:
                if items[0] == -1 or iitem in items:
                    istep, iatoms, icell = __process_item__(il, lines, spc_symbs, frac)
                    print ("Reading step # %s ..." % istep)
                    rsteps.append(istep)
                    ratoms.append(iatoms)
                    rcell.append(icell)
                iitem += 1
    return rsteps, ratoms, rcell

