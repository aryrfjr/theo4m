# -*- coding: utf-8 -*-
"""
Loads information from Quantum Espresso files.
"""
from numpy import array
from numpy import dot
from ase import Atom
from ase import Atoms

__all__ = [
        "read_qe",
        ]

def read_qe(qefile, task):
    """Reads the pw.x input/output files, returning an ase.Atoms object."""
    fileobj = open(qefile)
    lines = fileobj.readlines()
    fileobj.close()
    if task == "PW_INP": # Reading a pw.x input file
        for i, line in enumerate(lines):
            if "nat" in line:
            # Reading the number of atoms in the cell
                if "," in line.split()[2]:
                    nat = int(line.split()[2][:len(line.split()[2])-1])
                else:
                    nat = int(line.split()[2])
            elif "ntyp" in line:
                if "," in line.split()[2]:
                    ntypat = int(line.split()[2][:len(line.split()[2])-1])
                else:
                    ntypat = int(line.split()[2])
            elif "CELL_PARAMETERS" in line:
            # Reading the cell vectors
                cell = [x.split()[0:3] for x in lines[i + 1:i + 4]]
                cell = array([[float(col) for col in row] for row in cell])
            elif "ATOMIC_POSITIONS" in line:
                if "crystal" in line:
                # Reading the atoms and creating a collection of ase.Atoms objects
                    geom_start = i + 1
                    geom_stop = geom_start + nat
                    species = [line.split()[0] for line in lines[geom_start:geom_stop]]
                    geom = dot(array([[float(col) for col in line.split()[1:4]]
                        for line in lines[geom_start:geom_stop]]), cell)
                else:
                # Reading the atoms and creating a collection of ase.Atoms objects
                    geom_start = i + 1
                    geom_stop = geom_start + nat
                    species = [line.split()[0] for line in lines[geom_start:geom_stop]]
                    geom = array([[float(col) for col in line.split()[1:4]]
                        for line in lines[geom_start:geom_stop]])
        # Returning the input structure
        rstrc = Atoms(
            cell=cell,
            pbc=True,
            positions=geom,
            symbols="".join(species))
        return rstrc
    elif task == "PW_OUT_RELAX": # Reading a pw.x output file for a calculation = "relax"
        status = "NONE"
        rstrcs = []
        rtotEs = []
        rtotFs = []
        rforces = []
        rstress = []
        for i, line in enumerate(lines):
            # Initial information related to the input cell
            if "number of atoms/cell" in line:
            # Reading the number of atoms in the cell
                nat = int(line.split()[4])
            elif "number of atomic types" in line:
                ntypat = int(line.split()[5])
            elif "crystal axes: (cart. coord. in units of alat)" in line:
            # Reading the cell vectors
                cell = [x.split()[3:6] for x in lines[i + 1:i + 4]]
                cell = array([[float(col) for col in row] for row in cell])
            elif "Crystallographic axes" in line:
            # Reading the input coordinates and creating a collection of ase.Atoms objects
                geom_start = i + 3
                geom_stop = geom_start + nat
                species = [line.split()[1] for line in lines[geom_start:geom_stop]]
                geom = dot(array([[float(col) for col in line.split()[6:9]]
                    for line in lines[geom_start:geom_stop]]), cell)
                tstrc = Atoms(
                    cell=cell,
                    pbc=True,
                    positions=geom,
                    symbols="".join(species))
                rstrcs.append(tstrc)
                #print ("Appending coordinates (first)")
            # Now, just after each SCF cycle
            # Reading total energy
            elif "Forces acting on atoms" in line:
                forces_start = i + 2
                forces_stop = forces_start + nat
                try:
                    forces = array([[float(col) for col in line.split()[6:9]]
                        for line in lines[forces_start:forces_stop]])
                    #print ("Appending forces")
                    rforces.append(forces)
                except ValueError:
                    # expected to occur when forces are too big
                    # and so incompatible with the format used in QE
                    # for instance:
                    # atom    3 type  2   force =   674.57999165  312.30521069-1079.69944125
                    print ("Rerror reading forces in file:")
                    print (qefile)
                    #print ("Appending forces (empty)")
                    rforces.append([])
            elif "!    total energy" in line:
                rtotEs.append(float(line.split()[4]))
                #print ("Appending energy")
            elif "total   stress  (Ry/bohr**3)" in line:
            # Reading the stress tensor
                stress = [x.split()[0:3] for x in lines[i + 1:i + 4]]
                stress = array([[float(col) for col in row] for row in stress])
                rstress.append(stress)
                #print ("Appending stress")
            elif "Total force" in line:
                rtotFs.append(float(line.split()[3]))
                #print ("Appending total forces")
            elif "ATOMIC_POSITIONS (alat)" in line:
            # Reading the relaxed and creating a collection of ase.Atoms objects
                geom_start = i + 1
                geom_stop = geom_start + nat
                species = [line.split()[0] for line in lines[geom_start:geom_stop]]
                geom = array([[float(col) for col in line.split()[1:4]]
                    for line in lines[geom_start:geom_stop]])
                tstrc = Atoms(
                    cell=cell,
                    pbc=True,
                    positions=geom,
                    symbols="".join(species))
                rstrcs.append(tstrc)
                #print ("Appending coordinates")
            elif "convergence NOT achieved after 100 iterations: stopping" in line:
            # Removing the last item the vector with structures
                status = "SCF_NOT_CONVERGED"
                rstrcs.pop()
                #print ("Removing coordinates")
        # Checking if no even the first SCF started
        if len(rtotEs) == 0 and status == "NONE":
            status = "CRASH"
            rstrcs.pop()
            #print ("Removing coordinates")
        # Checking if the SCF has not been finished because of timeout
        if len(rstrcs) > len(rtotEs) and status == "NONE":
            status = "TIMEOUT_OR_CRASH"
            rstrcs.pop()
            #print ("Removing coordinates")
        # Checking if the BFGS has been finished
        if status == "TIMEOUT_OR_CRASH" and "JOB DONE" in lines[len(lines)-2]:
            status = "FINISHED"
        # Returning a collection of cells and properties
        return status, rstrcs, rtotEs, rtotFs, rforces, rstress

