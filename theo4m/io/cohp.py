# -*- coding: utf-8 -*-
"""
Loads information from COHP (Crystal Orbital Hamilton Populations) files.
"""

from theo4m.atom.interaction import Interaction

__all__ = [
        'read_cohp',
        ]

def read_cohp(cohpfile, task, inter):
    fileobj = open(cohpfile)
    lines = fileobj.readlines()
    fileobj.close()
    if task == "ICOHPLIST":
        """Reads the ICOHPLIST.lobster output file, returning a collection of Interaction objects."""
        bonds = {}
        for i in range(1, len(lines)):
            bond = Interaction()
            bond.set_id(int(lines[i].split()[0]))
            bond.set_symbA(lines[i].split()[1])
            bond.set_symbB(lines[i].split()[2])
            bond.set_distance(float(lines[i].split()[3]))
            bond.set_transx(int(lines[i].split()[4]))
            bond.set_transy(int(lines[i].split()[5]))
            bond.set_transz(int(lines[i].split()[6]))
            bond.set_info(float(lines[i].split()[7]))
            bonds[bond.get_symbA()+"-"+bond.get_symbB()] = bond
        return bonds
    elif task == "COHPCAR":
        """Reads the COHPCAR.lobster output file, returning a COHPPlot object."""
        ncols = int(lines[1].split()[0])
        nspin = int(lines[1].split()[1])
        npts = int(lines[1].split()[2])
        iint = -1
        energy = []
        cohp = []
        for i in range(3, 3+ncols-1):
            if ":"+inter+"(" in lines[i]:
                iint = i - 3
                break
        if iint == -1: return None
        for i in range(ncols+2,ncols+2+npts):
            energy.append(float(lines[i].split()[0]))
            cohp.append(-float(lines[i].split()[(iint*2)+3]))
        return cohp, energy
    else:
        return None

