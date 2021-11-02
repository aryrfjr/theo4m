# -*- coding: utf-8 -*-
"""

    =============================================  ================================
     Package                                        Method
    =============================================  ================================
     Deprecated! Use read_lammps_dump instead       read_lammps
     LAMMPS .dump files                             read_lammps_dump
     LAMMPS .lmp.out files                          read_lammps_out
     LAMMPS .data files                             read_lammps_data
     Quantum Espresso input/output files            read_qe
     COHP (Crystal Orbital Hamilton Populations)    read_cohp
     Extended .xyz files                            read_extended_xyz
    =============================================  ================================

"""

__all__ = [
        "read_lammps", 
        "read_lammps_dump", 
        "read_lammps_out", 
        "read_lammps_data", 
        "read_qe", 
        "read_cohp", 
        "read_extended_xyz"
        ]

def read_lammps(lmpoutput = "", items = [-1], spc_symbs = [], frac = False):
    from theo4m.io.lammps import read_lammps
    return read_lammps(lmpoutput, items, spc_symbs, frac)

def read_lammps_dump(lmpdump_file = "", items = [-1], spc_symbs = [], frac = False):
    from theo4m.io.lammps import read_lammps_dump
    return read_lammps_dump(lmpdump_file, items, spc_symbs, frac)

def read_lammps_out(lmpout_file = ""):
    from theo4m.io.lammps import read_lammps_out
    return read_lammps_out(lmpout_file)

def read_lammps_data(lmpdata_file = "", spc_symbs = [], frac = False):
    from theo4m.io.lammps import read_lammps_data
    return read_lammps_data(lmpdata_file, spc_symbs, frac)

def read_qe(qefile = "", task = "PW_OUT"):
    from theo4m.io.qe import read_qe
    return read_qe(qefile, task)

def read_cohp(cohpfile = "", task = "ICOHPLIST", inter = ""):
    from theo4m.io.cohp import read_cohp
    return read_cohp(cohpfile, task, inter)

def read_extended_xyz(xyzfile = ""):
    from theo4m.io.xyz import read_extended_xyz
    return read_extended_xyz(xyzfile)

