# -*- coding: utf-8 -*-
"""

    Read GIPAWAtoms object(s) from file.

    =======================================  ================
     Package                                  Method
    =======================================  ================
     Quantum Espresso GIPAW output files      read_qe
    =======================================  ================

"""

__all__ = ['read_qe']

def read_qe(rundir, prefix, tasks, nuclide_label = "", core_relax = ""):
    from theo4m.nmr.gipaw.io.qe import read_qe
    return read_qe(rundir, prefix, tasks, nuclide_label = nuclide_label, core_relax = core_relax)

