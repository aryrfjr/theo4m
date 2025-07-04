# -*- coding: utf-8 -*-
"""

    Read LAPWAtoms object(s) from file.

    =======================================  ================
     Package                                  Method
    =======================================  ================
     Elk                    output files      read_elk
    =======================================  ================

"""

__all__ = ['read_elk']

def read_elk(rundir, prefix, tasks, nuclides):
    from theo4m.nmr.lapw.io.elk import read_elk
    return read_elk(rundir, prefix, tasks, nuclides)

