#!/usr/bin/python

# python ~/dev-python/utils/PDOS/read_atomic.py 

import os.path
import matplotlib as mtp
import matplotlib.pyplot as plt
import sys

#
# The script
#

rundir = sys.argv[1]

# 0 - logarithmic derivatives
# 1 - wavefunctions and projectors for each momentum channel
# 2 - densities
# 3 - the projectors in reciprocal space
# 4 - the potential
# 5 - exchange potential
# ? - Fourier transforms of the kinetic energy and potential of the occupied states
wts = int(sys.argv[2])

multich = [0,1]

if wts in multich:
    nch = int(sys.argv[3])
    nrows = nch / 2 + nch % 2
    #plt.figure(1)
else:
    nch = 1

if wts == 0:
    plt.figure(num="LDs "+rundir)
    fileobj = open(rundir+"/ld1.dlog")
    linesae = fileobj.readlines()
    fileobj.close()
    fileobj = open(rundir+"/ld1ps.dlog")
    linesps = fileobj.readlines()
    fileobj.close()
    for ild in range(nch):
        ene = []
        ldae = []
        ldps = []
        for i, line in enumerate(linesae):
            ene.append(float(linesae[i].split()[0]))
            ldae.append(float(linesae[i].split()[ild+1]))
            ldps.append(float(linesps[i].split()[ild+1]))
        plt.subplot(nrows, 2, ild + 1)
        plt.plot(ene, ldae, "k-")
        plt.plot(ene, ldps, "r--")
elif wts == 1:
    plt.figure(num="PWs "+rundir)
    fileobj = open(rundir+"/ld1ps.wfc")
    lineswfcs = fileobj.readlines()
    fileobj.close()
    fileobj = open(rundir+"/beta")
    linesbeta = fileobj.readlines()
    fileobj.close()
    nwfcs = int(sys.argv[3])
    for ild in range(nch):
        r = []
        pwae = []
        pwps = []
        proj1 = []
        proj2 = []
        for i, line in enumerate(lineswfcs):
            if i == 0: continue
            r.append(float(lineswfcs[i].split()[0]))
            pwps.append(float(lineswfcs[i].split()[1+ild]))
            pwae.append(float(lineswfcs[i].split()[3+ild]))
            proj1.append(float(linesbeta[i].split()[1+(ild*2)]))
            proj2.append(float(linesbeta[i].split()[2+(ild*2)]))
        plt.subplot(nrows, 2, ild + 1)
        plt.plot(r, pwae, "k-")
        plt.plot(r, pwps, "r-")
        if int(sys.argv[4]):
            plt.plot(r, proj1, "r--")
            plt.plot(r, proj2, "r--")
elif wts == 2:
        plt.figure(num="density "+rundir)
        fileobj = open(rundir+"/density")
        lines = fileobj.readlines()
        fileobj.close()
        r = []
        coreps = []
        val = []
        core = []
        for i, line in enumerate(lines):
            r.append(float(lines[i].split()[0]))
            coreps.append(float(lines[i].split()[1]))
            val.append(float(lines[i].split()[2]))
            core.append(float(lines[i].split()[3]))
        plt.subplot(1, 1, 1)
        plt.plot(r, core, "k-")
        plt.plot(r, val, "r-")
        plt.plot(r, coreps, "k--")
        plt.axis([float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])])
plt.show()

