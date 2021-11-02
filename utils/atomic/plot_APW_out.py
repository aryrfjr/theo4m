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

# What to show (https://docs.abinit.org/tutorial/paw3/)
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
    for ild in range(nch):
        fileobj = open(rundir+"/logderiv."+str(ild))
        lines = fileobj.readlines()
        fileobj.close()
        ene = []
        ldae = []
        ldps = []
        for i, line in enumerate(lines):
            ene.append(float(lines[i].split()[0]))
            ldae.append(float(lines[i].split()[1]))
            ldps.append(float(lines[i].split()[2]))
        plt.subplot(nrows, 2, ild + 1)
        plt.plot(ene, ldae, "k-")
        plt.plot(ene, ldps, "r--")
elif wts == 1:
    plt.figure(num="PWs "+rundir)
    for ild in range(nch):
        fileobj = open(rundir+"/wfn"+str(ild+1))
        lines = fileobj.readlines()
        fileobj.close()
        r = []
        pwae = []
        pwps = []
        proj = []
        for i, line in enumerate(lines):
            if i == 0: continue
            r.append(float(lines[i].split()[0]))
            pwae.append(float(lines[i].split()[1]))
            pwps.append(float(lines[i].split()[2]))
            proj.append(float(lines[i].split()[3]))
        plt.subplot(nrows, 2, ild + 1)
        plt.plot(r, pwae, "k-")
        plt.plot(r, pwps, "r-")
        if int(sys.argv[4]): plt.plot(r, proj, "r--")
elif wts == 2:
        plt.figure(num="density "+rundir)
        fileobj = open(rundir+"/density")
        lines = fileobj.readlines()
        fileobj.close()
        r = []
        core = []
        val = []
        coreps = []
        valps = []
        for i, line in enumerate(lines):
            r.append(float(lines[i].split()[0]))
            core.append(float(lines[i].split()[1]))
            val.append(float(lines[i].split()[2]))
            coreps.append(float(lines[i].split()[3]))
            valps.append(float(lines[i].split()[4]))
        plt.subplot(1, 1, 1)
        plt.plot(r, core, "k-")
        plt.plot(r, val, "r-")
        plt.plot(r, coreps, "k--")
        plt.plot(r, valps, "r--")
        plt.axis([float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])])
elif wts == 3:
        plt.figure(num="proj. in rec. space "+rundir)
        maxnw = int(sys.argv[3])
        plt.subplot(1, 1, 1)
        colors = ["k-","r-","g-","b-","m-"]
        ic = 0
        for iw in range(maxnw):
            if os.path.isfile(rundir+"/tprod."+str(iw+1)):
                fileobj = open(rundir+"/tprod."+str(iw+1))
                lines = fileobj.readlines()
                fileobj.close()
                r = []
                tp = []
                for i, line in enumerate(lines):
                    r.append(float(lines[i].split()[0]))
                    tp.append(float(lines[i].split()[1]))
                plt.plot(r, tp, colors[ic])
                ic += 1
elif wts == 4:
        plt.figure(num="potential "+rundir)
        fileobj = open(rundir+"/potential")
        lines = fileobj.readlines()
        fileobj.close()
        r = []
        aep = []
        psp = []
        for i, line in enumerate(lines):
            r.append(float(lines[i].split()[0]))
            aep.append(float(lines[i].split()[1]))
            psp.append(float(lines[i].split()[2]))
        plt.subplot(1, 1, 1)
        plt.plot(r, aep, "k-")
        plt.plot(r, psp, "r--")
        plt.axis([float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])])
elif wts == 5:
        plt.figure(num="rVx "+rundir)
        fileobj = open(rundir+"/rVx")
        lines = fileobj.readlines()
        fileobj.close()
        r = []
        aep = []
        psp = []
        for i, line in enumerate(lines):
            r.append(float(lines[i].split()[0]))
            aep.append(float(lines[i].split()[1]))
            psp.append(float(lines[i].split()[2]))
        plt.subplot(1, 1, 1)
        plt.plot(r, aep, "k-")
        plt.plot(r, psp, "r--")
        plt.axis([float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])])
plt.show()

