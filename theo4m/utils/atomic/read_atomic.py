#!/usr/bin/python

# python ~/dev-python/utils/PDOS/read_atomic.py 

import matplotlib as mtp
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/aryjr/dev-python")

from atomic.io import read_ld1_output
from atomic.atomicrun import AtomicRun

#
# The script
#

# What to show
wts = sys.argv[2]

# 
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

# Loading data
ar = read_ld1_output(sys.argv[1], sys.argv[3], sys.argv[4])

if wts == 'ps.wfc':
    # Plotting pseudo and AE KS orbitals
    labels, data = ar.get_wfc_ps()
    plt.figure(1)
    plt.subplot(1, 1, 1)
    #plt.text(labelsp[i][0], labelsp[i][1], labels[i])
    ic = 0
    lst = "--"
    tp = "-PS"
    for ip in range(len(labels)):
        plt.plot(data[0], data[ip+1], colors[ic]+lst, label=labels[ip]+tp)
        if ic < (len(labels)/2)-1:
            ic += 1
        else:
            ic = 0
            lst = "-"
            tp = "-AE"
    plt.legend()
    #plt.axis([0,2,-2,2])
    plt.show()
elif wts == 'dlog':
    # Plotting the logarithmic derivatives
    LDs = ar.get_log_deriv()
    data = LDs.get_dlog()
    dataps = LDs.get_dlog_ps()
    plt.figure(1)
    plt.subplot(1, 1, 1)
    for ip in range(len(data)-1):
        plt.plot(data[0], data[ip+1], "k--")
        plt.plot(data[0], dataps[ip+1], "r-")    
    plt.show()
elif wts == 'charge':
    aec = ar.get_ae_charge()
    psc = ar.get_ps_charge()
    plt.figure(1)
    plt.subplot(1, 2, 1)
    plt.plot(aec[0], aec[1], "k-")
    plt.plot(psc[0], psc[1], "k--")
    plt.axis([0.0,2.0,0.0,30.0])
    plt.subplot(1, 2, 2)
    plt.plot(aec[0], aec[2], "r-")
    plt.plot(aec[0], aec[3], "b-")
    plt.plot(psc[0], psc[2], "r--")
    plt.plot(psc[0], psc[3], "b--")
    plt.axis([0.0,2.0,0.0,30.0])
    plt.show()
elif wts == 'test':
    test = ar.get_test()
    els = test.get_energy_levels()
    bcs = test.get_bessel_cutoffs()
    bss = test.get_bessel_states_set()
    for i in range(len(els)):
        print "%s %.5f" % (els[i].get_nl(), els[i].get_AE() - els[i].get_PS())
    print ''
    wbss = len(bss) - 1
    for lb in range(3): # for l = 0, 1, and 2
        for i in range(len(els)):
            if (els[i].get_l() == lb):
                print "energy level %s = %.5f" % (els[i].get_nl(), els[i].get_PS())
        print "l = %d: %.5f %.5f %.5f" % (lb, bss[wbss][lb][0], bss[wbss][lb][1], bss[wbss][lb][2])
        print ''
    for i in range(1,len(bss)):
        print bcs[i], bss[i][0], bss[i][1], bss[i][2]

