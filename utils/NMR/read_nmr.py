#!/usr/bin/python

import sys
sys.path.append("/home/aryjr/dev-python")

from numpy import array
from theo4m.nmr.gipaw.io import read_qe
from theo4m.nmr.gipaw.gipawrun import GIPAWRun
from theo4m.nmr.gipaw.gipawatom import GIPAWAtom
from theo4m.nmr.lapw.io import read_elk
from theo4m.nmr.lapw.lapwrun import LAPWRun
from theo4m.nmr.lapw.lapwatom import LAPWAtom
from theo4m.nmr.nuclide import Nuclide

#
# The script
#

#
# What to do?
# 0 - Show everything
# 1 - Show only the experimental value followed by the different XC values
# 2 - Show Tables with shifts
# 3 - Show Tables with EFG parameters
wtd = int(sys.argv[2])

if wtd == 0:
    wss = int(sys.argv[3])
elif wtd == 1:
    #
    # What to show
    # 0 - K_iso
    # 1 - Cq
    # 2 - sigma_s
    # 3 - eta
    wts = int(sys.argv[3])
    #
    # How density was computed
    # 0 - nart
    # 1 - art
    dm = int(sys.argv[4])
    #
    # sites to show
    hms = int(sys.argv[5]) # how many sites
    sts = array([col for col in sys.argv[6:hms+6]])
    if wts == 2:
        wxc = sys.argv[hms+6] # which XC

#
# Reading the input file
ifile = sys.argv[1]
fileobj = open(ifile)
lines = fileobj.readlines()
fileobj.close()

prefix = lines[0].split()[0]
nl = 1
nn = int(lines[nl].split()[0]) # number of nuclides
nl += 1
nuclides = [None]*nn
for inn in range(nn):
    nuclides[inn] = Nuclide()
    nuclides[inn].set_symbol(lines[nl].split()[0])
    nuclides[inn].set_quadrupolar_moment(float(lines[nl].split()[1]))
    nuclides[inn].set_nuclear_spin(float(lines[nl].split()[2]))
    nl += 1
if lines[nl].split()[0] == "QE-GIPAW":
    nsites = int(lines[nl].split()[1]) # number of different sites of this nuclide and their ids
    sites_qe = array([int(col) for col in lines[nl].split()[2:nsites+2]]) # number of different sites of this nuclide and their ids
    nl += 1
if lines[nl].split()[0] == "Elk":
    sites_elk = array([int(col) for col in lines[nl].split()[1:nsites+1]]) # sites of this nuclide and their ids
    nl += 1
# experimental values of Cq, eta, and Kiso
Cq_exp = [0.0]*nsites
eta_exp = [0.0]*nsites
Kiso_exp = [0.0]*nsites
for ist in range(nsites):
    Cq_exp[ist] = float(lines[nl].split()[0])
    eta_exp[ist] = float(lines[nl].split()[1])
    Kiso_exp[ist] = float(lines[nl].split()[2])
    nl += 1
# ref prefix, id and its experimental chemical shift
ref_prefix = lines[nl].split()[0]
ref_id = int(lines[nl].split()[1])
ref_delta_exp = float(lines[nl].split()[2])
nl += 1
nsim = int(lines[nl].split()[0]) # number of different simulations
nl += 1
sim_label = [""]*nsim
sim_type = [""]*nsim
sim_paw = [""]*nsim
sim_crm = [""]*nsim
sim_ss_dir = [""]*nsim
sim_so_dir = [""]*nsim
sim_ref_dir = [""]*nsim
crms = [] # the core_relax_method
for ism in range(nsim):
    sim_label[ism] = lines[nl].split()[0]
    sim_type[ism] = lines[nl].split()[1]
    sim_ss_dir[ism] = lines[nl+1].split()[0]
    if sim_type[ism] == "QE-GIPAW":
        sim_paw[ism] = lines[nl].split()[2]
        sim_crm[ism] = int(lines[nl].split()[3])
        if not sim_crm[ism] in crms: crms.append(sim_crm[ism])
        sim_so_dir[ism] = lines[nl+2].split()[0]
        sim_ref_dir[ism] = lines[nl+3].split()[0]
        nl += 4
    elif sim_type[ism] == "Elk":
        nl += 2

#
# Now loading the data for each simulation
sim_ss_run = [None]*nsim
sim_so_run = [None]*nsim
sim_ref_run = [None]*nsim
for ism in range(nsim):
    if sim_type[ism] == "QE-GIPAW":
        sim_ss_run[ism] = read_qe(sim_ss_dir[ism], prefix, ["hyperfine_art","hyperfine_nart", "efg"], nuclides, crms)
        sim_ss_run[ism].compute_efg_parameters() # not I will pass a Nuclide object
        sim_ss_run[ism].compute_sigma_s(True, True, crms)
        sim_so_run[ism] = read_qe(sim_so_dir[ism], prefix, ["nmr"], nuclides, crms)
        sim_ref_run[ism] = read_qe(sim_ref_dir[ism], ref_prefix, ["nmr"], nuclides, crms)
    elif sim_type[ism] == "Elk":
        sim_ss_run[ism] = read_elk(sim_ss_dir[ism], prefix, ["hyperfine", "efg"], nuclides)
        sim_ss_run[ism].compute_efg_parameters()
        sim_ss_run[ism].compute_sigma_s()

#
# And finally the output
if wtd == 0:
    # Kiso
    for ist in range(nsites):
        if wss == -1 or wss == sites_qe[ist]:
            site = sim_ss_run[0].get_atom(sites_qe[ist])
            print "K_iso:", site.symbol+"("+str(sites_qe[ist])+")"
            print "%50s %30s %15s" % ("nart", "art", "exp.")
            for ism in range(nsim):
                if sim_type[ism] == "QE-GIPAW":
                    sref = sim_ref_run[ism].get_atom(ref_id).get_sigma_o()
                    sigma_sn = sim_ss_run[ism].get_atom(sites_qe[ist]).get_sigma_s_nart(sim_crm[ism] - 1)
                    sigma_sa = sim_ss_run[ism].get_atom(sites_qe[ist]).get_sigma_s_art(sim_crm[ism] - 1)
                    sigma_o = sim_so_run[ism].get_atom(sites_qe[ist]).get_sigma_o()
                    Kison = (sref + ref_delta_exp) - (sigma_sn + sigma_o)
                    Kisoa = (sref + ref_delta_exp) - (sigma_sa + sigma_o)
                    print "%10s %15s %3d %15.2f (%9.2f)  %15.2f (%9.2f) %10.2f " % (sim_label[ism], 
                               sim_paw[ism], sim_crm[ism], Kison, (Kiso_exp[ist] - Kison), Kisoa, (Kiso_exp[ist] - Kisoa), Kiso_exp[ist])
            print ""
    
    # Cq
    for ist in range(nsites):
        if wss == -1 or wss == sites_qe[ist]:
            site = sim_ss_run[0].get_atom(sites_qe[ist])
            print "Cq:", site.symbol+"("+str(sites_qe[ist])+")"
            for ism in range(nsim): 
                if sim_type[ism] == "QE-GIPAW":
                    sites = sites_qe
                elif sim_type[ism] == "Elk":
                    sites = sites_elk
                print "%10s %10s %15s %10.2f (%6.2f) " % (sim_type[ism], sim_label[ism], sim_paw[ism], 
                                           sim_ss_run[ism].get_atom(sites[ist]).get_Cq(), 
                                           (Cq_exp[ist] - abs(sim_ss_run[ism].get_atom(sites[ist]).get_Cq())))
            print ""
    
    # eta
    for ist in range(nsites):
        if wss == -1 or wss == sites_qe[ist]:
            site = sim_ss_run[0].get_atom(sites_qe[ist])
            print "eta:", site.symbol+"("+str(sites_qe[ist])+")"
            for ism in range(nsim): 
                if sim_type[ism] == "QE-GIPAW":
                    sites = sites_qe
                elif sim_type[ism] == "Elk":
                    sites = sites_elk
                print "%10s %10s %15s %10.2f (%6.2f) " % (sim_type[ism], sim_label[ism], sim_paw[ism], 
                                           sim_ss_run[ism].get_atom(sites[ist]).get_eta(), 
                                           (eta_exp[ist] - abs(sim_ss_run[ism].get_atom(sites[ist]).get_eta())))
            print ""
    
    # sigma_s
    for ist in range(nsites):
        if wss == -1 or wss == sites_qe[ist]:
            site = sim_ss_run[0].get_atom(sites_qe[ist])
            print "sigma_s:", site.symbol+"("+str(sites_qe[ist])+")"
            print "%60s %30s" % ("nart", "art")
            for ism in range(nsim):
                if sim_type[ism] == "QE-GIPAW":
                    # checking if there is a Elk counterpart
                    ielk = -1
                    for isme in range(nsim):
                        if sim_type[isme] == "Elk" and sim_label[ism] == sim_label[isme]:
                            ielk = isme
                            break
                    if ielk > -1:
                        sigma_sn_qe = sim_ss_run[ism].get_atom(sites_qe[ist]).get_sigma_s_nart(sim_crm[ism] - 1)
                        sigma_sa_qe = sim_ss_run[ism].get_atom(sites_qe[ist]).get_sigma_s_art(sim_crm[ism] - 1)
                        sigma_s_elk = sim_ss_run[ielk].get_atom(sites_elk[ist]).get_sigma_s()
                        print "%10s %10s %15s %3d %15.2f (%9.2f) %15.2f (%9.2f) " % (sim_type[ism], sim_label[ism], sim_paw[ism], sim_crm[ism], 
                                       sigma_sn_qe, ((sigma_s_elk - sigma_sn_qe)), sigma_sa_qe, ((sigma_s_elk - sigma_sa_qe)))
                    else:
                        print "%10s %10s %15s %3d %15.2f %15.2f" % (sim_type[ism], sim_label[ism], sim_paw[ism], sim_crm[ism], 
                                           sim_ss_run[ism].get_atom(sites_qe[ist]).get_sigma_s_nart(sim_crm[ism] - 1),
                                           sim_ss_run[ism].get_atom(sites_qe[ist]).get_sigma_s_art(sim_crm[ism] - 1))
                elif sim_type[ism] == "Elk":
                    print "%10s %10s %19s %15.2f" % (sim_type[ism], sim_label[ism], "", 
                                       sim_ss_run[ism].get_atom(sites_elk[ist]).get_sigma_s())
            print ""
    
    # Kax
    for ist in range(nsites):
        if wss == -1 or wss == sites_qe[ist]:
            site = sim_ss_run[0].get_atom(sites_qe[ist])
            print "Kax:", site.symbol+"("+str(sites_qe[ist])+")"
            print "%55s %15s" % ("nart", "art")
            for ism in range(nsim):
                if sim_type[ism] == "QE-GIPAW":
                    print "%10s %10s %15s %3d %15.2f %15.2f" % (sim_type[ism], sim_label[ism], sim_paw[ism], sim_crm[ism], 
                                       sim_ss_run[ism].get_atom(sites_qe[ist]).get_Kax_nart(sim_crm[ism] - 1),
                                       sim_ss_run[ism].get_atom(sites_qe[ist]).get_Kax_art(sim_crm[ism] - 1))
            print ""

    # sigma_o
    for ist in range(nsites):
        if wss == -1 or wss == sites_qe[ist]:
            site = sim_ss_run[0].get_atom(sites_qe[ist])
            print "sigma_o:", site.symbol+"("+str(sites_qe[ist])+")"
            for ism in range(nsim):
                if sim_type[ism] == "QE-GIPAW":
                    sref = sim_ref_run[ism].get_atom(ref_id).get_sigma_o()
                    sigma_o = sim_so_run[ism].get_atom(sites_qe[ist]).get_sigma_o()
                    dso = sigma_o # (sref + ref_delta_exp) - sigma_o
                    print "%10s %15s %3d %15.2f" % (sim_label[ism], 
                               sim_paw[ism], sim_crm[ism], dso)
            print ""

    # sigma_s_bare
    for ist in range(nsites):
        if wss == -1 or wss == sites_qe[ist]:
            site = sim_ss_run[0].get_atom(sites_qe[ist])
            print "sigma_s_bare:", site.symbol+"("+str(sites_qe[ist])+")"
            for ism in range(nsim):
                if sim_type[ism] == "QE-GIPAW":
                    print "%10s %10s %15s %3d %15.2f" % (sim_type[ism], sim_label[ism], sim_paw[ism], sim_crm[ism], 
                                       sim_ss_run[ism].get_atom(sites_qe[ist]).get_sigma_s_bare(sim_crm[ism] - 1))
            print ""

    # sigma_s_GIPAW
    for ist in range(nsites):
        if wss == -1 or wss == sites_qe[ist]:
            site = sim_ss_run[0].get_atom(sites_qe[ist])
            print "sigma_s_GIPAW:", site.symbol+"("+str(sites_qe[ist])+")"
            print "%55s %15s" % ("nart", "art")
            for ism in range(nsim):
                if sim_type[ism] == "QE-GIPAW":
                   print "%10s %10s %15s %3d %15.2f %15.2f" % (sim_type[ism], sim_label[ism], sim_paw[ism], sim_crm[ism], 
                                      sim_ss_run[ism].get_atom(sites_qe[ist]).get_sigma_s_GIPAW_nart(sim_crm[ism] - 1),
                                      sim_ss_run[ism].get_atom(sites_qe[ist]).get_sigma_s_GIPAW_art(sim_crm[ism] - 1))
            print ""

    # sigma_s_core_relax
    for ist in range(nsites):
        if wss == -1 or wss == sites_qe[ist]:
            site = sim_ss_run[0].get_atom(sites_qe[ist])
            print "sigma_s_core_relax:", site.symbol+"("+str(sites_qe[ist])+")"
            print "%55s %15s" % ("nart", "art")
            for ism in range(nsim):
                if sim_type[ism] == "QE-GIPAW":
                   print "%10s %10s %15s %3d %15.2f %15.2f" % (sim_type[ism], sim_label[ism], sim_paw[ism], sim_crm[ism], 
                                      sim_ss_run[ism].get_atom(sites_qe[ist]).get_sigma_s_core_relax_nart(sim_crm[ism] - 1),
                                      sim_ss_run[ism].get_atom(sites_qe[ist]).get_sigma_s_core_relax_art(sim_crm[ism] - 1))
            print ""

    # f-sum rule
    print "f-sum rule:"
    print "%58s %15s %15s" % ("expected", "computed", "difference")
    for ism in range(nsim):
        if sim_type[ism] == "QE-GIPAW":
            #print "%10s %10s %15s %3d %15.2f %15.2f %15.2f" % (sim_type[ism], sim_label[ism], sim_paw[ism], sim_crm[ism], 
            #                                            sim_so_run[ism].get_f_sum_rule_exp(), 
            #                                            sim_so_run[ism].get_f_sum_rule_comp(), 
            #                                            sim_so_run[ism].get_f_sum_rule_exp()-sim_so_run[ism].get_f_sum_rule_comp())
            print " & %10s & %15.2f & %15.2f & %15.2f\\\\" % (sim_label[ism], 
                                                        sim_so_run[ism].get_f_sum_rule_exp(), 
                                                        sim_so_run[ism].get_f_sum_rule_comp(), 
                                                        sim_so_run[ism].get_f_sum_rule_exp()-sim_so_run[ism].get_f_sum_rule_comp())
elif wtd == 1:
    # For now, assuming that the first four are QE-GIPAW PBE, PW91, LDA, and PBEsol
    if wts == 0:
        for ists in range(hms): # how many sites to show
            sts_avg = sts[ists].split("-")
            for isq in range(len(sites_qe)):
                if int(sts_avg[0]) == sites_qe[isq]:
                    expi = isq
                    break
            linetp = "%15.2f " % Kiso_exp[expi]
            for ism in range(nsim): # all the simulations
                if sim_type[ism] == "QE-GIPAW":
                    sref = sim_ref_run[ism].get_atom(ref_id).get_sigma_o()
                    sigma_o = sim_so_run[ism].get_atom(int(sts_avg[0])).get_sigma_o()
                    Kiso = 0.0
                    for iavg in range(len(sts_avg)):
                        if dm == 0:
                            sigma_s = sim_ss_run[ism].get_atom(int(sts_avg[iavg])).get_sigma_s_nart(sim_crm[ism] - 1)
                        else:
                            sigma_s = sim_ss_run[ism].get_atom(int(sts_avg[iavg])).get_sigma_s_art(sim_crm[ism] - 1)
                        Kiso += (sref + ref_delta_exp) - (sigma_s + sigma_o)
                    linetp += "%15.2f " % (Kiso / len(sts_avg))
            print linetp
    elif wts == 1: # Cq
        for ists in range(hms): # how many sites to show
            for isq in range(len(sites_qe)):
                if int(sts[ists]) == sites_qe[isq]:
                    expi = isq
                    break
            linetp = "%15.2f " % Cq_exp[expi]
            for ism in range(nsim): # all the simulations
                #if sim_type[ism] == "QE-GIPAW":
                linetp += "%15.2f " % abs(sim_ss_run[ism].get_atom(int(sts[ists])).get_Cq())
            print linetp[0:len(linetp)-1]
    elif wts == 2:
        # For now, assuming that the next four are the Elk PBE, PW91, LDA, and PBEsol
        for ists in range(hms): # how many sites to show
            sts_avg = sts[ists].split("-")
            for ism in range(nsim): # all the simulations
                if sim_type[ism] == "QE-GIPAW" and sim_label[ism] == wxc:
                    if prefix == "VAl3":
                        ielk = int(sts_avg[0]) + 1
                    else:
                        ielk = int(sts_avg[0])
                    if dm == 0:
                        ssavg = 0.0
                        for iavg in range(len(sts_avg)):
                            ssavg += sim_ss_run[ism].get_atom(int(sts_avg[iavg])).get_sigma_s_nart(sim_crm[ism] - 1)
                        print "%15.2f %15.2f" % (sim_ss_run[ism+4].get_atom(ielk).get_sigma_s(), ssavg/len(sts_avg))
                    else:
                        ssavg = 0.0
                        for iavg in range(len(sts_avg)):
                            ssavg += sim_ss_run[ism].get_atom(int(sts_avg[iavg])).get_sigma_s_art(sim_crm[ism] - 1)
                        print "%15.2f %15.2f" % (sim_ss_run[ism+4].get_atom(ielk).get_sigma_s(), ssavg/len(sts_avg))
    elif wts == 3: # eta
        for ists in range(hms): # how many sites to show
            for isq in range(len(sites_qe)):
                if int(sts[ists]) == sites_qe[isq]:
                    expi = isq
                    break
            linetp = "%15.2f " % eta_exp[expi]
            for ism in range(nsim): # all the simulations
                #if sim_type[ism] == "QE-GIPAW":
                linetp += "%15.2f " % abs(sim_ss_run[ism].get_atom(int(sts[ists])).get_eta())
            print linetp[0:len(linetp)-1]
elif wtd == 2: # Show Tables with shifts
    last_site = -1
    show_prefix = True
    for ist in range(nsites):
        for ism in range(nsim):
            if last_site != ist:
                site = sim_ss_run[0].get_atom(sites_qe[ist])
                site_label = site.symbol+"("+str(sites_qe[ist])+")"
                site_Kiso_exp = Kiso_exp[ist]
                last_site = ist
            else:
                site_label = ""
                site_Kiso_exp = -100000000
            if sim_type[ism] == "QE-GIPAW":
                sigma_s_elk = 0.0
                for isme in range(nsim):
                    if sim_type[isme] == "Elk" and sim_label[ism] == sim_label[isme]:
                        sigma_s_elk = sim_ss_run[isme].get_atom(sites_elk[ist]).get_sigma_s()
                        break
                sref = sim_ref_run[ism].get_atom(ref_id).get_sigma_o()
                sigma_sn = sim_ss_run[ism].get_atom(sites_qe[ist]).get_sigma_s_nart(sim_crm[ism] - 1)
                sigma_sa = sim_ss_run[ism].get_atom(sites_qe[ist]).get_sigma_s_art(sim_crm[ism] - 1)
                sigma_o = sim_so_run[ism].get_atom(sites_qe[ist]).get_sigma_o()
                Kison = (sref + ref_delta_exp) - (sigma_sn + sigma_o)
                Kisoa = (sref + ref_delta_exp) - (sigma_sa + sigma_o)
                #print "%10s & %10s & %15.2f & %15.2f & %15.2f & %15.2f & %15.2f & %15.2f & %15.2f\\\\" % (
                #          site_label, sim_label[ism], sigma_s_elk, 
                #          sigma_sn, sigma_sa, sigma_o, Kison, Kisoa, site_Kiso_exp)-
                if show_prefix:
                    prefix_str = prefix
                    show_prefix = False
                else:
                    prefix_str = ""
                #print "%10s & %10s & %10s & %15.2f & %15.2f & %15.2f & %15.2f & %15.2f & %15.2f & %15.2f\\\\" % (prefix_str, 
                #          site_label, sim_label[ism], sigma_s_elk, 
                #          sigma_sn, sigma_sa, sigma_o, Kison, Kisoa, site_Kiso_exp)
                print "%10s & %10s & %10s & %15.2f & %15.2f & %15.2f & %15.2f & %15.2f & %15.2f\\\\" % (prefix_str, 
                          site_label, sim_label[ism], 
                          sigma_sn, sigma_sa, sigma_o, Kison, Kisoa, site_Kiso_exp)
elif wtd == 3: # Show Tables with EFG parameters
    last_site = -1
    show_prefix = True
    for ist in range(nsites):
        for ism in range(nsim):
            if last_site != ist:
                site = sim_ss_run[0].get_atom(sites_qe[ist])
                site_label = site.symbol+"("+str(sites_qe[ist])+")"
                site_Cq_exp = Cq_exp[ist]
                site_eta_exp = eta_exp[ist]
                last_site = ist
            else:
                site_label = ""
                site_Cq_exp = -100000000
                site_eta_exp = -100000000
            if sim_type[ism] == "QE-GIPAW":
                Cq_elk = 0.0
                eta_elk = 0.0
                for isme in range(nsim):
                    if sim_type[isme] == "Elk" and sim_label[ism] == sim_label[isme]:
                        Cq_elk = sim_ss_run[isme].get_atom(sites_elk[ist]).get_Cq()
                        eta_elk = sim_ss_run[isme].get_atom(sites_elk[ist]).get_eta()
                        break
                Cq = sim_ss_run[ism].get_atom(sites_qe[ist]).get_Cq()
                eta = sim_ss_run[ism].get_atom(sites_qe[ist]).get_eta()
                if show_prefix:
                    prefix_str = prefix
                    show_prefix = False
                else:
                    prefix_str = ""
                print "%10s & %10s & %10s & %15.2f & %15.2f & %15.2f & %15.2f & %15.2f & %15.2f\\\\" % (prefix_str, 
                          site_label, sim_label[ism], Cq_elk, 
                          Cq, abs(site_Cq_exp), eta_elk, eta, site_eta_exp)

