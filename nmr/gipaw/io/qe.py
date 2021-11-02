# -*- coding: utf-8 -*-
"""
Reads the Quantum Espresso and QE-GIPAW output files.
"""

# http://docs.scipy.org/doc/numpy/
import numpy as np
from numpy import array, dot
# from this same package
from theo4m.nmr.gipaw.gipawrun import GIPAWRun
from theo4m.nmr.gipaw.gipawatom import GIPAWAtom

__all__ = [
        'read_qe',
        ]

def read_qe(rundir, prefix, tasks, nuclide_label = "", core_relax = ""):
    """
    Reads the pw.x output file and the gipaw.x output files for 
    the 'nmr', 'hyperfine' (art and nart), and 'efg' tasks, 
    returning a GIPAWRun object.
    """
    MHz2dens = 1165.627050372 # conversion factor
    pwfilei = rundir + '/'+prefix+'.pw.inp'
    pwfileo = rundir + '/'+prefix+'.pw.out'
    nmrfile = rundir + '/'+prefix+'-nmr.gipaw.out'
    # the HF and EFG tasks file names may vary
    hfartfile = rundir + '/'+prefix
    hfnartfile = rundir + '/'+prefix
    efgfile = rundir + '/'+prefix
    if nuclide_label != "":
        hfartfile += '-'+nuclide_label
        hfnartfile += '-'+nuclide_label
        efgfile += '-'+nuclide_label
    hfartfile += '-hfi-art'
    hfnartfile += '-hfi-nart'
    #print core_relax
    if core_relax != "":
        hfartfile += '-'+core_relax
        hfnartfile += '-'+core_relax
    hfartfile += '.gipaw.out'
    hfnartfile += '.gipaw.out'
    efgfile += '-efg.gipaw.out'
    bndfile = rundir + '/'+prefix+'-bands.dat'
    gipaw = GIPAWRun()
    # parsing the SCF input file
    fileobj = open(pwfilei)
    lines = fileobj.readlines()
    fileobj.close()
    for i, line in enumerate(lines):
        if 'ntyp' in line:
        # Reading the number of atomic species
            tmp_s = line.split()[2]
            if ',' in tmp_s: tmp_s = tmp_s[0:len(tmp_s)-1]
            ntyp = int(tmp_s)
        if 'tot_magnetization' in line:
        # Reading the total magnetization
            tmp_s = line.split()[2]
            if ',' in tmp_s: tmp_s = tmp_s[0:len(tmp_s)-1]
            gipaw.set_total_magnetization(float(tmp_s))
    # parsing the SCF output file
    fileobj = open(pwfileo)
    lines = fileobj.readlines()
    fileobj.close()
    for i, line in enumerate(lines):
        if 'number of atoms/cell' in line:
        # Reading the number of atoms in the cell
            nat = int(line.split()[4])
        elif 'crystal axes: (cart. coord. in units of alat)' in line:
        # Reading the cell vectors
            cell = [x.split()[3:6] for x in lines[i + 1:i + 4]]
            cell = array([[float(col) for col in row] for row in cell])
            gipaw.set_cell(cell)
        elif 'Crystallographic axes' in line:
        # Reading the atoms and creating a collection
            geom_start = i + 3
            geom_stop = geom_start + nat
            symbols = [line.split()[1] for line in lines[geom_start:geom_stop]]
            geom = dot(array([[float(col) for col in line.split()[6:9]]
                for line in lines[geom_start:geom_stop]]), cell)
            atoms = [None]*len(symbols)
            for ia in range(len(symbols)):
                atom = GIPAWAtom(symbol = symbols[ia], position = geom[ia])
                atoms[ia] = atom
                atoms[ia].tag = ia + 1
            gipaw.set_atoms(atoms)
        elif 'number of k points=' in line:
        # Reading the number of k-poins
            gipaw.set_k_points_number(int(line.split()[4]))
        elif 'number of electrons' in line:
        # Reading the number electrons
            gipaw.set_nelec(float(line.split()[4]))
        elif 'the spin up/dw Fermi energies are' in line:
        # Reading the Fermi energies
            gipaw.set_fermi_energies(float(line.split()[6]),float(line.split()[7]))
        elif 'the Fermi energy is' in line:
        # Reading the Fermi energy
            gipaw.set_fermi_energy(float(line.split()[4]))
        elif line.startswith('!'):
        # Reading the total energy
            gipaw.set_total_energy(float(line.split()[4]))
    # parsing the 'hyperfine' task output for art
    if 'hyperfine_art' in tasks:
        fileobj = open(hfartfile)
        lines = fileobj.readlines()
        fileobj.close()
        for i, line in enumerate(lines):
            if '----- spin-densities in bohrradius^-3 -----' in line:
                for ia in range(nat):
                    atoms[ia].set_rho_s_bare(float(lines[i+2+ia].split()[2]))
                    atoms[ia].set_rho_s_GIPAW_art(float(lines[i+2+ia].split()[3]))
                    atoms[ia].set_rho_s_core_relax_art(float(lines[i+2+ia].split()[4]))
                    atoms[ia].set_rho_s_total_art(float(lines[i+2+ia].split()[5]))
            if 'PRINCIPAL AXIS OF THE DIPOLAR COUPLINGS:' in line:
                il = i
                for ia in range(nat):
                    atoms[ia].set_rho_dip_xx_art(float(lines[il+1].split()[3])/MHz2dens)
                    atoms[ia].set_rho_dip_yy_art(float(lines[il+2].split()[3])/MHz2dens)
                    atoms[ia].set_rho_dip_zz_art(float(lines[il+3].split()[3])/MHz2dens)
                    il += 4
    # parsing the 'hyperfine' task output for nart
    if 'hyperfine_nart' in tasks or 'hyperfine' in tasks:
        fileobj = open(hfnartfile)
        lines = fileobj.readlines()
        fileobj.close()
        for i, line in enumerate(lines):
            #if 'the spin majority channel is' in line:
            #    print line.split()[5]
            if '----- spin-densities in bohrradius^-3 -----' in line:
                for ia in range(nat):
                    atoms[ia].set_rho_s_bare(float(lines[i+2+ia].split()[2]))
                    atoms[ia].set_rho_s_GIPAW_nart(float(lines[i+2+ia].split()[3]))
                    atoms[ia].set_rho_s_core_relax_nart(float(lines[i+2+ia].split()[4]))
                    atoms[ia].set_rho_s_total_nart(float(lines[i+2+ia].split()[5]))
            if 'PRINCIPAL AXIS OF THE DIPOLAR COUPLINGS:' in line:
                il = i
                for ia in range(nat):
                    atoms[ia].set_rho_dip_xx_nart(float(lines[il+1].split()[3])/MHz2dens)
                    atoms[ia].set_rho_dip_yy_nart(float(lines[il+2].split()[3])/MHz2dens)
                    atoms[ia].set_rho_dip_zz_nart(float(lines[il+3].split()[3])/MHz2dens)
                    il += 4
    # parsing the 'efg' task
    if 'efg' in tasks:
        fileobj = open(efgfile)
        lines = fileobj.readlines()
        fileobj.close()
        for i, line in enumerate(lines):
            if '----- total EFG (symmetrized) -----' in line:
                il = i + 1
                for ia in range(nat):
                    efg = array([[float(col) for col in line.split()[2:5]] for line in lines[il:il + 3]])
                    atoms[ia].set_efg(efg)
                    il = il + 4
    # parsing the 'nmr' task
    if 'nmr' in tasks:
        fileobj = open(nmrfile)
        lines = fileobj.readlines()
        fileobj.close()
        for i, line in enumerate(lines):
            if 'f-sum rule (should be' in line:
                vfsr = line.split()[4]
                gipaw.set_f_sum_rule_exp(float(vfsr[0:len(vfsr)-2]))
                gipaw.set_f_sum_rule_comp((float(lines[i+1].split()[0]) + 
                                               float(lines[i+2].split()[1]) + 
                                               float(lines[i+3].split()[2])) / 3.0)
            if 'Total NMR chemical shifts in ppm: ---------------------------------------' in line:
                il = i + 3
                for ia in range(nat):
                    atoms[ia].set_sigma_o(float(lines[il].split()[10]))
                    il = il + 10
    # parsing the 'bands' task
    if 'bands' in tasks:
        fileobj = open(bndfile)
        lines = fileobj.readlines()
        fileobj.close()
        strt = lines[0].split()[2]
        nbnd = int(strt[0:len(strt)-1]) # removing the comma
        nks = int(lines[0].split()[4])
        kpts = [0.0,0.0,0.0]*nks
        bandst = []
        nlbs = int(nbnd/10) + 1
        ikpt = 0
        il = 1
        while il < len(lines):
            kpts[ikpt] = [float(col) for col in lines[il].split()]
            il += 1
            bandst.append([0.0]*nbnd)
            for j in range(nlbs):
                if j < nlbs - 1:
                    bandst[ikpt][j*10:(j+1)*10] = [float(col) - gipaw.get_fermi_energy() for col in lines[il].split()]
                    #bandst[ikpt][j*10:(j+1)*10] = [float(col) for col in lines[il].split()]
                else:
                    bandst[ikpt][j*10:(j*10)+(nbnd%10)] = [float(col) - gipaw.get_fermi_energy() for col in lines[il].split()]
                    #bandst[ikpt][j*10:(j*10)+(nbnd%10)] = [float(col) for col in lines[il].split()]
                il += 1
            ikpt += 1
        bands = []
        for ib in range(nbnd):
            bands.append([0.0]*nks)
            for ik in range(nks):            
                bands[ib][ik] = bandst[ik][ib]
        gipaw.set_kptsb(kpts)
        gipaw.set_bands(bands)
    return gipaw

