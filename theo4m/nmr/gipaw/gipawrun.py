"""
Definition of the GIPAWRun class.

Contains information related to a GIPAW simulation.
"""

from numpy import linalg as la

class GIPAWRun(object):
    """
    GIPAWRun object.
    """

    def __init__(self):
        self.atoms = None
        self.total_energy = None
        self.fermi_up = None
        self.fermi_down = None
        self.nkpts = None
        self.tmag = None
        self.cell = None
        self.Bext = None
        self.kptsb = None
        self.bands = None
        self.nelec = None
        self.f_sum_rule_exp = None
        self.f_sum_rule_comp = None

    def set_atoms(self, a):
        self.atoms = a

    def set_total_energy(self, te):
        self.total_energy = te

    def set_fermi_energies(self, up, down):
        self.fermi_up = up
        self.fermi_down = down

    # for non-spin-polarized runs
    def set_fermi_energy(self, f):
        self.fermi_up = f

    def set_k_points_number(self, n):
        self.nkpts = n

    def set_total_magnetization(self, tm):
        self.tmag = tm

    def set_cell(self, c):
        self.cell = c

    def set_kptsb(self, kptsb):
        self.kptsb = kptsb

    def set_bands(self, bands):
        self.bands = bands

    def set_nelec(self, nelec):
        self.nelec = nelec

    def set_f_sum_rule_exp(self, fsre):
        self.f_sum_rule_exp = fsre

    def set_f_sum_rule_comp(self, fsrc):
        self.f_sum_rule_comp = fsrc

    def get_atoms(self):
        return self.atoms

    def get_atom(self, aid):
        return self.atoms[aid]

    def get_total_energy(self):
        return self.total_energy

    def get_fermi_energies(self):
        return self.fermi_up, self.fermi_down

    # for non-spin-polarized runs
    def get_fermi_energy(self):
        return self.fermi_up

    def get_k_points_number(self):
        return self.nkpts

    def get_total_magnetization(self):
        return self.tmag

    def get_cell(self):
        return self.cell

    def get_Bext(self):
        return self.Bext

    def get_kptsb(self):
        return self.kptsb

    def get_bands(self):
        return self.bands

    def get_nelec(self):
        return self.nelec

    def get_f_sum_rule_exp(self):
        return self.f_sum_rule_exp

    def get_f_sum_rule_comp(self):
        return self.f_sum_rule_comp

    def compute_sigma_s(self, art, nart):
        Bmag = 9.274009682e-24 # Bohr magneton in J/T
        eV2J = 1.60217733e-19
        eV2Ry = 0.073498618
        fpi = 12.5663708 # 4*pi
        mu0 = fpi*0.0000001 # permeability of a vacuum in H/m
        a0 = 0.0000000000529177 # Bohr radius in m
        ooa03 = 1.0/(a0**3.0) # 1/(a0^3)
        ge = 2.002319 # Electron spin g-factor
        self.Bext = ((self.fermi_up - self.fermi_down)*eV2J)/(2.0*Bmag)
        for atom in self.atoms:
            Bhf = (2.0/3.0)*0.5*mu0*ooa03*ge*Bmag*atom.get_rho_s_bare()*1000.0
            atom.set_sigma_s_bare(-Bhf/self.Bext*1000.0)
            if art:
                Bhf = (2.0/3.0)*0.5*mu0*ooa03*ge*Bmag*atom.get_rho_s_GIPAW_art()*1000.0
                atom.set_sigma_s_GIPAW_art(-Bhf/self.Bext*1000.0)
                Bhf = (2.0/3.0)*0.5*mu0*ooa03*ge*Bmag*atom.get_rho_s_core_relax_art()*1000.0
                atom.set_sigma_s_core_relax_art(-Bhf/self.Bext*1000.0)
                Bhf = (2.0/3.0)*0.5*mu0*ooa03*ge*Bmag*(atom.get_rho_s_total_art())*1000.0
                atom.set_sigma_s_art(-Bhf/self.Bext*1000.0)
                atom.set_Bhf(Bhf/1000.0)
                Kax = atom.get_rho_dip_zz_art() - atom.get_rho_dip_xx_art()
                Bhf = (2.0/3.0)*0.5*mu0*ooa03*ge*Bmag*(Kax)*1000.0
                atom.set_Kax_art(abs(-Bhf/self.Bext*1000.0))
            if nart:
                Bhf = (2.0/3.0)*0.5*mu0*ooa03*ge*Bmag*atom.get_rho_s_GIPAW_nart()*1000.0
                atom.set_sigma_s_GIPAW_nart(-Bhf/self.Bext*1000.0)
                Bhf = (2.0/3.0)*0.5*mu0*ooa03*ge*Bmag*atom.get_rho_s_core_relax_nart()*1000.0
                atom.set_sigma_s_core_relax_nart(-Bhf/self.Bext*1000.0)
                Bhf = (2.0/3.0)*0.5*mu0*ooa03*ge*Bmag*atom.get_rho_s_total_nart()*1000.0
                atom.set_sigma_s_nart(-Bhf/self.Bext*1000.0)
                atom.set_Bhf(Bhf/1000.0)
                Kax = atom.get_rho_dip_zz_nart() - atom.get_rho_dip_xx_nart()
                Bhf = (2.0/3.0)*0.5*mu0*ooa03*ge*Bmag*(Kax)*1000.0
                atom.set_Kax_nart(abs(-Bhf/self.Bext*1000.0))

    def get_atoms_by_symbol(self, symbol):
        absymb = []
        for ia in range(len(self.atoms)):
            if self.atoms[ia].symbol == symbol:
                absymb.append(self.atoms[ia])
        return absymb

    def compute_efg_parameters(self, nuclide):
        # Constants
        HARTREE_SI = 4.35974394E-18 # J
        ELECTRONVOLT_SI = 1.602176487E-19 # J  
        AUTOEV = HARTREE_SI / ELECTRONVOLT_SI
        RYTOEV = AUTOEV / 2.0
        BOHR_RADIUS_SI = 0.52917720859E-10 # m
        BOHR_RADIUS_CM = BOHR_RADIUS_SI * 100.0
        BOHR_RADIUS_ANGS = BOHR_RADIUS_CM * 1.0E8
        ANGSTROM_AU = 1.0/BOHR_RADIUS_ANGS
        ELECTRON_SI = 1.602176487E-19 # C
        for atom in self.atoms:
            tp = self.__diag_tensor__(atom.get_efg(), sort=True)
            atom.set_Vxx(tp[0][1]) # the EFG Vxx eigenvalue
            atom.set_Vyy(tp[1][1]) # the EFG Vyy eigenvalue
            atom.set_Vzz(tp[2][1]) # the EFG Vzz eigenvalue
            atom.set_vVxx(tp[0][2:5]) # the EFG Vxx eigenvector
            atom.set_vVyy(tp[1][2:5]) # the EFG Vyy eigenvector
            atom.set_vVzz(tp[2][2:5]) # the EFG Vzz eigenvector
            cte = (RYTOEV * 2.0 * (ANGSTROM_AU ** 2) * ELECTRONVOLT_SI * 1E20) / 6.62620
            atom.set_Cq(atom.get_Vzz() * nuclide.get_quadrupolar_moment() * cte / 100.0)
            if abs(atom.get_Vzz()) > 0.0:
                atom.set_eta((atom.get_Vxx() - atom.get_Vyy()) / atom.get_Vzz())
            else:
                atom.set_eta(0.0)
            # The quadrupolar frequency also in MHz
            # http://anorganik.uni-tuebingen.de/klaus/nmr/index.php?p=conventions/efg/quadtools
            twoi = 2.0 * nuclide.get_nuclear_spin()
            atom.set_nuQ((3.0 * atom.get_Cq()) / (twoi * (twoi - 1.0)))

    def __diag_tensor__(self, tns, sort=False):
        """Compute a 3x3 tensor eigenvalues and eigenvectors and sort it:"""
        sol = la.eig(tns)
        evl = sol[0]
        evc = sol[1]
        ind = [0,1,2]
        # Sorting: |Vzz| > |Vyy| > |Vxx|
        t = [(abs(evl[0]),evl[0],evc[0,0],evc[1,0],evc[2,0],ind[0]),
             (abs(evl[1]),evl[1],evc[0,1],evc[1,1],evc[2,1],ind[1]),
             (abs(evl[2]),evl[2],evc[0,2],evc[1,2],evc[2,2],ind[2])]
        if sort: t = sorted(t, key=lambda t: t[0])
        return t

