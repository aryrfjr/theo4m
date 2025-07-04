"""
Definition of the GIPAWAtom class.

This module defines an extension to the the central object in 
the ASE package: the Atom object.
"""

from ase.atom import Atom

class GIPAWAtom(Atom):
    """
    GIPAWAtom object.

    The GIPAWAtom object inherits all the information about an atomic 
    site of a given system defined in the parent class ase.atom.Atom. 
    Additionally it has some extra information relevant in 
    ab-initio NMR simulations.

    Notes:

    The decompositions of the TOTAL Magnetic Shielding (ms) and the 
    EFG (efg) tensors are theory-dependent.
    
    Parameters:

    code: str 
        Code used to compute spectral parameters, can be:
        'QE-GIPAW' -> [Phys. Rev. B 63, 245101 (2001)], [Physical Review B 76, 024401 (2007)]
    """

    def __init__(self, symbol=None,
                 position=None, tag=None, 
                 momentum=None, mass=None,
                 magmom=None, charge=None):
        Atom.__init__(self, symbol,
                      position, tag, momentum, 
                      mass, magmom, charge)
        self.sigma_o = None
        self.sigma_para = None
        self.sigma_s_art = None # one for each core_relax_method
        self.sigma_s_nart = None
        self.rho_s_bare = None
        self.sigma_s_bare = None
        self.rho_s_GIPAW_art = None
        self.sigma_s_GIPAW_art = None
        self.rho_s_core_relax_art = None
        self.sigma_s_core_relax_art = None
        self.rho_s_GIPAW_nart = None
        self.sigma_s_GIPAW_nart = None
        self.rho_s_core_relax_nart = None
        self.sigma_s_core_relax_nart = None
        self.rho_s_total_art = None
        self.rho_s_total_nart = None
        self.rho_dip_xx_art = None
        self.rho_dip_yy_art = None
        self.rho_dip_zz_art = None
        self.rho_dip_xx_nart = None
        self.rho_dip_yy_nart = None
        self.rho_dip_zz_nart = None
        self.Kax_art = None
        self.Kax_nart = None
        self.efg = None # the total EFG tensor
        self.Vxx = None # eigenvalue
        self.Vyy = None 
        self.Vzz = None 
        self.vVxx = None # eigenvector
        self.vVyy = None 
        self.vVzz = None 
        self.Cq = None
        self.eta = None
        self.nuQ = None
        self.Bhf = None

    def set_sigma_o(self, so):
        self.sigma_o = so

    def set_sigma_para(self, sp):
        self.sigma_para = sp

    def set_sigma_s_art(self, ssa):
        self.sigma_s_art = ssa

    def set_sigma_s_nart(self, ssn):
        self.sigma_s_nart = ssn

    def set_rho_s_bare(self, sb):
        self.rho_s_bare = sb

    def set_sigma_s_bare(self, sb):
        self.sigma_s_bare = sb

    def set_rho_s_GIPAW_art(self, sga):
        self.rho_s_GIPAW_art = sga

    def set_sigma_s_GIPAW_art(self, sga):
        self.sigma_s_GIPAW_art = sga

    def set_rho_s_core_relax_art(self, sca):
        self.rho_s_core_relax_art = sca

    def set_sigma_s_core_relax_art(self, sca):
        self.sigma_s_core_relax_art = sca

    def set_rho_s_GIPAW_nart(self, sgn):
        self.rho_s_GIPAW_nart = sgn

    def set_sigma_s_GIPAW_nart(self, sgn):
        self.sigma_s_GIPAW_nart = sgn

    def set_rho_s_core_relax_nart(self, scn):
        self.rho_s_core_relax_nart = scn

    def set_sigma_s_core_relax_nart(self, scn):
        self.sigma_s_core_relax_nart = scn

    def set_rho_s_total_art(self, sta):
        self.rho_s_total_art = sta

    def set_sigma_s_total_art(self, sta):
        self.sigma_s_total_art = sta

    def set_rho_s_total_nart(self, stn):
        self.rho_s_total_nart = stn

    def set_rho_dip_xx_art(self, val):
        self.rho_dip_xx_art = val

    def set_rho_dip_yy_art(self, val):
        self.rho_dip_yy_art = val

    def set_rho_dip_zz_art(self, val):
        self.rho_dip_zz_art = val

    def set_rho_dip_xx_nart(self, val):
        self.rho_dip_xx_nart = val

    def set_rho_dip_yy_nart(self, val):
        self.rho_dip_yy_nart = val

    def set_rho_dip_zz_nart(self, val):
        self.rho_dip_zz_nart = val

    def set_Kax_art(self, val):
        self.Kax_art = val

    def set_Kax_nart(self, val):
        self.Kax_nart = val

    def set_efg(self, efg):
        self.efg = efg

    def set_Vxx(self, Vxx):
        self.Vxx = Vxx

    def set_Vyy(self, Vyy):
        self.Vyy = Vyy

    def set_Vzz(self, Vzz):
        self.Vzz = Vzz

    def set_vVxx(self, vVxx):
        self.vVxx = vVxx

    def set_vVyy(self, vVyy):
        self.vVyy = vVyy

    def set_vVzz(self, vVzz):
        self.vVzz = vVzz

    def set_Cq(self, Cq):
        self.Cq = Cq

    def set_eta(self, eta):
        self.eta = eta

    def set_nuQ(self, nuQ):
        self.nuQ = nuQ

    def set_Bhf(self, Bhf):
        self.Bhf = Bhf

    def get_sigma_o(self):
        return self.sigma_o

    def get_sigma_para(self):
        return self.sigma_para

    def get_sigma_s_art(self):
        return self.sigma_s_art

    def get_sigma_s_nart(self):
        return self.sigma_s_nart

    def get_rho_s_bare(self):
        return self.rho_s_bare

    def get_sigma_s_bare(self):
        return self.sigma_s_bare

    def get_rho_s_GIPAW_art(self):
        return self.rho_s_GIPAW_art

    def get_sigma_s_GIPAW_art(self):
        return self.sigma_s_GIPAW_art

    def get_rho_s_core_relax_art(self):
        return self.rho_s_core_relax_art

    def get_sigma_s_core_relax_art(self):
        return self.sigma_s_core_relax_art

    def get_rho_s_GIPAW_nart(self):
        return self.rho_s_GIPAW_nart

    def get_sigma_s_GIPAW_nart(self):
        return self.sigma_s_GIPAW_nart

    def get_rho_s_core_relax_nart(self):
        return self.rho_s_core_relax_nart

    def get_sigma_s_core_relax_nart(self):
        return self.sigma_s_core_relax_nart

    def get_rho_s_total_art(self):
        return self.rho_s_total_art

    def get_rho_s_total_nart(self):
        return self.rho_s_total_nart

    def get_rho_dip_xx_art(self):
        return self.rho_dip_xx_art

    def get_rho_dip_yy_art(self):
        return self.rho_dip_yy_art

    def get_rho_dip_zz_art(self):
        return self.rho_dip_zz_art

    def get_rho_dip_xx_nart(self):
        return self.rho_dip_xx_nart

    def get_rho_dip_yy_nart(self):
        return self.rho_dip_yy_nart

    def get_rho_dip_zz_nart(self):
        return self.rho_dip_zz_nart

    def get_Kax_art(self):
        return self.Kax_art

    def get_Kax_nart(self):
        return self.Kax_nart

    def get_efg(self):
        return self.efg

    def get_Vxx(self):
        return self.Vxx

    def get_Vyy(self):
        return self.Vyy

    def get_Vzz(self):
        return self.Vzz

    def get_vVxx(self):
        return self.vVxx

    def get_vVyy(self):
        return self.vVyy

    def get_vVzz(self):
        return self.vVzz

    def get_Cq(self):
        return self.Cq

    def get_eta(self):
        return self.eta

    def get_nuQ(self):
        return self.nuQ

    def get_Bhf(self):
        return self.Bhf

