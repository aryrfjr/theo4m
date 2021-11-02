"""
Definition of the LAPWAtom class.

This module defines an extension to the the central object in 
the ASE package: the Atom object.
"""

from ase.atom import Atom

class LAPWAtom(Atom):
    """
    LAPWAtom object.

    The LAPWAtom object inherits all the information about an atomic 
    site of a given system defined in the parent class ase.atom.Atom. 
    Additionally it has some extra information relevant in 
    ab-initio NMR simulations.

    Notes:

    The decompositions of the TOTAL Magnetic Shielding (ms) and the 
    EFG (efg) tensors are theory-dependent.
    
    Parameters:

    """

    def __init__(self, symbol=None,
                 position=None, tag=None, 
                 momentum=None, mass=None,
                 magmom=None, charge=None):
        Atom.__init__(self, symbol,
                      position, tag, momentum, 
                      mass, magmom, charge)
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
        self.sigma_s = None

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

    def set_sigma_s(self, sigma_s):
        self.sigma_s = sigma_s

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

    def get_sigma_s(self):
        return self.sigma_s

