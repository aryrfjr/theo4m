"""
Definition of the LAPWRun class.

Contains information related to a LAPW simulation.
"""

from numpy import linalg as la

class LAPWRun(object):
    """
    LAPWRun object.
    """

    def __init__(self):
        self.atoms = None
        self.nuclides = None
        self.total_energy = None
        self.nkpts = None
        self.cell = None
        self.Bext = None
        self.kptsb = None
        self.bands = None
        self.ncore = None
        self.nval = None
        self.nempty = None
        self.spinpol = None

    def set_atoms(self, a):
        self.atoms = a

    def set_nuclides(self, n):
        self.nuclides = n

    def set_total_energy(self, te):
        self.total_energy = te

    def set_k_points_number(self, n):
        self.nkpts = n

    def set_cell(self, c):
        self.cell = c

    def set_Bext(self, b):
        self.Bext = b

    def set_kptsb(self, kptsb):
        self.kptsb = kptsb

    def set_bands(self, bands):
        self.bands = bands

    def set_ncore(self, ncore):
        self.ncore = ncore

    def set_nval(self, nval):
        self.nval = nval

    def set_nempty(self, nempty):
        self.nempty = nempty

    def set_spinpol(self, spinpol):
        self.spinpol = spinpol

    def get_atoms(self):
        return self.atoms

    def get_atom(self, aid):
        return self.atoms[aid]

    def get_nuclides(self):
        return self.nuclides

    def get_total_energy(self):
        return self.total_energy

    def get_k_points_number(self):
        return self.nkpts

    def get_cell(self):
        return self.cell

    def get_Bext(self):
        return self.Bext

    def get_kptsb(self):
        return self.kptsb

    def get_bands(self):
        return self.bands

    def get_ncore(self):
        return self.ncore

    def get_nval(self):
        return self.nval

    def get_nempty(self):
        return self.nempty

    def get_spinpol(self):
        return self.spinpol

    def get_atoms_by_symbol(self, symbol):
        absymb = []
        for ia in range(len(self.atoms)):
            if self.atoms[ia].symbol == symbol:
                absymb.append(self.atoms[ia])
        return absymb

    def compute_sigma_s(self):
        for atom in self.atoms:
            atom.set_sigma_s(-(atom.get_Bhf() / self.Bext) * 1000000.0)

    def compute_efg_parameters(self):
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
            nuclide = self.__get_nuclide__(atom.symbol)
            atom.set_Cq(atom.get_Vzz() * nuclide.get_quadrupolar_moment() * cte / 100.0)
            if abs(atom.get_Vzz()) > 0.0:
                atom.set_eta((atom.get_Vxx() - atom.get_Vyy()) / atom.get_Vzz())
            else:
                atom.set_eta(0.0)
            # The quadrupolar frequency also in MHz
            # http://anorganik.uni-tuebingen.de/klaus/nmr/index.php?p=conventions/efg/quadtools
            twoi = 2.0 * nuclide.get_nuclear_spin()
            #atom.set_nuQ((3.0 * atom.get_Cq()) / (twoi * (twoi - 1.0)))

    def __get_nuclide__(self, symb):
        for nuclide in self.nuclides:
            if nuclide.get_symbol() == symb:
                return nuclide

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

