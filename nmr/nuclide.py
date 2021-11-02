"""
An atomic specie characterized by the specific constitution of its nucleus.

https://en.wikipedia.org/wiki/Nuclide
"""

class Nuclide(object):
    """
    Nuclide object.
    """

    def __init__(self):
        self.symbol = None
        self.ngf = None
        self.qm = None
        self.ns = None
        self.mn = None

    def set_symbol(self, s):
        self.symbol = s

    def get_symbol(self):
        return self.symbol

    def set_quadrupolar_moment(self, qm):
        self.qm = qm

    def get_quadrupolar_moment(self):
        return self.qm

    def set_nuclear_spin(self, ns):
        self.ns = ns

    def get_nuclear_spin(self):
        return self.ns

    def set_mass_number(self, mn):
        self.mn = mn

    def get_mass_number(self):
        return self.mn

    def set_gyromangnetic_ratio(self, gr):
        self.gr = gr
        # See TUM KS diaries on 18/02/2016-(10)
        self.ngf = gr * 0.20879370371676311495

    def get_gyromangnetic_ratio(self):
        return self.gr

    def set_nuclear_g_factor(self, g):
        self.ngf = g

    def get_nuclear_g_factor(self):
        return self.ngf

