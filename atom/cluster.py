"""

Definition of the Cluster class.

This module defines an additional set of information related to 
an atomic site associated to an Atom object from ASE package.

"""

import numpy as np
import math

class Cluster(object):
    """

    Cluster object.

    The Cluster object contains the geometrical information 
    associated to an atomic site and their neighbors in a given shell.
    
    Parameters:

    masses: list of float
        Atomic masses in atomic units.

    """

    def __init__(self):
        """ """
        self.ind=-1
        self.sind=-1
        self.symbol=None
        self.position=None
        self.ninds=[]
        self.nsinds=[]
        self.nsymbols=[]
        self.npositions=[]
        self.ndistances=[]

    def set_central_site(self, ind, sind, symbol, position):
        """ """
        self.ind = ind
        self.sind = sind
        self.symbol = symbol
        self.position = position

    def get_central_site_index(self):
        """ """
        return self.ind

    def get_central_site_secondary_index(self):
        """ """
        return self.sind

    def get_central_site_symbol(self):
        """ """
        return self.symbol

    def get_central_site_position(self):
        """ """
        return self.position

    def add_neighbor(self, ind, sind, symbol, position, distance):
        """ """
        self.ninds.append(ind)
        self.nsinds.append(sind)
        self.nsymbols.append(symbol)
        self.npositions.append(position)
        self.ndistances.append(distance)

    def get_neighbor_index(self, ind):
        """ """
        return self.ninds[ind]

    def get_neighbor_secondary_index(self, ind):
        """ """
        return self.nsinds[ind]

    def get_neighbors_indexes(self):
        """ """
        return self.ninds

    def get_neighbors_secondary_indexes(self):
        """ """
        return self.nsinds

    def get_neighbor_symbol(self, ind):
        """ """
        return self.nsymbols[ind]

    def get_neighbors_symbols(self):
        """ """
        return self.nsymbols

    def get_distance(self, ind):
        """ """
        return self.ndistances[ind]

    def get_number_of_neighbors(self):
        """ """
        return len(self.ninds)

    def to_xyz(self):
        print(len(self.nsymbols) + 1)
        print('')
        print(self.symbol, self.position[0], self.position[1], self.position[2])
        for i in range(len(self.nsymbols)):
            print(self.nsymbols[i], self.npositions[i][0], self.npositions[i][1], self.npositions[i][2])

