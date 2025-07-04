"""

Definition of the AtomicSites class.

This module defines an extension to the the central object in 
the ASE package: the Atoms object.

"""

import numpy as np
from math import sqrt
from ase.atoms import Atom
from ase.atoms import Atoms
from ase.calculators.neighborlist import NeighborList
from theo4m.atom.cluster import Cluster

class AtomicSites(Atoms):
    """

    AtomicSites object.

    The AtomicSites object inherits all the information about the atoms 
    of a given system defined in the parent class ase.atoms.Atoms. 
    Additionally it also has some extra information treating the atoms 
    as crystallographic sites.

    """

    def __init__(self, symbols=None,
                 positions=None, numbers=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 scaled_positions=None,
                 cell=None, pbc=None, celldisp=None,
                 constraint=None,
                 calculator=None,
                 info=None):
        super(AtomicSites, self).__init__(symbols,
                 positions, numbers,
                 tags, momenta, masses,
                 magmoms, charges,
                 scaled_positions,
                 cell, pbc, celldisp,
                 constraint,
                 calculator,
                 info)
        self.local_arrays = {}
        self.local_arrays["ANONYMOUS_CLUSTERS"] = []
        self.nl = None

    def get_1st_coord(self, nb_cutoff, central_atom):
        #
        # https://github.com/qsnake/ase/blob/master/ase/calculators/emt.py
        #
        # using half of the distance set by user
        self.update_neighbor_list(nb_cutoff)
        neighbors, offsets = self.nl.get_neighbors(central_atom)
        offsets = np.dot(offsets, self.cell)
        # setting the cluster
        cl = Cluster()
        cl.set_central_site(central_atom, central_atom, \
                  self.__getitem__(central_atom).symbol, \
                  self.__getitem__(central_atom).position)
        for nbi, offset in zip(neighbors, offsets):
            d = self.__getitem__(nbi).position + offset - self.__getitem__(central_atom).position
            r = sqrt(np.dot(d, d))
            if r > 0 and r <= nb_cutoff and nbi != 0:
                nbpos = self.__getitem__(nbi).position + offset
                cl.add_neighbor(nbi, nbi, self.__getitem__(nbi).symbol, nbpos, r)
        return cl

    def update_neighbor_list(self, nb_cutoff):
        # the neighbors must be checked with PBC
        cr = []
        for i in range(self.get_number_of_atoms()):
            cr.append(nb_cutoff)
        self.nl = NeighborList(cr, bothways=True)
        self.nl.update(self)

    def find_clusters(self, name_clusters, central_specie, coord_species, r0, drm, drp):
        """Find clusters of atoms based on the passed parameters"""
        self.__find_clusters__(name_clusters, central_specie, coord_species, r0, drm, drp)

    def find_cluster(self, central_specie_id, coord_species, r0, drm, drp, secondary_id = False):
        """Find a cluster of atoms based on the passed parameters"""
        cl = self.__find_cluster__(central_specie_id, coord_species, r0, drm, drp, secondary_id)
        self.local_arrays["ANONYMOUS_CLUSTERS"].append(cl)
        return cl

    def get_clusters(self, name_clusters):
        """Returns a collection of clusters.
        """
        return self.local_arrays[name_clusters]

    def get_site_cluster(self, name_clusters, id_site, secondary_id = False):
        """Returns a cluster by its central site index.
        """
        cls = self.local_arrays[name_clusters]
        if secondary_id:
            for i in range(len(cls)):
                if cls[i].get_central_site_index() == id_site:
                    return cls[i]
        else:
            for i in range(len(cls)):
                if cls[i].get_central_site_secondary_index() == id_site:
                    return cls[i]
        return None

    def get_site_cluster_as_fragment(self, name_clusters, id_site, va=10.0, vb=10.0, vc=10.0, secondary_id = False):
        """ """
        cl = self.get_site_cluster(name_clusters, id_site, secondary_id)
        a=[]
        a.append(self.__getitem__(cl.get_central_site_index()))
        for i in range(cl.get_number_of_neighbors()):
            vs = cl.get_neighbor_symbol(i)
            vp = cl.get_central_site_position() + cl.get_distance(i)
            vt = cl.get_neighbor_secondary_index(i)
            a.append(Atom(vs,vp, tag = vt))
        f = Atoms(a, cell=(va,vb,vc))
        return f

    def __find_clusters__(self, name_clusters, central_specie, coord_species, r0, drm, drp):
        cls = [] # a collection of clusters with a name
        # Creating a collection of Cluster objects
        for i in range(self.get_number_of_atoms()):
            if self.__getitem__(i).symbol == central_specie: # central atom
                # add the current cluster to the collection
                cls.append(self.__find_cluster__(i, coord_species, r0, drm, drp))
        # Set the collection of Cluster objects as an array
        self.local_arrays[name_clusters] = cls

    def __find_cluster__(self, central_specie_id, coord_species, r0, drm, drp, secondary_id = False):
        # building a list of neighbors around the current atom
        if secondary_id:
            # for example, each atom in LAMMPS has an id
            for i in range(self.get_number_of_atoms()):
                if self.__getitem__(i).tag == central_specie_id:
                    cat = self.__getitem__(i) # central atom
                    cat_pos = i # the central atom position in the atom list
                    break
        else:
            cat = self.__getitem__(central_specie_id) # central atom
            cat_pos = central_specie_id # the central atom position in the atom list
        if self.nl == None:
            # Not using ase.neighborlist.NeighborList
            nbs = self.get_chemical_symbols()
            nbc = self.get_positions()
            nbsi = self.get_tags()
        else:
            # Using ase.neighborlist.NeighborList
            indices, offsets = self.nl.get_neighbors(central_specie_id)
            nbs = []
            nbc = []
            nbsi = []
            for inb, offset in zip(indices, offsets):
                nbs.append(self.__getitem__(inb).symbol)
                nbc.append(self.positions[inb] + np.dot(offset, self.get_cell()))
                nbsi.append(self.__getitem__(inb).tag)
        # setting the cluster
        cl = Cluster()
        cl.set_central_site(cat_pos, central_specie_id, cat.symbol, cat.position)
        for j in range(len(nbs)):
            if cat_pos != j and self.__is_coord__(nbs[j], coord_species):
                dist_v = nbc[j] - nbc[cat_pos]
                d = abs(np.linalg.norm(dist_v))
                if d <= r0 + drp and d >= r0 - drm:
                    cl.add_neighbor(j, nbsi[j], nbs[j], dist_v)
        return cl

    def __is_coord__(self, symbol, coord_species):
        """Says if a symbol in in the array coord_species
        """
        isc = False
        for i in range(len(coord_species)):
            if (coord_species[i] == symbol):
                isc = True
                break
        return isc

