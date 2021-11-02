"""

Definition of the SOAPList class.

"""
# https://libatoms.github.io/GAP/quippy-descriptor-tutorial.html#A-many-body-descriptor:-SOAP
import quippy
import ase
from quippy import descriptors

class SOAPList(object):

    """

    SOAPList object.

    A list of SOAP descriptor vectors (the power 
    spectrum) associated to a collection of ase.Atoms 
    objects. It can manage a set of frames (supercells) 
    as well as the deriving set of clusters 
    representing atomic local environments.
    
    Parameters:

    atoms: ase.Atoms object

    """
    def __init__(self, atoms, verb = False):
        self.verb = verb
        # loading frames (cells)
        self.frames = atoms
        if self.verb: print ("Computing SOAPs for an ase.Atoms object with %d cells (frames)." % (len(self.frames)))
        self.asoaps = []
        self.asoapsargs = []
        self.pasoaps = []
        self.pacats = []
        self.pasoapsargs = []

    def get_frames(self):
        return self.frames

    """ 

    This method computes the average SOAP vectors for a set 
    of frames. Adding them to a set of descriptors computed 
    with different arguments (array asoapsargs).

    """
    def compute_average(self, 
                       cutoff = "1.0", 
                       l_max = "4",
                       n_max = "4",
                       atom_sigma = "0.5",
                       n_Z = "1", 
                       Z = "{1}", 
                       n_species = "1", 
                       species_Z = "{1}"):
        # setting the options of the SOAP descriptor:
        # https://libatoms.github.io/QUIP/descriptors.html#quippy.descriptors.Soap
        # see also:
        # http://libatoms.github.io/QUIP/Tutorials/Introduction.html
        # http://libatoms.github.io/QUIP/Tutorials/quippy-descriptor-tutorial.html#A-many-body-descriptor:-SOAP
        argsd = ("soap cutoff="+cutoff+
                " l_max="+l_max+
                " n_max="+n_max+
                " atom_sigma="+atom_sigma+
                " n_Z="+n_Z+
                " Z="+Z+
                " n_species="+n_species+
                " species_Z="+species_Z+" average=True")
        if self.verb: print (argsd) # verbosity
        # now creating a collection of SOAP vectors for each frame
        desc = descriptors.Descriptor(argsd)
        soaps = []
        desc_calcs = desc.calc(self.frames, grad=False)
        if type(desc_calcs) is list:
            for incf in range(len(self.frames)):
                # recalling that the index [0] above is because
                # there is a single average SOAP vector per frame
                soaps.append(desc_calcs[incf]["data"][0])
        else: # it's a dictionary
            # recalling that the index [0] above is because
            # there is a single average SOAP vector per frame
            soaps.append(desc_calcs["data"][0])
        self.asoaps.append(soaps)
        self.asoapsargs.append(argsd)

    """ 

    This method computes the per-atom SOAP vectors for a set 
    of frames. Adding them to a set of descriptors computed 
    with different arguments (array pasoapsargs).

    """
    def compute_per_atom(self, 
                       cutoff = "1.0", 
                       l_max = "4",
                       n_max = "4",
                       atom_sigma = "0.5",
                       n_Z = "1", 
                       Z = "{1}", 
                       n_species = "1", 
                       species_Z = "{1}"):
        # setting the options of the SOAP descriptor:
        # https://libatoms.github.io/QUIP/descriptors.html#quippy.descriptors.Soap
        # see also:
        # http://libatoms.github.io/QUIP/Tutorials/Introduction.html
        # http://libatoms.github.io/QUIP/Tutorials/quippy-descriptor-tutorial.html#A-many-body-descriptor:-SOAP
        argsd = ("soap cutoff="+cutoff+
                " l_max="+l_max+
                " n_max="+n_max+
                " atom_sigma="+atom_sigma+
                " n_Z="+n_Z+
                " Z="+Z+
                " n_species="+n_species+
                " species_Z="+species_Z+" average=False")
        if self.verb: print (argsd) # verbosity
        # now creating a collection of SOAP vectors for each frame
        desc = descriptors.Descriptor(argsd)
        soaps = []
        cats = []
        desc_calcs = desc.calc(self.frames, grad=False)
        if type(desc_calcs) is list:
            for incf in range(len(self.frames)):
                # the SOAP vectors
                soaps.append(desc_calcs[incf]["data"])
                # the respective central atoms
                cats.append(desc_calcs[incf]["ci"])
        else: # it's a dictionary
            # the SOAP vectors
            soaps.append(desc_calcs["data"])
            # the respective central atoms
            cats.append(desc_calcs["ci"])
        self.pasoaps.append(soaps)
        self.pacats.append(cats)
        self.pasoapsargs.append(argsd)

    """

    Return the average SOAP vectors of a given instance based on 
    the arguments (arg_index) of a certain frame.

    """
    def get_average_soap(self, arg_index, frame_index):
        return self.asoaps[arg_index][frame_index]

    """

    Return the string with arguments of a average SOAP (arg_index).

    """
    def get_asoaps_args(self, arg_index):
        return self.asoapsargs[arg_index]

    """

    Return the per-atom SOAP vectors of a given instance based on 
    the arguments (arg_index) of a certain frame. The second returned
    array are the corresponding central atoms (descriptor_index_0based).

    """
    def get_per_atom_soaps(self, arg_index, frame_index):
        return self.pasoaps[arg_index][frame_index], self.pacats[arg_index][frame_index]

    """

    Return the string with arguments of a per-atom SOAP (arg_index).

    """
    def get_pasoaps_args(self, arg_index):
        return self.pasoapsargs[arg_index]

    """

    Return the number of frames.

    """
    def get_n_frames(self):
        return len(self.frames)

