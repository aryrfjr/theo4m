"""

Definition of the Interaction class.

"""

class Interaction(object):
    """

    Interaction object.

    The Cluster object contains information related to 
    an interaction between two atoms A and B.
    
    """

    def __init__(self):
        """ """
        self.id=-1
        self.symbA=None
        self.symbB=None
        self.distance=0.0
        self.transx=0
        self.transy=0
        self.transz=0
        self.info=0.0 # any information related to the interaction (i.e., -ICOHP)

    def set_id(self, param):
        """ """
        self.id = param

    def get_id(self):
        """ """
        return self.id

    def set_symbA(self, param):
        """ """
        self.symbA = param

    def get_symbA(self):
        """ """
        return self.symbA

    def set_symbB(self, param):
        """ """
        self.symbB = param

    def get_symbB(self):
        """ """
        return self.symbB

    def set_distance(self, param):
        """ """
        self.distance = param

    def get_distance(self):
        """ """
        return self.distance

    def set_transx(self, param):
        """ """
        self.transx = param

    def get_transx(self):
        """ """
        return self.transx

    def set_transy(self, param):
        """ """
        self.transy = param

    def get_transy(self):
        """ """
        return self.transy

    def set_transz(self, param):
        """ """
        self.transz = param

    def get_transz(self):
        """ """
        return self.transz

    def set_info(self, param):
        """ """
        self.info = param

    def get_info(self):
        """ """
        return self.info

