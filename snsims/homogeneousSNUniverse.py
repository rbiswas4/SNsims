"""
A concrete base class to represent supernovae such that
- their spatial location on the 2 sphere are uniformly distributed 
- the number density of SN (SN per comoving volume per observer frame time matches the number density predicted by a SN rate
- SN are represented by SALT parameters, and the parameters are drawn from a multivariate probability distribution without any dependence on enviroment
"""
from .universe import HomogeneousSNUniverse

class HomogeneousSN(HomogeneousSNUniverse):

    def __init__(self, surveyDuration):
        self.surveyDuration = surveyDuration
