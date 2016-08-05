"""
Class for describing a collection of SN
"""
from __future__ import absolute_import, print_function
from future.utils import with_metaclass
import abc

__all__ = ['Universe', 'HomogeneousSNUniverse']

class Universe(with_metaclass(abc.ABCMeta, object)):
    """
    Class to represent a collection of SN. 

    .. notes: This does not need to know how the parameters or SN objects were
    obtained.
    """
    @abc.abstractproperty
    def snParams(self):
        """
        """
        pass

    @abc.abstractproperty
    def SN(self):
        """
        """
        pass

class HomogeneousSNUniverse(with_metaclass(abc.ABCMeta, Universe)):

    @abc.abstractproperty
    def rates(self):
        """
        """
        pass
    
    @abc.abstractproperty
    def paramDistribution(self):
        pass


class RateDistributions(with_metaclass(abc.ABCMeta, Object)):

    @abc.abstractmethod
    def zSamples(self, zbinEdges, skyFraction):
        """
        Given a collection of edges of redshift bins, and a skyFraction
        return a collection of expected (possibly non-integral) numbers of SN.
        Since this is an expected number, the number is a float and it is
        perfectly fine to get 3.45 SN in a redshift bin.
        """
        pass

