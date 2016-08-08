"""
Abstract Classes for describing a collection of SN to be simulated.
"""
from __future__ import absolute_import, print_function
from future.utils import with_metaclass
import abc

__all__ = ['Universe', 'HomogeneousSNUniverse', 'RateDistributions']

class Universe(with_metaclass(abc.ABCMeta, object)):
    """
    Class to represent a collection of SN. 

    .. notes: This does not need to know how the parameters or SN objects were
    obtained.
    """
    @abc.abstractproperty
    def randomState(self):
        """
        Random state characterizing the simulation
        """
        pass

    @abc.abstractproperty
    def snParams(self):
        """
        Generator of parameters for the SN involved. These parameters should be
        sufficient to completely describe the intrinsic properties of the SN.

        .. note : These parameters can be read off a file from some pre-existing
        samples, or be generated from a random number generator. In the latter
        case, we need to instantiate it with a random state characterizing the
        simulation. This is given by the abstract property randomState
        """
        pass

    @abc.abstractproperty
    def SN(self):
        """
        Generator for supernova instances. Taking the snParams associated with
        a single SN, this uses a supernova model to produce all known
        information about the SN. In other words, combined with a set of
        observations (observation MJD, a description of Observational
        Conditions, and one can get model valued light curves. In order to add
        scatter due to sky noise, additional random number seeds are needed.
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

class RateDistributions(with_metaclass(abc.ABCMeta, object)):

    @abc.abstractmethod
    def zSampleSize(self, zbinEdges, skyFraction):
        """
        Given a collection of edges of redshift bins, and a skyFraction
        return a collection of expected (possibly non-integral) numbers of SN.
        Since this is an expected number, the number is a float and it is
        perfectly fine to get 3.45 SN in a redshift bin.
        """
        pass
