"""
Abstract Classes for describing a collection of SN to be simulated.
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

        The data structure for the snParams is not defined, and can be chosen
        appropriately in any concrete implementation. The basic requirement is
        that the property SN works with snParams
        """
        pass

    @abc.abstractmethod
    def sn(self, id):
        """
        Generator for supernova instances. Taking the snParams associated with
        a single SN, this uses a supernova model to produce all known
        information about the SN. In other words, combined with a set of
        observations (observation MJD, a description of Observational
        Conditions, and one can get model valued light curves. In order to add
        scatter due to sky noise, additional random number seeds are needed.
        """
        pass

    @abc.abstractmethod
    def modelFluxValue(self, id, mjd, bands):
        """
        Return a flux value in units of maggies for the object with id id
        at a time mjd, and in the band indicated by the string.

        Parameters
        ----------
        id : int or string, mandatory
            identifying index of the object
        mjd : `np.ndarray` of dtype float

        .. note : it is expected that this will be done by a method in the
            object instance
        """

        pass
    @abc.abstractmethod
    def lc(self, id, deviations=None, randomState=None):
        """
        returns an instance of a `LightCurve` for object corresponding to the id.
        """

class HomogeneousSNUniverse(with_metaclass(abc.ABCMeta, Universe)):

    @abc.abstractproperty
    def rates(self):
        """
        """
        pass
    
    @abc.abstractproperty
    def paramDistribution(self):
        pass

