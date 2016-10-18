"""
A module with the abstract classes for implementing populaions of SN
"""
from __future__ import absolute_import, print_function, division
from future.utils import with_metaclass
import abc
import numpy as np
from .tessellations import Tiling


__all__ = ['SNParamDistribution', 'RateDistributions', 'SALT2Parameters',
           'PositionSamples', 'DisjointParams']

class SNParamDistribution(with_metaclass(abc.ABCMeta, object)):
    """

    """
    @abc.abstractproperty
    def randomState(self):
        pass

    @abc.abstractmethod
    def get_randomState(self):
        pass

    @abc.abstractproperty
    def set_randomState(self):
        pass


    @abc.abstractproperty
    def sn_params(self):
        pass

class RateDistributions(with_metaclass(abc.ABCMeta, object)):

    @abc.abstractproperty
    def randomState(self):
        pass

    @abc.abstractmethod
    def zSampleSize(self):
        """
        Given a collection of edges of redshift bins, and a skyFraction
        return a collection of expected (possibly non-integral) numbers of SN.
        Since this is an expected number, the number is a float and it is
        perfectly fine to get 3.45 SN in a redshift bin.
        """
        pass

    @abc.abstractproperty
    def zSamples(self):
        pass

class SALT2Parameters(with_metaclass(abc.ABCMeta, object)):

    @abc.abstractmethod
    def __init__(self, numSN, zSamples, snids=None, cSigma=0.1, x1Sigma=1.0,
                 rng=None):
        
        self._numSN = numSN
        pass

    @property
    def numSN(self):
        return self.numSN

    @abc.abstractproperty
    def randomState(self):
        pass

    @abc.abstractproperty
    def paramSamples(self):
        pass

class PositionSamples(with_metaclass(abc.ABCMeta, Tiling)):

    @abc.abstractproperty
    def randomState(self):
        pass

    def numPositions(self, tileID):
        return len(self.positions[0])

    @abc.abstractproperty
    def positions(self, tileID):
        # return ra, dec
        pass
    @staticmethod
    def samplePatchOnSphere(phi, theta, delta, size, rng):
        """
        Uniformly distributes samples on a patch on a sphere between phi \pm delta,
        and theta \pm delta on a sphere. Uniform distribution implies that the
        number of points in a patch of sphere is proportional to the area of the
        patch. Here, the coordinate system is the usual
        spherical coordinate system but with the azimuthal angle theta going from
        90 degrees at the North Pole, to -90 degrees at the South Pole, through
        0. at the equator. 
        
        This function is not equipped to handle wrap-around the ranges of theta
        phi and therefore does not work at the poles.
     
        Parameters
        ----------
        phi: float, mandatory, degrees
    	center of the spherical patch in ra with range 
        theta: float, mandatory, degrees
        delta: float, mandatory, degrees
        size: int, mandatory
            number of samples
        seed : int, optional, defaults to 1
            random Seed used for generating values
        Returns
        -------
        tuple of (phivals, thetavals) where phivals and thetavals are arrays of 
            size size in degrees.
        """
        u = rng.uniform(size=size)
        v = rng.uniform(size=size)
        phi = np.radians(phi)
        theta = np.radians(theta)
        delta = np.radians(delta)
    
        phivals = 2. * delta * u + (phi - delta)
        phivals = np.where ( phivals >= 0., phivals, phivals + 2. * np.pi)
        
        # use conventions in spherical coordinates
        theta = np.pi/2.0 - theta
     
        thetamax = theta + delta
        thetamin = theta - delta
    
        if thetamax > np.pi or thetamin < 0. :
            raise ValueError('Function not implemented to cover wrap around poles')
    
        # Cumulative Density Function is cos(thetamin) - cos(theta) / cos(thetamin) - cos(thetamax)
        a = np.cos(thetamin) - np.cos(thetamax)
        thetavals = np.arccos(-v * a + np.cos(thetamin))
    
        # Get back to -pi/2 to pi/2 range of decs
        thetavals = np.pi/2.0 - thetavals 
        return np.degrees(phivals) , np.degrees(thetavals)

class DisjointParams(with_metaclass(abc.ABCMeta, SNParamDistribution,
                                    SALT2Parameters,
                                    PositionSamples,
                                    RateDistributions)):
    pass

