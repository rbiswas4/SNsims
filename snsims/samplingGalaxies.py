from __future__ import absolute_import, division, print_function
import numpy as np
from numpy.random import gamma
from scipy.special import gammainc

__all__ = ['SersicSamples']
class SersicSamples(object):
    """
    Class for sampling sersic profiles in CatSim
    """
    def __init__(self, rng):
        """
        Parameters
        ---------
        rng : instance of `numpy.random.RandomState`

        """
        self.rng = rng
        self.fp4 = np.arange(0., 200., 0.01)
        self.xp4 = gammainc(8, self.fp4)
        self.fp1 = np.arange(0., 20., 0.001)
        self.xp1 = gammainc(2, self.fp1)

    @staticmethod
    def sampleAngles(a, b, numSamples=1, rng=np.random.RandomState()):
        """
        return a sample of the angle with respect to a position angle
        in units of degrees. For a single float, and `numSamples=1` the
        answer will still be in an array of length 1.

        Parameters
        ----------
        a : np.float, or `np.ndarray`
            semi-major axis of ellipe. If array, must have len of `numSamples`
        b : np.float, or `np.ndarray`
            semi-minor axis of ellipe. If array, must have len of `numSamples`
        numSamples: int, defaults to 1
            number of samples desired. must match the length of a and b if a
            and b are arrays.
        rng : `np.random.RandomState` instance, defaults to no argument

        Returns
        -------
        `np.ndarray` of length `numSamples` in units of degrees.

        .. note:: The two lengh parameters a and b need to have the same unit.
            The binary parameter is used to distribute the result in all four
            quadrants from the two that are in the range of arctan.
        """
        if isinstance(a, np.float):
            assert isinstance(b, np.float)
            if numSamples >=1:
                a = np.ones(numSamples)*a
                b = np.ones(numSamples)*b
        if len(a) != len(b) or len(b) != numSamples:
            raise ValueError('a, b, numSamples must have same lengths')
        u = rng.uniform(0., 2.0*np.pi, size=numSamples)
        # Use the binary parameter to handle other quadrants from arctan
        binary = np.random.choice([0, 1], size=numSamples, p=[0.5, 0.5])
        return np.degrees(np.arctan(b*np.tan(u)/a) + binary*np.pi)

    def sampleRadius(self, halfLightRadius, numSamples=1, sersicIndex=1):
        """
        return samples of the position, given the halfLightRadius and the
        sersicIndex. The answers are for sersic index of 1 and 4 only.  

        Parameters
        ----------
        halfLightRadius : `np.float` or `np.ndarray` of dtype `np.float`
            half light radius of the galaxy/bulge/disk
        numSamples : int, defaults to 1
            number of samples desired
        sersicIndix : np.float, defaults to 1
            sersic index, works for bulges and disks only
        Returns
        -------
        `np.ndarray` of radial distances from the center of the galaxy in the
        same units as halfLightRadius
        """
        if isinstance(halfLightRadius, np.float) and numSamples >= 1:
            halfLightRadius = np.ones(numSamples) * halfLightRadius
        elif numSamples != len(halfLightRadius):
            raise ValueError('The lengths must match')
        u = self.rng.uniform(size=numSamples)
        if sersicIndex == 1:
            x = np.interp(u, self.xp1, self.fp1)
            b1 = 1.678
            return  halfLightRadius * x / b1
        if sersicIndex == 4:
            x = np.interp(u, self.xp4, self.fp4)
            b4 = 7.669
            return halfLightRadius  * (x / b4) **4

