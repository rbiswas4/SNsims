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
    def sampleAngles(numSamples, a, b, rng):
        u = rng.uniform(0., 2.0*np.pi, size=numSamples)
        binary = np.random.choice([0, 1], size=numSamples, p=[0.5, 0.5])
        return np.degrees(np.arctan(b*np.tan(u)/a) + binary*np.pi)

    def sampleRadius(self, numSamples, halfLightRadius, sersicIndex=1):
        if isinstance(halfLightRadius, np.float) and numSamples > 1:
            halfLightRadius = np.ones(numSamples) * halfLightRadius
        if numSamples != len(halfLightRadius):
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
