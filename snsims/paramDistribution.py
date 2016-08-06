"""
Module with concrete implementations of functions important for sampling
parameters
"""
from __future__ import absolute_import, print_function
from .universe import RateDistributions
from astropy.cosmology import Planck15
import numpy as np


__all__ = ['PowerLawRates']

class PowerLawRates(RateDistributions):

    def __init__(self, alpha=2.6e-5, beta=1.5, cosmo=Planck15):
        """
        Parameters
        ----------
        cosmo : Instance of `astropy.cosmology` class, optional, defaults to Planck15
            data structure specifying the cosmological parameters 
        alpha : float, optional, defaults to 2.6e-5
        beta : float, optional, defaults to 1.5
            parameter in expression  
        """
        self.alpha = alpha
        self.beta = beta
        self.cosmo = cosmo

    def snRate(self, z):
        """
        The rate of SN at a redshift z in units of number of SN/ comoving
        volume in Mpc^3/yr in earth years according to the commonly used
        power-law expression 

        .. math:: rate(z) = \alpha (h/0.7)^3 (1.0 + z)^\beta
        
        Parameters
        ----------
        z : array-like, mandatory 
            redshifts at which the rate is evaluated

        Examples
        --------
        """
        res = self.alpha * (1.0 + z)**self.beta 
        res *= ((self.cosmo.h / 0.7) **3.)  
        return res

    def zSampleSize(self, zbinEdges, DeltaT, skyFraction=None, fieldArea=None,
                    zlower=None, zhigher=None, numBins=None):
        """
        Parameters
        ----------
        zbinEdges : `nunpy.ndarray` of edges of zbins, defaults to None
            Should be of the form np.array([z0, z1, z2]) which will have
            zbins (z0, z1) and (z1, z2)
        skyFraction : np.float, optional, 
        fieldArea : optional, units of degrees squared
            area of sky considered.
        zlower : float, optional, defaults to None
            lower edge of z range
        zhigher : float, optional, defaults to None
            higher edge of z range
        numBins : int, optional, defaults to None
           if not None, overrides zbinEdges
        """
        if fieldArea is None and skyFraction is None:
            raise ValueError('both fieldArea and skyFraction cannot be None')
        elif fieldArea is not None and skyFraction is not None:
            raise ValueError('both fieldArea and skyFraction cannot be given')
        if fieldArea is not None: 
            skyFraction = fieldArea * np.radians(1.)**2.0  / 4.0 / np.pi

        if zbinEdges is None:
            if any(_x is None for _x in (zlower, zhigher, numBins)):
                raise valueError('zbinEdges or zlower, zhigher, numbins must'
                                 'be given')
            else:
                zbinEdges = np.linspace(zlower, zhigher, numBins + 1)

        # midpoints of z bins where the rate is evaluated
        z_mids = 0.5 * (zbinEdges[1:] + zbinEdges[:-1])
        snpervolume = self.snRate(z_mids) 

        # Comoving volume of the univere in between zlower and zhigher
        vols = self.cosmo.comoving_volume(zbinEdges)

        # Comoving volume in each bin
        vol = vols[1:] - vols[:-1]
        vol *= skyFraction
        
        numSN = vol * snpervolume * DeltaT
        return numSN.value

    
    def numSN(zlower, zhigher, numBins, cosmo, fieldArea, DeltaT, snrate):
        """
        Return the number of expected supernovae in time DeltaT, with a rate snrate
        in a redshift range zlower, zhigher divided into numBins equal redshift
        bins. The variation of the rate within a bin is ignored.
    
        Parameters
        ----------
        zlower : mandatory, float
            lower limit on redshift range
        zhigher : mandatory, float
            upper limit on redshift range
        numBins : mandatory, integer
            number of bins
        cosmo : `astropy.cosmology` instance, mandatory
            cosmological parameters
        fieldArea : mandatory, units of radian square
            sky area considered
        """
        z_edges = np.linspace(zlower, zhigher, numBins + 1)
        z_bins = 0.5 * (z_edges[1:] + z_edges[:-1])
        
        # Comoving volume of the univere in between zlower and zhigher
        vol = cosmo.comoving_volume(z_edges[1:]) - \
              cosmo.comoving_volume(z_edges[:-1])
        
        # fullvol = vol.sum()
        vol *= fieldArea / 4.0 / np.pi
        sn = snrate(z_bins, cosmo)
        
        # normalize the time window to rest frame time
        numSN = DeltaT * sn * vol / (1. + z_bins)
        return numSN.value# , fullvol

