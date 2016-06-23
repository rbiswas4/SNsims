import numpy as np

__all__ = ['snRate', 'numSN']

def snRate(z, cosmo, alpha=2.6e-5, beta=1.5):
    """
    The rate of SN in units of number of SN/ comoving volume in Mpc^3/yr
    in earth years according to the commonly used power-law expression 
    .. math:: rate(z) = \alpha (h/0.7)^3 (1.0 + z)^\beta
    
    Parameters
    ----------
    z : array-like, mandatory 
        redshift
    cosmo : Instance of `astropy.cosmology` class, mandatory
        data structure specifying the cosmological parameters 
    alpha : float, optional, defaults to 2.6e-5
        
    beta : float, optional, defaults to 1.5
        parameter in expression  

    Example :
    ------

    """
    res = alpha *(1.0 + z)**beta 
    res *= ((cosmo.h/0.7)**3.)  
    return res

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

