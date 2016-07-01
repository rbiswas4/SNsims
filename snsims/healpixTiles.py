"""
Implement a concrete `Tiling` class on the basis of healpix tiles.
"""
from __future__ import absolute_import, print_function, division
import healpy as hp

from .tessellations import Tiling

__all__ = ['HealpixTiles']

class HealpixTiles(Tiling):
    """
    A concrete Tiling class based on Healpix Tiles. The user is
    allowed to choose the following parameters:
    NSIDE:


    Attributes
    ----------

    nside : int, power of 2, defaults to 256
        healpix nside parameter

    """
    def __init__(self,
                 nside=256,
                 dbname=None,
                 dbConn=None):
        """
        nside : int, power of 2, defaults to 256
            nside parameter of healpix. determines the size of the tiles
            so that there are 12 * nside **2 equally sized tiles covering
            the sphere.
        """
        self.nside = nside
        self.npix = hp.nside2npix(nside)
        self._tileArea = hp.nside2pixarea(nside)

    @property
    def tileIDSequence(self):
        return xrange(self.npix)

    def area(self, tileID):
        if tileID not in self.tileIDSequence:
            raise ValueError('parameter tileID not in set of healpixIDs')
        return self._tileArea

    def tileIDsForSN(self, ra, dec):
        """
        Parameters
        ----------
        ra : `numpyp.ndarray` or float, degrees, mandatory
        dec : `numpy.ndarray` or float, degrees, mandatory
        """
        # If scalar float or list, convert to numpy array
        dec = np.ravel(dec)
        ra = np.ravel(ra)

        # Convert to usual spherical coordinates
        theta = - np.radians(dec) + np.pi/2.
        phi = np.radians(ra)

        inds = hp.ang2pix(nside=self.nside, theta=theta, phi=phi, nest=True)
        return inds

    def pointingSequenceForTile(self, tileID, allPointings, **kwargs):
        return None
                
