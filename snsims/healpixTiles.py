"""
Implement a concrete `Tiling` class on the basis of healpix tiles. The
heavy lifting is done by the package OpSimSummary.
"""
from __future__ import absolute_import, print_function, division
import healpy as hp
import numpy as np
import opsimsummary  as oss

from .tessellations import Tiling

__all__ = ['HealpixTiles']


class HealpixTiles(Tiling):
    """
    A concrete Tiling class based on Healpix Tiles. The user is
    allowed to choose the following parameters:

    Attributes
    ----------
    nside : int, power of 2, defaults to 256
        healpix nside parameter

    """
    def __init__(self,
                 nside=256,
                 healpixelizedOpSim=None,
                 preComputedMap=None):
        """
        nside : int, power of 2, defaults to 256
            nside parameter of healpix. determines the size of the tiles
            so that there are 12 * nside **2 equally sized tiles covering
            the sphere.
        """
        self.nside = nside
        self.npix = hp.nside2npix(nside)
        self._tileArea = hp.nside2pixarea(nside)
        self.hpOpSim = healpixelizedOpSim
        self.preComputedMap = preComputedMap
        if self.hpOpSim is None and self.preComputedMap is None:
            raise ValueError('hpOpSim and preComputedMap cannot both be None')

    @property
    def preComputedEngine(self):
        engine = self.preComputedEngine
        if engine is not None:
            engine = create_engine(self.preComputedMap, echo=False)
        return engine



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

    def _pointingFromPrecomputedDB(self, tileID, tableName='simlib'):

        tName = tableName
        sql = 'SELECT obsHistID FROM @tName WHERE ipix == {}'.format(tileID)
        return pd.read_sql_query(sql, con=self.preComputedEngine).values

    def _pointingFromHpOpSim(self, tileID):
        return self.hpOpSim.obsHistIdsForTile(tileID)


    def _tileFromHpOpSim(self, pointing):
        return self.hpOpSim.set_index('obsHistID').ix(pointing)['hids']

    def _tileFromPreComputedDB(self, pointing, tableName='simlib'):
        tName = tableName
        sql = 'SELECT ipix FROM @tName WHERE obsHistID == {}'.format(pointing)
        return pd.read_sql_query(sql, con=self.preComputedEngine).values

    def tilesForPointing(self, pointing, alltiles=None, **kwargs):
        """
        return a maximal sequence of tile ID s for a particular OpSim pointing
        """
        if self.preComputedMap is not None:
            return _tileFromPreComputedDB(self, pointing, tableName='simlib')
        elif self.hpOpSim is not None:
            return _tileFromHpOpSim(self, pointing)
        else:
            raise ValueError('both attributes preComputedMap and hpOpSim cannot'
                             ' be None')

    def pointingSequenceForTile(self, tileID, allPointings, **kwargs):
        """
        return a maximal sequence of pointings for a particular tileID
        """
        if self.preComputedMap is not None:
            return self._pointingFromPrecomputedDB(tileID, tableName='simlib')
        elif self.hpOpSim is not None:
            return _pointingFromHpOpSim(self, tileID)
        else:
            raise ValueError('both attributes preComputedMap and hpOpSim cannot'
                             ' be None')
                
