"""
Implement a concrete `Tiling` class on the basis of healpix tiles. The
heavy lifting is done by the package OpSimSummary.
"""
from __future__ import absolute_import, print_function, division
import healpy as hp
import numpy as np
import opsimsummary  as oss
import pandas as pd
from .tessellations import Tiling
from sqlalchemy import create_engine

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
        self._preComputedMap = preComputedMap
        if self.hpOpSim is None and self.preComputedMap is None:
            raise ValueError('hpOpSim and preComputedMap cannot both be None')
        self._preComputedEngine = None

    @property
    def preComputedMap(self):
        if self._preComputedMap is not None:
            if not self._preComputedMap.startswith('sqlite'):
                self._preComputedMap = 'sqlite:///' + self._preComputedMap

        return self._preComputedMap

    @property
    def preComputedEngine(self):
        engine = self._preComputedEngine
        if engine is None:
            engine = create_engine(self.preComputedMap, echo=False)
        return engine

    @property
    def tileIDSequence(self):
        return xrange(self.npix)

    def area(self, tileID):
        if tileID not in self.tileIDSequence:
            raise ValueError('parameter tileID not in set of healpixIDs')
        return self._tileArea * np.degrees(1.) * np.degrees(1.)

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

        sql = 'SELECT obsHistID FROM {0} WHERE ipix == {1}'\
            .format(tableName, tileID)
        return pd.read_sql_query(sql, con=self.preComputedEngine)\
            .values.flatten()

    def _pointingFromHpOpSim(self, tileID):
        return self.hpOpSim.obsHistIdsForTile(tileID)


    def _tileFromHpOpSim(self, pointing):
        return self.hpOpSim.set_index('obsHistID').ix(pointing)['hids']

    def _tileFromPreComputedDB(self, pointing, tableName='simlib'):
        sql = 'SELECT ipix FROM {0} WHERE obsHistID == {1}'\
            .format(tableName, pointing)
        return pd.read_sql_query(sql, con=self.preComputedEngine)\
            .values.flatten()

    def tilesForPointing(self, pointing, alltiles=None, **kwargs):
        """
        return a maximal sequence of tile ID s for a particular OpSim pointing
        """
        if self.preComputedMap is not None:
            return self._tileFromPreComputedDB(self, pointing, tableName='simlib')
        elif self.hpOpSim is not None:
            return self._tileFromHpOpSim(self, pointing)
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
            return self._pointingFromHpOpSim(tileID)
        else:
            raise ValueError('both attributes preComputedMap and hpOpSim cannot'
                             ' be None')
                

    def _angularSamples(self, phi_c, theta_c, radius, numSamples, tileID, rng):
        phi, theta = super(self.__class__, self).samplePatchOnSphere(phi=phi_c,
								     theta=theta_c,
                                                                     delta=radius, 
                                                                     size=numSamples,
                                                                     degrees=False
                                                                     rng=rng)
        tileIds = hp.ang2pix(nside=self.nside, theta=np.radians(theta),
			     phi=np.radians(phi), nest=True)
        inTile = tileIds == tileID
        return phi[inTile], theta[inTile]
        

    def positions(self, tileID, numSamples, rng=None):
        """
        Return a tuple of (res_phi, res_theta) where res_phi and res_theta are
        spatially uniform samples  of positions of size numSamples.

        Parameters
        ---------
        tileID : int, mandatory

        numSamples : number of positions required

        rng : instance of `np.random.RandomState`


        Returns
        -------
        """
        # set the random number generator seed
        if rng is None:
            rng = np.random.RandomState(tileID)

        # zero return arrays
        res_theta = np.zeros(numSamples)
        res_phi = np.zeros(numSamples)

        # Set the center of the patch at the vertex of the tile
        theta_c, phi_c = hp.pix2ang(nside=self.nside, ipix=tileID, nest=True)
        radius = 2 * np.sqrt(self.area(tileID) / np.pi)

        # number of 
        num_already = 0

        while numSamples > 0:
            phi, theta = self._angularSamples(phi_c, theta_c, radius=radius,
                                              numSamples=numSamples,
                                              tileID=tileID,
                                              rng=rng)
            # print(self.nside, tileID, self.area(tileID))
            s = rng.get_state()
            print(s)
            num_obtained = len(phi)
            res_phi[num_already:num_obtained + num_already] = phi
            res_theta[num_already:num_obtained + num_already] = theta
            num_already += num_obtained
            numSamples -= num_obtained

        return res_phi, res_theta

