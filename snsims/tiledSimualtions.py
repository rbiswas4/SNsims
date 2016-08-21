#!/usr/bin/env python
from __future__ import absolute_import, print_function, division
from future.utils import with_metaclass
import abc

from .tessellations import Tiling
from .universe import Universe

__all__ = ['SimulationTile']

class SimulationTile(snsims.Universe):
    def __init__(self, paramDist, rate, NSIDE, tileID, hpOpSim):
        self._randomState = None
        self.Tiling = snsims.HealpixTiles(nside=NSIDE, healpixelizedOpSim=hpOpSim)
        self.tileID = tileID
        self.fieldArea = self.Tiling.area(tileID)
        self.zdist = rate(rng=self.randomState, fieldArea=self.fieldArea)
        self.zsamples = self.zdist.zSamples
        self.numSN = len(self.zsamples)
        self.positions = self.Tiling.positions(self.tileID, self.numSN, rng=self.randomState)
        self._snParamTable = None
    @property
    def snParamTable(self):
        if self._snParamTable is None:
            self.snParams()
        return self._snParamTable
    @property
    def randomState(self):
        if self._randomState is None:
            self._randomState = np.random.RandomState(self.tileID)
        return self._randomState
    def SN(self):
        pass
    def snParams(self):
        zsamples = self.zdist.zSamples
        numSN = len(zsamples)
        positions = self.Tiling.positions(self.tileID, numSN, rng=self.randomState)
        ra = self.positions[0]
        dec = - self.positions[1] + 45.0 
        # Why do we need numSN
        sp = snsims.SimpleSALTDist(numSN=numSN, rng=self.randomState, zSamples=self.zsamples).paramSamples
        sp['ra'] = self.positions[0]
        sp['dec'] = self.positions[1]
        sp['snid'] = np.left_shift(self.tileID, 20) + np.arange(numSN)
        sp.set_index('snid', inplace=True)
        self._snParamTable = sp
        return sp

