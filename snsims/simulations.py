#!/usr/bin/env python
from __future__ import absolute_import, print_function, division
from future.utils import with_metaclass
import abc

from .tessellations import Tiling
from .universe import Universe
from .paramDistribution import SimpleSALTDist
from .healpixTiles import HealpixTiles
import os
import numpy as np
import pandas as pd
from lsst.sims.photUtils import BandpassDict
from lsst.sims.catUtils.supernovae import SNObject

__all__ = ['SimulationTile']
class SimulationTile(Universe):
    def __init__(self,
                 paramDist,
                 rate,
                 NSIDE,
                 tileID,
                 hpOpSim,
                 allPointings=None,
                 timeRange=None):

        self._randomState = None
        self.Tiling = HealpixTiles(nside=NSIDE, preComputedMap=hpOpSim)
        self.tileID = tileID
        self.fieldArea = self.Tiling.area(tileID)
        self.zdist = rate(rng=self.randomState, fieldArea=self.fieldArea)
        self.zsamples = self.zdist.zSamples
        self.numSN = len(self.zsamples)
        self.positions = self.Tiling.positions(self.tileID, self.numSN,
                                               rng=self.randomState)
        self._snParamTable = None
        self.columns = ('expMJD', 'filter', 'fieldID', 'fiveSigmaDepth')
        self.tilePointings = self.Tiling.pointingSequenceForTile(self.tileID, 
                                                                 allPointings=allPointings,
                                                                 columns=self.columns)
        self._timeRange = timeRange
        self.bandPasses = BandpassDict.loadTotalBandpassesFromFiles()
        
    @property
    def minPeakTime(self):
        if  self._timeRange is None:
            minTime = self.tilePointings.expMJD.min()
        else:
            minTime = self._timeRange[0]
        return minTime
    
    @property
    def maxPeakTime(self):
        if  self._timeRange is None:
            maxTime = self.tilePointings.expMJD.max()
        else:
            maxTime = self._timeRange[1]
        return maxTime
    
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

    def snParams(self):
        zsamples = self.zdist.zSamples
        numSN = len(zsamples)
        positions = self.Tiling.positions(self.tileID, numSN,
                                          rng=self.randomState)
        ra = self.positions[0]
        dec = - self.positions[1] + 45.0 
        # Why do we need numSN
        sp = SimpleSALTDist(numSN=numSN, rng=self.randomState,
                            zSamples=self.zsamples).paramSamples
        sp['ra'] = self.positions[0]
        sp['dec'] = self.positions[1]
        sp['snid'] = np.left_shift(self.tileID, 20) + np.arange(numSN)
        sp.set_index('snid', inplace=True)
        self._snParamTable = sp
        if self.minPeakTime is None or self.maxPeakTime is None:
            pass
        else:
            sp['t0'] = self.minPeakTime + \
                       (self.maxPeakTime - self.minPeakTime) * sp['t0']
        return sp
    
    @staticmethod
    def getSNCosmoParamDict(odict, SNCosmoModel):
        mydict = dict()
        param_names = SNCosmoModel.param_names

        for param in odict.index.values:
            if param in param_names:
                mydict[param] = odict[param]
        return mydict
                
    def SN(self, snid, timeRange='model'):
        mySNParams = self.snParamTable.ix[snid]
        sn = SNObject(ra=mySNParams.ra, dec=mySNParams.dec)
        #print mySNParams
        sncosmo_params = self.getSNCosmoParamDict(mySNParams, sn)
        #print(sncosmo_params)
        sn.set(**sncosmo_params)
        z = sn.get('z')
        t0 = sn.get('t0')
        # lcMinTime = t0 - 20. * (1.0 + z)
        # lcMaxTime = t0 + 50. * (1.0 + z )
        return sn

    @staticmethod
    def modelFlux(snid, times, bands):

        assert len(times) == len(bands)
        flux = np.zeros(len(times))

        for i, band in enumerate(bands):
            bp = self.bandPasses[band]
            flux[i] = sn.catsimBandFlux(bandpassobject=bp, time=times[i])
        return flux



    def lc(self, snid):
        lcMinTime = self.SN(snid, timeRange='model').mintime()
        lcMaxTime = self.SN(snid, timeRange='model').maxtime()
        if lcMinTime is None or lcMaxTime is None:
            df = self.tilePointings
        else:
            df = self.tilePointings.query('expMJD < @lcMaxTime and expMJD > @lcMinTime')
        df['snid'] = snid
        df['ModelFlux'] = self.modelFlux(snid=snid, times=df.expMJD.values, bands=df.filter.values)
        fluxerr = np.zeros(len(df))
        for row in df.iterrows():
            # print(row['expMJD'], row['filter'], row['fiveSigmaDepth'])
            bp = self.bandPasses[row['filter']]
            fluxerr[i] = sn.catsimBandFluxError(time=row['expMJD'],
                                             bandpassobject=bp,
                                             fluxinMaggies=row['ModelFlux'],
                                             m5=row['fiveSigmaDepth'])

        rng = self.randomState
        df['fluxerr'] = fluxerr
        deviations = rng.normal(size=len(df)) 
        df['deviations'] = deviations
        df['flux'] = df['ModelFlux'] + df['deviations'] * df['fluxerr']
        df['zp'] = 0.
        df['zpsys']= 'ab'
        lc = df[['expMJD', 'filter', 'ModelFlux', 'fieldID', 'flux', 'fluxerr',
                 'zp', 'zpsys', 'fieldID']]
        return LightCurve(lc)

    def writeTile(self, fileName, timeRange='model'):
        """
        """
        count = 0
        for snid in self.snParamTable.index.values:
            self.writeSN(snid, fileName, timeRange=timeRange)
            if count % 50 == 0:
                if count == 0:
                    pass
                print('another 50', snid)
            count += 1

    def writeSN(self, snid, fileName, timeRange='model'):
        """
        Write light curve of SN to disc

        Parameters
        ----------
        snid : int/string
            SN id of SN 
        fileName : string, mandatory

        timeRange : string, optional, defaults to model
            time range over which the light curve is written to disk

        """
        sn, df = self.SN(snid, timeRange)
        df['filter'] = df['filter'].astype(str)
        with pd.get_store(fileName) as store:
            store.append('tile_{}'.format(self.tileID), df)
