#!/usr/bin/env python
from __future__ import absolute_import, print_function, division
from future.utils import with_metaclass
import abc

from .tessellations import Tiling
from .universe import Universe

__all__ = ['SimulationTile']
class SimulationTile(snsims.Universe):
    def __init__(self, paramDist, rate, NSIDE, tileID, hpOpSim, allPointings=None, timeRange=None):
        self._randomState = None
        self.Tiling = snsims.HealpixTiles(nside=NSIDE, preComputedMap=hpOpSim)
        self.tileID = tileID
        self.fieldArea = self.Tiling.area(tileID)
        self.zdist = rate(rng=self.randomState, fieldArea=self.fieldArea)
        self.zsamples = self.zdist.zSamples
        self.numSN = len(self.zsamples)
        self.positions = self.Tiling.positions(self.tileID, self.numSN, rng=self.randomState)
        self._snParamTable = None
        self.columns = ('expMJD', 'filter', 'fieldID', 'fiveSigmaDepth')
        self.tilePointings = self.Tiling.pointingSequenceForTile(self.tileID, 
                                                                 allPointings=allPointings, columns=self.columns)
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
        
        if self.minPeakTime is None or self.maxPeakTime is None:
            pass
        else:
            sp['t0'] = self.minPeakTime + (self.maxPeakTime - self.minPeakTime) * sp['t0']
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
        lcMinTime = t0 - 20. * (1.0 + z)
        lcMaxTime = t0 + 50. * (1.0 + z )
        df = self.tilePointings.query('expMJD < @lcMaxTime and expMJD > @lcMinTime')
        df['snid'] = snid
        fluxes = []
        fluxerrs = []
        for rows in df.iterrows():
            row = rows[1]
            # print(row['expMJD'], row['filter'], row['fiveSigmaDepth'])
            bp = self.bandPasses[row['filter']]
            flux = sn.catsimBandFlux(bandpassobject=bp, time=row['expMJD'])
            fluxerr = sn.catsimBandFluxError(time=row['expMJD'], bandpassobject=bp, fluxinMaggies=flux,
                                             m5=row['fiveSigmaDepth'])
            fluxes.append(flux)
            fluxerrs.append(fluxerr)
        df['flux'] = fluxes
        df['fluxerrs'] = fluxerrs
        return sn, df
    def writeTile(self, fileName, timeRange='model'):
        count = 0
        for snid in self.snParamTable.index.values:
            self.writeSN(snid, fileName, timeRange=timeRange)
            if count % 50 == 0:
                print('another 50', snid)
            count += 1
    def writeSN(self, snid, fileName, timeRange='model'):
        sn, df = self.SN(snid, timeRange)
        # with pd.get_store(fileName) as store:
        #     store.append('tile_{}'.format(self.tileID), df)
        df.to_hdf(fileName, key='tile_{}'.format(self.tileID), mode='a', format='t')
