#!/usr/bin/env python
from __future__ import absolute_import, print_function, division
from future.utils import with_metaclass
import abc

from .tessellations import Tiling
from .universe import Universe
from .paramDistribution import SimpleSALTDist
from .healpixTiles import HealpixTiles
from .lightcurve import LightCurve
import os
import numpy as np
import pandas as pd
from lsst.sims.photUtils import BandpassDict
from lsst.sims.catUtils.supernovae import SNObject

__all__ = ['SimulationTile', 'EntireSimulation', 'TiledSimulation']

class EntireSimulation(Universe):
    """
    Simulation of a set of SN from a set of telescope pointings
    and a set of SN. The simulation is perfectly reproducible if
    both the pointings, paramsDF are the same (in terms of ordering)
    
    Parameters
    -----------
    rng : instance of `numpy.random.RandomState` 
    pointings: instance of `pd.DataFrame` 
    	dataFrame with a minimal set of columns
	[`expMJD`, `filter`, `fiveSigmaDepth`]
    paramsDF : `pd.DataFrame`
	the minimal set of columns are 
	[`snid`, `x0`, `t0`, `x1` , `c` , `snra`, `sndec`]

    Attributes
    ----------
    randomState : `numpy.random.RandomState`
    snParams : `pd.DataFrame`

    """
    def __init__(self, rng, pointings, paramsDF, angularUnits='radians'):
        self.pointings = pointings
        self._paramsDf = paramsDF
        self._rng = rng
        self.angularUnits = angularUnits
        self.bandPasses = BandpassDict.loadTotalBandpassesFromFiles()
    
    @property
    def randomState(self):
        return self._rng

    @property
    def snParams(self):
        return self._paramsDf

    @staticmethod
    def getSNCosmoParamDict(odict, SNCosmoModel):
        mydict = dict()
        param_names = SNCosmoModel.param_names

        for param in odict.index.values:
            if param in param_names:
                mydict[param] = odict[param]
        return mydict
    
    def SN(self, snid, timeRange='model'):
        mySNParams = self.snParams.ix[snid]
        if self.angularUnits == 'radians':
            myra = np.radians(mySNParams.snra)
            mydec = np.decdians(mySNPadecms.sndec)
        elif self.angularUnits = 'degrees':
            myra = mySNParams.snra
            mydec = mySNPadecms.sndec
        sn = SNObject(ra=myra, dec=mydec)
        sncosmo_params = self.getSNCosmoParamDict(mySNParams, sn)
        sn.set(**sncosmo_params)
        return sn
    
    def lc(self, snid):
        sn = self.SN(snid, timeRange='model')
        lcMinTime = sn.mintime()
        lcMaxTime = sn.maxtime()
        # lcMinTime = self.SN(snid, timeRange='model').mintime()
        # lcMaxTime = self.SN(snid, timeRange='model').maxtime()
        if lcMinTime is None or lcMaxTime is None:
            df = self.pointings.copy()
        else:
            df = self.pointings.query('expMJD < @lcMaxTime and expMJD > @lcMinTime').copy()
        df['snid'] = snid
        fluxerr = np.zeros(len(df))
        modelFlux = np.zeros(len(df))
        for i, rowtuple in enumerate(df.iterrows()):
            row = rowtuple[1]
            # print(row['expMJD'], row['filter'], row['fiveSigmaDepth'])
            bp = self.bandPasses[row['filter']]
            modelFlux[i] = self.staticModelFlux(sn, row['expMJD'],
                                                bandpassobject=bp)
            fluxerr[i] = sn.catsimBandFluxError(time=row['expMJD'],
                                                bandpassobject=bp,
                                                fluxinMaggies=modelFlux[i],
                                                m5=row['fiveSigmaDepth'])

        rng = self.randomState
        df['fluxerr'] = fluxerr
        deviations = rng.normal(size=len(df))
        df['deviations'] = deviations
        df['zp'] = 0.
        df['ModelFlux'] = modelFlux
        df['flux'] = df['ModelFlux'] + df['deviations'] * df['fluxerr']
        df['zpsys']= 'ab'
        lc = df[['snid', 'expMJD', 'filter', 'ModelFlux', 'fieldID', 'flux', 'fluxerr',
                 'zp', 'zpsys', 'fieldID']]
        return LightCurve(lc)
    
    @staticmethod
    def staticModelFlux(sn, time, bandpassobject):
          return sn.catsimBandFlux(bandpassobject=bandpassobject,
                                   time=time)
    def modelFlux(self, snid, time, bandpassobject):
        # assert len(times) == len(bands)
        # flux = np.zeros(len(times))
	sn = self.SN(snid)
	return self.staticModelFlux(sn, time=time, bandpassobject=bandpassobject)
        # return self.SN(snid).catsimBandFlux(bandpassobject=bandpassobject,
        #                                    time=time)


    def writeSNParams(self, paramFileName, IDVal=0):
	"""
	Write the dataframe `self.snParams` to a file 

	Parameters
	----------
	paramFileName : Instance of string
	    paramFileName
	IDVal : integer
	    used as a key to write a group
	"""

        if paramFileName.endswith('.hdf'):
            self.snParams.to_hdf(paramFileName, key='{}'.format(IDVal))
        else:
            raise NotImplementedError('Only methods to write to hdf files'
                                      'implemented')

    def writeSN(self, snid, fileName, IDVal=0, timeRange='model'):
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
        lc = self.lc(snid)
        df = lc.lightCurve 
        df['band'] = df['band'].astype(str)
        with pd.get_store(fileName) as store:
            store.append('tile_{}'.format(IDVal), df)

class TiledSimulation(EntireSimulation):

    def __init__(self,
		 paramDF,
                 NSIDE,
                 tileID,
                 hpOpSim,
		 rng=None,
                 allPointings=None,
                 timeRange=None):
	"""
	Parameters
	----------
	paramDF
	"""
	self.tileID = tileID
        self._randomState = rng
	if self._randomState is None:
	    self._randomState = np.random.RandomState(self.tileID)
        self.Tiling = HealpixTiles(nside=NSIDE, preComputedMap=hpOpSim)
        self.fieldArea = self.Tiling.area(self.tileID)
        self.columns = ('expMJD', 'filter', 'fieldID', 'fiveSigmaDepth')
        self.tilePointings = self.Tiling.pointingSequenceForTile(self.tileID, 
                                                                 allPointings=allPointings,
                                                                 columns=self.columns)
	super(TiledSimulation, self).__init__(rng=self._randomState,
					      pointings=self.tilePointings,
					      paramsDF=paramDF)

class SimulationTile(Universe):
    def __init__(self,
                 paramDist,
                 rate,
                 NSIDE,
                 tileID,
                 hpOpSim,
                 allPointings=None,
                 timeRange=None,
                 angularUnits='radians'):

        self._randomState = np.random.RandomState(self.tileID)
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
        sncosmo_params = self.getSNCosmoParamDict(mySNParams, sn)
        sn.set(**sncosmo_params)
        z = sn.get('z')
        t0 = sn.get('t0')
        # lcMinTime = t0 - 20. * (1.0 + z)
        # lcMaxTime = t0 + 50. * (1.0 + z )
        return sn

    @staticmethod
    def staticModelFlux(sn, time, bandpassobject):
        
        return sn.catsimBandFlux(bandpassobject=bandpassobject,
                                 time=time) 

    def modelFlux(self, snid, time, bandpassobject):

        # assert len(times) == len(bands)
        # flux = np.zeros(len(times))


        # flux = np.asarray(list(self.SN(snid).catsimBandFlux(bandpassobject=self.bandPasses[bands[i]],
        # time=times[i]) for i in range(len(times)))) 
        #for i, band in enumerate(bands):
        #    bp = self.bandPasses[band]
        #    flux[i] = self.SN(snid).catsimBandFlux(bandpassobject=bp, time=times[i])
        # print(len(flux), len(times))
        return self.SN(snid).catsimBandFlux(bandpassobject=bandpassobject,
                                                        time=time) 



    def lc(self, snid):
        sn = self.SN(snid, timeRange='model')
        lcMinTime = sn.mintime()
        lcMaxTime = sn.maxtime()
        # lcMinTime = self.SN(snid, timeRange='model').mintime()
        # lcMaxTime = self.SN(snid, timeRange='model').maxtime()
        if lcMinTime is None or lcMaxTime is None:
            df = self.tilePointings.copy()
        else:
            df = self.tilePointings.query('expMJD < @lcMaxTime and expMJD > @lcMinTime').copy()
        df['snid'] = snid
        fluxerr = np.zeros(len(df))
        modelFlux = np.zeros(len(df))
        for i, rowtuple in enumerate(df.iterrows()):
            row = rowtuple[1]
            # print(row['expMJD'], row['filter'], row['fiveSigmaDepth'])
            bp = self.bandPasses[row['filter']]
            modelFlux[i] = self.staticModelFlux(sn, row['expMJD'],
                                                bandpassobject=bp)
            fluxerr[i] = sn.catsimBandFluxError(time=row['expMJD'],
                                                bandpassobject=bp,
                                                # fluxinMaggies=row['ModelFlux'],
                                                fluxinMaggies=modelFlux[i],
                                                m5=row['fiveSigmaDepth'])

        rng = self.randomState
        df['fluxerr'] = fluxerr
        deviations = rng.normal(size=len(df)) 
        df['deviations'] = deviations
        df['zp'] = 0.
        df['ModelFlux'] = modelFlux
        df['flux'] = df['ModelFlux'] + df['deviations'] * df['fluxerr']
        df['zpsys']= 'ab'
        lc = df[['snid', 'expMJD', 'filter', 'ModelFlux', 'fieldID', 'flux', 'fluxerr',
                 'zp', 'zpsys', 'fieldID']]
        return LightCurve(lc)

    def writeTile(self, fileName, timeRange='model', paramFileName=None):
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
        if paramFileName is None:
            filename_parts = fileName.split('.')
            filename_parts[-2] += '_params'
            paramFileName = '.'.join(filename_parts)
        self.writeSNParams(paramFileName)

    def writeSNParams(self, paramFileName):
        if paramFileName.endswith('.hdf'):
            self.snParamTable.to_hdf(paramFileName, key='{}'.format(self.tileID))
        else:
            raise NotImplementedError('Only methods to write to hdf files'
                                      'implemented')

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
        lc = self.lc(snid)
        df = lc.lightCurve 
        df['band'] = df['band'].astype(str)
        with pd.get_store(fileName) as store:
            store.append('tile_{}'.format(self.tileID), df)


