import os
import numpy as np
import pandas as pd
import time
import opsimsummary as oss

import snsims
import healpy as hp


zdist = snsims.PowerLawRates(rng=np.random.RandomState(1), 
                             fieldArea=hp.nside2pixarea(nside=4, degrees=True), 
                             zbinEdges=np.arange(0.0001, 1.1, 0.1))


print (snsims.__VERSION__)
print(oss.__VERSION__)
opsimOut = oss.OpSimOutput.fromOpSimHDF('../../../data/LSST/OpSimData/minion_1016.hdf',
                                        subset='combined')
NSIDE = 256

def func(i):
    tileID = i
    simTile = snsims.SimulationTile(snsims.SimpleSALTDist,
                                    NSIDE=NSIDE,
                                    tileID=tileID,
                                    hpOpSim='../../OpSimSummary/scripts/healpixelized_MINION_1016_256.db',
                                    rate=snsims.PowerLawRates,
                                    allPointings=opsimOut)
    simTile.writeTile(fileName='simTiles_{}.hdf'.format(i))
tstart = time.time()
func(0)
tend = time.time()
print(tstart, tend, tend-tstart)
