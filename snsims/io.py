#!/usr/bin/env python
"""
As this software develops, we will probably add several classes of writers and
readers. The objective of the io module is to try to abstract away the
differences in the implementation of these different writers.

Thus, all writers should have certain functionality:
    - read the output product of the simulations to provide SN light curves as
        instances of `snsims.LightCurve` class. This should be possible by
        providing it an output object and an SNID and should work transparently
        whether it is a SQL database, or an hdf file or a csv file. 
"""
from __future__ import absolute_import, division, print_function
import os
import pandas as pd
import glob
import numpy as np

__all__ = ['OpSimField']
class OpSimField(object):

    def __init__(self, fieldID, dirname, doneTiles=[], maxLog2ObsHistID=22):
        self.dirname = os.path.abspath(dirname)
        self.fieldID = fieldID
        self.doneTiles = set(doneTiles)
        self.thisRun = None
        self.validate()
        self.maxLog2ObsHistID = maxLog2ObsHistID

    def paramDF(self, toDoList=None):
        if toDoList is None:
            toDoList = self.thisRun
        hdfFiles = list(os.path.join(self.dirname,
                                     'simTiles_{}_params.hdf'.format(tileID))
                        for tileID in toDoList)
        dfs = list(pd.read_hdf(hdfFile) for hdfFile in hdfFiles)
        return pd.concat(dfs)
    @staticmethod
    def cleanup(df):
        colnames = df.columns.values
        fieldID_inds = np.where(colnames=='fieldID')[0] 
        fieldID = df.iloc[:, fieldID_inds].values
        good_cols = np.where(colnames!='fieldID')
        lcs = df[colnames[good_cols]].copy()
        lcs['fieldID'] = fieldID[:, 0] 
        return lcs

        return lcs
    def lcDF(self, toDoList=None):
        if toDoList is None:
            toDoList = self.thisRun
        lcFiles = list(os.path.join(self.dirname,
                                     'simTiles_{}.hdf'.format(tileID)) for tileID in self.thisRun)
        lclist = list(pd.read_hdf(lcFile).reset_index() for lcFile in lcFiles)
        lcs = pd.concat(lclist, ignore_index=True)
        lcs['indVals'] = np.left_shift(lcs.snid, self.maxLog2ObsHistID) + \
                         lcs.obsHistID
        lcs = self.cleanup(lcs)
        lcs.set_index('indVals', inplace=True)
        return lcs

    def writeFiles(self, toDoList=None, baseFileName=None):
        if toDoList is None:
            toDoList = self.thisRun
        if baseFileName is None:
            baseFileName = 'SimFielID_{}'.format(self.fieldID)
        lcdf = self.lcDF(toDoList)
        pdf = self.paramDF(toDoList)
        lcdf.to_hdf(baseFileName +'_lcs.hdf', key=str(self.fieldID))
        pdf.to_hdf(baseFileName +'_params.hdf', key=str(self.fieldID))

    @property
    def runTiles(self):

        hdfFiles = glob.glob(self.dirname + '/*params.hdf')
        tiles = list(fname.split('_')[2] for fname in hdfFiles)
        tiles = list(int(tile) for tile in tiles)
        return set(tiles)

    @property
    def fieldTiles(self):
        csvFile = os.path.join(self.dirname,
                               'tiles_{}.csv'.format(self.fieldID))
        assert os.path.exists(csvFile)
        tileIDs = set(pd.read_csv(csvFile, names=['tileId']).tileId.values.tolist())
        return tileIDs

    def validate(self):
        csvFile = os.path.join(self.dirname,
                               'tiles_{}.csv '.format(self.fieldID))
        toDo = self.fieldTiles - self.doneTiles
        assert toDo == self.runTiles
        self.thisRun = toDo




