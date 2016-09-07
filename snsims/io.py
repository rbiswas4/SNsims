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
import glob

class OpSimField(object):

    def __init__(fieldID, dirname, doneTiles=None, maxLog2ObsHistID=22):
        self.dirname = self.abspath(dirname)
        self.fieldID = fieldID
        self.doneTiles = set(doneTiles)
        self.thisRun = None
        self.validate()

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
        fieldID = lcs.iloc[:, fieldID_inds].values
        good_cols = np.where(colnames!='fieldID')
        lcs = lcs[good_cols].copy()
        lcs['fieldID'] = fieldID 

        return lcs
    def lcDF(self, toDoList=None):
        if toDoList is None:
            toDoList = self.thisRun
        lcFiles = list(os.path.join(self.dirname,
                                     'simTiles_{}.hdf'.format(tileID)))
        lclist = list(pd.read_hdf(lcFile).reset_index() for lcFile in lcFiles)
        lcs = pd.concat(lclist, ignore_index=True)
        lcs['indVals'] = np.left_shift(lcs.snid, maxLog2ObsHistID) + \
                         lcs.obsHistID
        lcs = self.cleanup(lcs)
        return lcs.set_index('indVals', inplace=True)

    def writeFiles(self, toDoList=None, baseFileName=None):
        if toDoList is None:
            toDoList = self.thisRun
        if baseFileName is None:
            baseFileName = 'SimFielID_{}'.format(self.fieldID)
        lcdf = self.lcDF(toDoList)
        pdf = self.paramDF(toDoList)
        lcdf.write(baseFileName +'_lcs.hdf')
        pdf.write(baseFileName +'_params.hdf')

    @property
    def runTiles(self):

        hdfFiles = glob.glob(self.dirname + '/*params.hdf')
        tiles = list(int(fname.split('.hdf')[0].split('params')[-1])
                     for fname in hdfFiles)
        return set(tiles)

    @property
    def fieldTiles(self):
        csvFile = os.path.join(self.dirname,
                               'tiles_{}.csv '.format(self.fieldID))
        assert os.path.file.exists(csvFile)
        tileIDs = set(pd.read_csv(csvFile).values.tolist())
        return tileIDs

    def validate(self):
        csvFile = os.path.join(self.dirname,
                               'tiles_{}.csv '.format(self.fieldID))
        toDo = tileIDs - self.doneTiles
        assert toDo == self.runTiles
        self.thisRun = toDo




