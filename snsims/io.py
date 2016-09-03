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

class CollectFiles(object):

    def __init__(tileID, dirname):
        self.dirname = self.abspath(dirname)
        self.tileID = tileID
        self.validate()

    def validateInput(self):
        csvFile = os.path.join(self.dirname, 'tiles_{}.csv '.format(self.tileID))
        assert os.path.file.exists(csvFile)
        tiles = pd.read_csv(csvFile).values
        return tiles


