#!/usr/bin/env python
import os
import opsimsummary as oss
import snsims
import pytest
import unittest

datadir = os.path.join(oss.__path__[0], 'example_data')
opsimdb = os.path.join(datadir, 'enigma_1189_micro.db')
precomp = os.path.join(datadir, 'healpixels_micro.db')

class test_healpixTiles(unittest.TestCase):
    def setUpClass(self):
        hpOpSim = oss.HealPixelizedOpSim.fromOpSimDB(opsimdb)

