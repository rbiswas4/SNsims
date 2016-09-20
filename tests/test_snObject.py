from lsst.sims.catUtils.supernovae import SNObject
def test_loadSN():
    sn = SNObject(ra=0., dec=-27.5)
    sn.set(z=0.5)
    d = sn.SNstate
    assert len(d) > 0
    assert d['z'] == 0.5
    assert d['x1'] == 0.0

