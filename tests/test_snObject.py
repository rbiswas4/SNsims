from lsst.sims.catUtils.supernovae import SNObject
def test_loadSN():
    sn = SNObject(ra=0., dec=-27.5)
    sn.set(z=0.5)
    d = sn.SNstate
    assert len(d) > 0
    assert d['z'] == 0.5
    assert d['x1'] == 0.0


def test_attributeDefaults(self):
    """
    Check the defaults and the setter properties for rectifySED and
    modelOutSideRange
    """
    snobj = SNObject(ra=30., dec=-60., source='salt2')
    self.assertEqual(snobj.rectifySED, True)
    self.assertEqual(snobj.modelOutSideTemporalRange, 'zero')

    snobj.rectifySED = False
    self.assertFalse(snobj.rectifySED, False)
    self.assertEqual(snobj.modelOutSideTemporalRange, 'zero')
