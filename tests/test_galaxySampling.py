import snsims
from snsims import SersicSamples
import numpy as np

def test_expProfile():
   rng = np.random.RandomState(0)
   ss = SersicSamples(rng)
   s1 = ss.sampleRadius(numSamples=1000000, halfLightRadius=1.0, sersicIndex=1) 
   np.isclose(np.median(s1), 1.0, atol=0.001)
def test_deVacProfile():
   rng = np.random.RandomState(0)
   ss = SersicSamples(rng)
   s1 = ss.sampleRadius(numSamples=1000000, halfLightRadius=1.0, sersicIndex=4) 
   np.isclose(np.median(s1), 1., atol=0.001)
