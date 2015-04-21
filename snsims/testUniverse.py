#!/usr/bin/env python


import sncosmo
import numpy as np
import scipy
import Sources
import astropy.modeling
import astropy.cosmology

"""
Each SN model has a rule for where/when it explodes and for the distribution
of internal model parameters
"""

class ModelInfo(object):
      """docstring for ModelInfo"""
      def __init__(self, model, modelpdf, nsne, t0, z, coords):
          self.model = model            #sncosmo model
          self.modelpdf=modelpdf                  #satifies numpy.random distribution 
          self.nsne = nsne
          self.t0=t0
          self.z = z
          self.coords = coords

""" 
Consider realizing SNe through galaxies.  The rules for the occurance of supernovae are governed by
galaxy propreties.
So for each galaxy give a description of what things explode inside of it.
"""

start_mjd= 0
end_mjd=100*365.25

types= ['SN IIL', 'SN IIP','SN Ib', 'SN Ib/c', 'SN Ic']
nmodels = 0
for t in types:
    nmodels += len(Sources.registry_sources_as_models(t, subclassName="sncosmo.TimeSeriesSource"))
rate = 0.05/365.25/nmodels # 2 per century distributed over all types
nsne = lambda galaxy: np.random.poisson(rate * (end_mjd-start_mjd)/(1+galaxy.z))
t0=  (np.random.uniform, {'low':start_mjd, 'high':end_mjd})
z = lambda galaxy: galaxy.z
def coords_gal(galaxy):
    theta_i = 90-np.random.exponential(scale=galaxy.size)
    phi_i = np.random.uniform(low = 0, high=1)*360.
    # need to work out coordinates particularly unit conventions
    return astropy.modeling.rotations.RotateNative2Celestial.evaluate(phi_i,theta_i,galaxy.ra,galaxy.dec,0)
modelpdf= (np.random.normal,{'loc':-19.3,'scale':0.3})


## In this simple example all sne have the same properties
sneInAGalaxy = []
for t in types:
    for model in Sources.registry_sources_as_models(t, subclassName="sncosmo.TimeSeriesSource"):
        sneInAGalaxy.append(
            ModelInfo(model, modelpdf, nsne, t0, z, coords_gal))

class Galaxy(object):
    """docstring for Galaxy"""
    def __init__(self, z, ra,dec,size, sninfo, seed):
        self.seed = seed
        self.z = z
        self.ra=ra
        self.dec=dec
        self.size = size
        self.sninfo=sninfo

"""
lets make a lot of galaxies all with the same rule for its constituant supernovae
"""

gals = []
for i in xrange(1000):
    gals.append(Galaxy(0.5, 0, 0, 2/3600.,sneInAGalaxy,i))

"""
Introduce a class that contains the different ways that a SN can be generated
"""

class SNGenerator(object):

    # a vectorized version of this will improve performance
    @staticmethod
    def sneInGalaxies(galaxies):
        out = []
        for  galaxy in galaxies:
            np.random.seed(seed=galaxy.seed)
            for info in galaxy.sninfo:
                #for each type determine whether there is a SN
                nsn=info.nsne(galaxy)
                for i in xrange(nsn):
                    # galaxy resdshift is first parameter
                    pars = [info.z(galaxy)]

                    pars.extend(info.coords(galaxy))

                    # date of explosion is second parameter
                    pars.append(info.t0[0](**info.t0[1]))

                    # parameters of the target
                    pars.append(info.modelpdf[0](**info.modelpdf[1]))

                    out.append((info.model,pars))
        return out

    @staticmethod
    def sneInUniverse(sneInUniverse, cosmology):
        out=[]
        for info in sneInUniverse:
            #for each type determine whether there is a SN
            nsn=info.nsne(cosmology)
            for i in xrange(nsn):
                # galaxy resdshift is first parameter
                pars = [info.z(cosmology)]

                pars.extend(info.coords())

                # date of explosion is second parameter
                pars.append(info.t0[0](**info.t0[1]))

                # parameters of the target
                pars.append(info.modelpdf[0](**info.modelpdf[1]))

                out.append((info.model,pars))

"""
Generate the supernovae for the galaxies
"""

# sne= SNGenerator.sneInGalaxies(gals)
# print len(sne)  
# for a,b in sne:
#     print a.source.name,b


"""
The alternative is to generate supernovae in a universe not directly associated with a galaxy.
Then the setup different models is a bit different.
"""

start_mjd= 0
end_mjd=1*365.25

zmin=0.03
zmax =1.2

types= ['SN IIL', 'SN IIP','SN Ib', 'SN Ib/c', 'SN Ic']
nmodels = 0
for t in types:
    nmodels += len(Sources.registry_sources_as_models(t, subclassName="sncosmo.TimeSeriesSource"))
rate = 0.25e-4/nmodels # per Mpc/yr

def nsne_universe(cosmology):
    return np.random.poisson(rate* 
        scipy.integrate.quad(lambda z: cosmology.differential_comoving_volume(z).value/(1+z),zmin,zmax)[0])

t0=  (np.random.uniform, {'low':start_mjd, 'high':end_mjd})
cosmology = astropy.cosmology.WMAP9
znorm = scipy.integrate.quad(lambda z: cosmology.differential_comoving_volume(z).value/(1+z),zmin,zmax)[0]

def z_universe(cosmology):
    ans = np.random.uniform(low = 0, high=1)
    ans = scipy.optimize.newton(lambda zp: scipy.integrate.quad(lambda z: 
        cosmology.differential_comoving_volume(z).value/(1+z),zmin,zp)[0].znorm - ans)
    return ans

def coords_universe():
    phi = np.random.uniform(low = 0, high=1)*360.
    theta = np.arccos(np.random.uniform(low = 0, high=np.pi))*360/np.pi
    return (theta,phi)
modelpdf= (np.random.normal,{'loc':-19.3,'scale':0.3})

## In this simple example all sne have the same properties
sneInAUniverse = []
for t in types:
    for model in Sources.registry_sources_as_models(t, subclassName="sncosmo.TimeSeriesSource"):
        sneInAUniverse.append(
            ModelInfo(model, modelpdf, nsne_universe, t0, z_universe, coords_universe))


sne= SNGenerator.sneInUniverse(sneInAUniverse, cosmology)

