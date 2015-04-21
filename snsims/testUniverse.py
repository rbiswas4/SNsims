#!/usr/bin/env python


import sncosmo
import numpy as np
import scipy
import Sources
import astropy.modeling
import astropy.cosmology

"""
Introduce Galaxy and Universe classes that may exist elsewhere
already but are included here for the demonstration of the 
proposed usage.  Note that the galaxy and universe are responsible
for knowing how supernovae explide in them.
"""
class Galaxy(object):
    """docstring for Galaxy"""
    def __init__(self, z, ra,dec,size, snrule, seed):
        self.seed = seed
        self.z = z
        self.ra=ra
        self.dec=dec
        self.size = size
        self.snrule=snrule

class Universe(object):
    def __init__(self, cosmology, snrule, seed):
        self.cosmology = cosmology
        self.seed = seed
        self.snrule = snrule

"""
The objective of this code is generate a set of supernovae that occur in a simulation.
The set is dercribed by a list of ntuples containing the model and its model parameters.
"""

"""
sncosmo has different models for supernovae.  Each one can have its own rules
for how they populate the universe, i.e. the pdf of their model parameters.
This class contain those rules.  The machine that generates supernovae expects and
knows how to use this class
"""

class ModelRules(object):
      """docstring for ModelRules"""
      def __init__(self, model, modelpdf, nsne, t0, z, coords):
          self.model = model            # sncosmo model
          self.modelpdf=modelpdf        # realization of source parameters 
          self.nsne = nsne              # realization of number of supernovae
          self.t0=t0                    # realization of date of explosion
          self.z = z                    # realization of redshift
          self.coords = coords          # realization of ra and dec

"""
There are two ways we anticipate realizing supernovae.

1. Given a galaxy what is the occurance of supernovae.
2. Given a universe what is the occurance of supernovae.

The rules of realization of model parameters are different depending on
whether the generation is galaxy or universe based.
"""

""" 
First consider galaxy-based generation.  The rules for the occurance of supernovae
are governed by galaxy propreties.  So for each galaxy give a description of what
things explode inside of it.

This container contains the information needed for specifying galaxy-based
SN generation.

This class is not meant to be reused, just an example here, though further
design could make it useful for other uses.
"""

class GalaxyRules(object):
    """docstring for GalaxySetup"""
    def __init__(self):
        
        self.start_mjd= 0
        self.end_mjd=10*365.25

        self.types= ['SN IIL', 'SN IIP','SN Ib', 'SN Ib/c', 'SN Ic']
        self.nmodels = 0
        for t in self.types:
            self.nmodels += len(Sources.registry_sources_as_models(t, subclassName="sncosmo.TimeSeriesSource"))

        self.rate = 0.05/365.25/self.nmodels # 2 per century distributed over all self.types

    def nsne(self,galaxy):
        return np.random.poisson(self.rate * (self.end_mjd-self.start_mjd)/(1+galaxy.z))

    def t0(self):
        return (np.random.uniform, {'low':self.start_mjd, 'high':self.end_mjd})

    def z(self, galaxy):
        return galaxy.z

    def coords_gal(self, galaxy):
        theta_i = 90-np.random.exponential(scale=galaxy.size)
        phi_i = np.random.uniform(low = 0, high=1)*360.
        # need to work out coordinates particularly unit conventions
        return astropy.modeling.rotations.RotateNative2Celestial.evaluate(phi_i,theta_i,galaxy.ra,galaxy.dec,0)

    def modelpdf(self):
        return (np.random.normal,{'loc':-19.3,'scale':0.3})

    #For all sncosmo model types, make an associated ModelRules
    def snRules(self):
        ## In this simple example all SN types have the same rules for generation
        sneInAGalaxy = []
        for t in self.types:
            for model in Sources.registry_sources_as_models(t, subclassName="sncosmo.TimeSeriesSource"):
                sneInAGalaxy.append(
                    ModelRules(model, self.modelpdf, self.nsne, self.t0, self.z, self.coords_gal))
        return sneInAGalaxy


"""
Now when we make a galaxy we specify rules for the occurance of the SN types.
Here we assume SNe go off in all galaxies the same way.
"""

gals = []
rules = GalaxyRules() # only one of these assuming all galaxies have the same properties

for i in xrange(1000):
    gals.append(Galaxy(0.5, 0, 0, 2/3600.,rules.snRules(),i))

"""
Introduce a class that contains the different ways that a SN can be generated.
"""

class SNGenerator(object):

    # a vectorized version of this will improve performance
    @staticmethod
    def sneInGalaxies(galaxies):
        out = []
        for  galaxy in galaxies:
            np.random.seed(seed=galaxy.seed)
            for info in galaxy.snrule:
                #for each type determine whether there is a SN
                nsn=info.nsne(galaxy)
                for i in xrange(nsn):
                    # galaxy resdshift is first parameter
                    pars = [info.z(galaxy)]

                    pars.extend(info.coords(galaxy))

                    # date of explosion is second parameter
                    pars.append(info.t0()[0](**info.t0()[1]))

                    # parameters of the target
                    pars.append(info.modelpdf()[0](**info.modelpdf()[1]))

                    out.append((info.model,pars))
        return out

    @staticmethod
    def sneInUniverse(universe):
        out=[]
        np.random.seed(seed=universe.seed)
        snRule = universe.snrule
        for info in snRule:
            #for each type determine whether there is a SN
            nsn=info.nsne(universe)
            for i in xrange(nsn):
                # galaxy resdshift is first parameter
                pars = [info.z(universe)]

                pars.extend(info.coords())

                # date of explosion is second parameter
                pars.append(info.t0()[0](**info.t0()[1]))

                # parameters of the target
                pars.append(info.modelpdf()[0](**info.modelpdf()[1]))
                print pars
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

class UniverseRules(object):
    """docstring for GalaxySetup"""
    def __init__(self):

        self.zmin=0.03
        self.zmax =1.2

        self.start_mjd= 0
        self.end_mjd=1*365.25

        #although outputs are in degrees inputs in radians
        self.rarange=np.array([10,20])*np.pi/180
        self.decrange =  np.array([10,20])*np.pi/180

        self.solidangle=self.solidAngle()

        self.types= ['SN IIL', 'SN IIP','SN Ib', 'SN Ib/c', 'SN Ic']
        self.nmodels = 0
        for t in self.types:
            self.nmodels += len(Sources.registry_sources_as_models(t, subclassName="sncosmo.TimeSeriesSource"))

        self.rate = 0.25e-4/self.nmodels # per Mpc/yrself.

        self._znorm=dict()

    def solidAngle(self):
        return (np.cos(np.pi/2-self.decrange[1])-
            np.cos(np.pi/2-self.decrange[0]))*(self.rarange[1]-self.rarange[0])

    def nsne(self, universe):
        return np.random.poisson(self.rate*self.solidangle/4/np.pi* (self.end_mjd-self.start_mjd)*
            self.znorm(universe.cosmology))

    def t0(self):
        return (np.random.uniform, {'low':self.start_mjd, 'high':self.end_mjd})

    def znorm(self,cosmology):
        if cosmology not in self._znorm:
            self._znorm[cosmology] = scipy.integrate.quad(
                lambda z: cosmology.differential_comoving_volume(z).value/(1+z),self.zmin,self.zmax)[0]
        return self._znorm[cosmology]

    def z(self,universe):
        ans = np.random.uniform(low = 0, high=1)
        ans = scipy.optimize.newton(lambda zp: scipy.integrate.quad(lambda z: 
            universe.cosmology.differential_comoving_volume(z).value/(1+z),self.zmin,zp)[0]/self.znorm(universe.cosmology) - ans,0.5)
        return ans

    def coords(self):
        phi = np.random.uniform(low = self.decrange[0], high=self.decrange[1])*180/np.pi
        theta = np.arccos(np.random.uniform(low = self.rarange[0], high=self.rarange[1]))*180/np.pi
        return (theta,phi)

    def modelpdf(self):
        return (np.random.normal,{'loc':-19.3,'scale':0.3})

    #For all sncosmo model types, make an associated ModelRules
    def snRules(self):
        ## In this simple example all sne have the same properties
        sneInAUniverse = []
        for t in self.types:
            for model in Sources.registry_sources_as_models(t, subclassName="sncosmo.TimeSeriesSource"):
                sneInAUniverse.append(
                    ModelRules(model, self.modelpdf, self.nsne, self.t0, self.z,
                    self.coords))
        return sneInAUniverse

"""
Given the rules for how supernovae go off un a universe, generate supernovae
"""

rules = UniverseRules()
universe = Universe(astropy.cosmology.WMAP9, rules.snRules(),0)
sne= SNGenerator.sneInUniverse(universe)

