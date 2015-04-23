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

# class ModelRules(object):
#       """docstring for ModelRules"""
#       def __init__(self, model, modelpdf, nsne, t0, z, coords):
#           self.model = model            # sncosmo model
#           self.modelpdf=modelpdf        # realization of source parameters 
#           self.nsne = nsne              # realization of number of supernovae
#           self.t0=t0                    # realization of date of explosion
#           self.z = z                    # realization of redshift
#           self.coords = coords          # realization of ra and dec

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

# class GalaxyRules(object):
#     """docstring for GalaxySetup"""
#     def __init__(self):
        
#         self.start_mjd= 0
#         self.end_mjd=10*365.25

#         self.types= ['SN IIL', 'SN IIP','SN Ib', 'SN Ib/c', 'SN Ic']
#         self.nmodels = 0
#         for t in self.types:
#             self.nmodels += len(Sources.registry_sources_as_models(t, subclassName="sncosmo.TimeSeriesSource"))

#         self.rate = 0.05/365.25/self.nmodels # 2 per century distributed over all self.types

#     def nsne(self,galaxy):
#         return np.random.poisson(self.rate * (self.end_mjd-self.start_mjd)/(1+galaxy.z))

#     def t0(self):
#         return (np.random.uniform, {'low':self.start_mjd, 'high':self.end_mjd})

#     def z(self, galaxy):
#         return galaxy.z

#     def coords_gal(self, galaxy):
#         theta_i = 90-np.random.exponential(scale=galaxy.size)
#         phi_i = np.random.uniform(low = 0, high=1)*360.
#         # need to work out coordinates particularly unit conventions
#         return astropy.modeling.rotations.RotateNative2Celestial.evaluate(phi_i,theta_i,galaxy.ra,galaxy.dec,0)

#     def modelpdf(self):
#         return (np.random.normal,{'loc':-19.3,'scale':0.3})

#     #For all sncosmo model types, make an associated ModelRules
#     def snRules(self):
#         ## In this simple example all SN types have the same rules for generation
#         sneInAGalaxy = []
#         for t in self.types:
#             for model in Sources.registry_sources_as_models(t, subclassName="sncosmo.TimeSeriesSource"):
#                 sneInAGalaxy.append(
#                     ModelRules(model, self.modelpdf, self.nsne, self.t0, self.z, self.coords_gal))
#         return sneInAGalaxy


"""
Now when we make a galaxy we specify rules for the occurance of the SN types.
Here we assume SNe go off in all galaxies the same way.
"""

# gals = []
# rules = GalaxyRules() # only one of these assuming all galaxies have the same properties

# for i in xrange(1000):
#     gals.append(Galaxy(0.5, 0, 0, 2/3600.,rules.snRules(),i))

"""
Introduce a class that contains the different ways that a SN can be generated.
"""

# class SNGenerator(object):

#     # a vectorized version of this will improve performance
#     @staticmethod
#     def sneInGalaxies(galaxies):
#         out = []
#         for  galaxy in galaxies:
#             np.random.seed(seed=galaxy.seed)
#             for info in galaxy.snrule:
#                 #for each type determine whether there is a SN
#                 nsn=info.nsne(galaxy)
#                 for i in xrange(nsn):
#                     # galaxy resdshift is first parameter
#                     pars = [info.z(galaxy)]

#                     pars.extend(info.coords(galaxy))

#                     # date of explosion is second parameter
#                     pars.append(info.t0()[0](**info.t0()[1]))

#                     # parameters of the target
#                     pars.append(info.modelpdf()[0](**info.modelpdf()[1]))

#                     out.append((info.model,pars))
#         return out

#     @staticmethod
#     def sneInUniverse(universe):
#         out=[]
#         np.random.seed(seed=universe.seed)
#         snRule = universe.snrule
#         for info in snRule:
#             #for each type determine whether there is a SN
#             nsn=info.nsne(universe)
#             for i in xrange(nsn):
#                 # galaxy resdshift is first parameter
#                 pars = [info.z(universe)]

#                 pars.extend(info.coords())

#                 # date of explosion is second parameter
#                 pars.append(info.t0()[0](**info.t0()[1]))

#                 # parameters of the target
#                 pars.append(info.modelpdf()[0](**info.modelpdf()[1]))
#                 print pars
#                 out.append((info.model,pars))

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

# class UniverseRules(object):
#     """docstring for GalaxySetup"""
#     def __init__(self):

#         self.zmin=0.03
#         self.zmax =1.2

#         self.start_mjd= 0
#         self.end_mjd=1*365.25

#         #although outputs are in degrees inputs in radians
#         self.rarange=np.array([10,20])*np.pi/180
#         self.decrange =  np.array([10,20])*np.pi/180

#         # self.solidangle=self.solidAngle()

#         self.types= ['SN IIL', 'SN IIP','SN Ib', 'SN Ib/c', 'SN Ic']
#         self.nmodels = 0
#         for t in self.types:
#             self.nmodels += len(Sources.registry_sources_as_models(t, subclassName="sncosmo.TimeSeriesSource"))

#         self.rate = 0.25e-4/self.nmodels # per Mpc/yrself.

#         self._znorm=dict()

#         # variables that control the realization of supernovae into time and spatial slices
#         # to allow efficient reproducable SN retrieval for the same universe 
#         self._tbin = 10. // years
#         self._rabin = np.pi/8.

#     # def solidAngle(self):
#     #     return (np.cos(np.pi/2-self.decrange[1])-
#     #         np.cos(np.pi/2-self.decrange[0]))*(self.rarange[1]-self.rarange[0])

    
#     def _nsnePerSlice(self,universe,tindex,raindex):
#         random_index = 0
#         #the reproducible seed that controls generation
#         numpy.random.RandomState.seed([universe.seed,tindex,raindex,random_index])
#         return np.random.poisson(self.rate*self._rabin/2/np.pi* self._tbin*
#             self.znorm(universe.cosmology))

#     def _t0PerSlice(self, universe,tindex raindex):
#         random_index=1
#         numpy.random.RandomState.seed([universe.seed,tindex,raindex,random_index])
#         return (np.random.uniform, {'low':tindex*self._tbin, 'high':self.(tindex+1)*self._tbin})

#     def _zPerSlice(self, universe,tindex raindex):
#         random_index=2
#         numpy.random.RandomState.seed([universe.seed,tindex,raindex,random_index])
#         ans = np.random.uniform(low = 0, high=1)
#         ans = scipy.optimize.newton(lambda zp: scipy.integrate.quad(lambda z: 
#             universe.cosmology.differential_comoving_volume(z).value/(1+z),self.zmin,zp)[0]/self.znorm(universe.cosmology) - ans,0.5)
#         return ans

#     def _coordsPerSlice(self, universe,tindex raindex):
#         random_index=3
#         numpy.random.RandomState.seed([universe.seed,tindex,raindex,random_index])
#         phi = np.random.uniform(low = -90, high = 90 )
#         theta = np.arccos(np.random.uniform(low = raindex*self._rabin,
#             high=(raindex+1)*self._rabin))*180/np.pi
#         return (theta,phi)

#     def _getSliceIndeces(self):
#         tindeces = numpy.arange(int(self.start_mjd/self._tbin/365.25),int(self.end_mjd/self._tbin/365.25),
#             self._tbin,dtype='int')
#         raindeces = numpy.arange(int(self.rarange[0]/self._rabin),int(self.rarange[1]/self._rabin),
#             self._rabin,dtype='int')
#         return tindeces,raindeces

#     def getSNe(self, universe):
#         out=[]
#         for tindex in tindeces:
#             for raindex in raindeces:
#                 pars = [info.z(universe)]

#                 pars.extend(info.coords())

#                 # date of explosion is second parameter
#                 pars.append(info.t0()[0](**info.t0()[1]))

#                 # parameters of the target
#                 pars.append(info.modelpdf()[0](**info.modelpdf()[1]))
#                 print pars
#                 out.append((info.model,pars))


#     def nsne(self, universe):
#         return np.random.poisson(self.rate*self.solidangle/4/np.pi* (self.end_mjd-self.start_mjd)*
#             self.znorm(universe.cosmology))

#     def t0(self, universe):
#         return (np.random.uniform, {'low':self.start_mjd, 'high':self.end_mjd})

#     def znorm(self,cosmology):
#         if cosmology not in self._znorm:
#             self._znorm[cosmology] = scipy.integrate.quad(
#                 lambda z: cosmology.differential_comoving_volume(z).value/(1+z),self.zmin,self.zmax)[0]
#         return self._znorm[cosmology]

#     def z(self,universe):
#         ans = np.random.uniform(low = 0, high=1)
#         ans = scipy.optimize.newton(lambda zp: scipy.integrate.quad(lambda z: 
#             universe.cosmology.differential_comoving_volume(z).value/(1+z),self.zmin,zp)[0]/self.znorm(universe.cosmology) - ans,0.5)
#         return ans

#     def coords(self):
#         phi = np.random.uniform(low = self.decrange[0], high=self.decrange[1])*180/np.pi
#         theta = np.arccos(np.random.uniform(low = self.rarange[0], high=self.rarange[1]))*180/np.pi
#         return (theta,phi)

#     def modelpdf(self):
#         return (np.random.normal,{'loc':-19.3,'scale':0.3})

#     #For all sncosmo model types, make an associated ModelRules
#     def snRules(self):
#         ## In this simple example all sne have the same properties
#         sneInAUniverse = []
#         for t in self.types:
#             for model in Sources.registry_sources_as_models(t, subclassName="sncosmo.TimeSeriesSource"):
#                 sneInAUniverse.append(
#                     ModelRules(model, self.modelpdf, self.nsne, self.t0, self.z,
#                     self.coords))
#         return sneInAUniverse

class ModelRealizer(object):
    """docstring for ModelRealizerForUniverse"""
    def __init__(self, model):
        super(ModelRealizer, self).__init__()
        self.model = model
        self.rate = 0.25e-4/10/365.25 # per Mpc/yr
        self._znorm=dict()

        self.zmin=0.03
        self.zmax =1.2


        #model mean and sigma magnitude
        self.M0=-19.3
        self.M0_sig = 0.3

        # variables that control the realization of supernovae into time and spatial slices
        # to allow efficient reproducable SN retrieval for the same universe 
        self._tbin = 365.25 # years
        self._rabin = np.pi/64.

        self._rs  = np.random.RandomState()
        
    def _nsnePerSlice(self,universe,tindex,raindex):
        random_index = 0
        #the reproducible seed that controls generation
        self._rs.seed(np.array([universe.seed,tindex,raindex,random_index]))
        return np.random.poisson(self.rate*self._rabin/2/np.pi* self._tbin*
            self.znorm(universe.cosmology))

    def _t0PerSlice(self, universe,tindex, raindex, nsne=1):
        random_index=1
        self._rs.seed(np.array([universe.seed,tindex,raindex,random_index]))
        return np.random.uniform(size=nsne,low=tindex*self._tbin, high=(tindex+1)*self._tbin)

    def _zPerSlice_universe(self, universe,tindex, raindex, nsne=1, subset='default'):
        if subset == 'default':
            calcind = xrange(0,nsne)
        else:
            calcind = np.where(subset)[0]

        random_index=2
        self._rs.seed(np.array([universe.seed,tindex,raindex,random_index]))
        anss = np.random.uniform(size=nsne,low = 0, high=1)

        for i in xrange(len(calcind)):
            anss[calcind[i]]=scipy.optimize.newton(lambda zp: scipy.integrate.quad(lambda z: 
            universe.cosmology.differential_comoving_volume(z).value/(1+z),self.zmin,zp)[0]/self.znorm(universe.cosmology) - anss[calcind[i]],
            0.5,fprime=lambda z:universe.cosmology.differential_comoving_volume(z).value/(1+z)/self.znorm(universe.cosmology))
        return anss

    def _coordsPerSlice_universe(self, universe,tindex, raindex, nsne=1):
        random_index=3
        self._rs.seed(np.array([universe.seed,tindex,raindex,random_index]))
        phi = np.random.uniform(size=nsne,low = -np.pi/2, high = np.pi/2 )
        theta = np.arccos(np.random.uniform(size=nsne,low = np.cos(raindex*self._rabin),
            high=np.cos((raindex+1)*self._rabin)))
        return (theta,phi)

    def _coordsPerSlice_galaxy(self, galaxy,tindex, nsne=1):
        random_index=3
        self._rs.seed(np.array([universe.seed,tindex,random_index]))
        theta_i = 90-np.random.exponential(scale=galaxy.size)
        phi_i = np.random.uniform(low = 0, high=1)*360.
        # need to work out coordinates particularly unit conventions
        return np.pi/180*astropy.modeling.rotations.RotateNative2Celestial.evaluate(phi_i,theta_i,galaxy.ra,galaxy.dec,0)

    def _modelParametersPerSlice(self, universe,tindex, raindex, nsne=1, subset='default'):
        if subset == 'default':
            calcind = xrange(0,nsne)
        else:
            calcind = np.where(subset)[0]
        random_index=4
        self._rs.seed(np.array([universe.seed,tindex,raindex,random_index]))
        ans=np.zeros(nsne)
        ans[calcind] = np.random.normal(size=nsne,loc=self.M0,scale=self.M0_sig)[calcind]
        return ans

    def _getSliceIndeces(self, trange, rarange=None):
        tmin = int(trange[0]/self._tbin)
        tmax = max(np.ceil(trange[1]/self._tbin),tmin+1)
        tindeces = np.arange(tmin,tmax,dtype='int')

        raindeces=None
        if rarange is not None:
            ramin = int(rarange[0]/self._rabin)
            ramax = max(np.ceil(rarange[1]/self._rabin),ramin+1)
            raindeces = np.arange(ramin,ramax,dtype='int')
        return tindeces,raindeces


    def znorm(self,cosmology):
        if cosmology not in self._znorm:
            self._znorm[cosmology] = scipy.integrate.quad(
                lambda z: cosmology.differential_comoving_volume(z).value/(1+z),self.zmin,self.zmax)[0]
        return self._znorm[cosmology]

    def getSNe_universe(self, universe, trange, rarange, decrange):
        tindeces, raindeces = self._getSliceIndeces(trange, rarange)
        allzs=[]
        allcoords=[]
        allt0=[]
        allmp=[]

        #this can be parallelized
        for tindex in tindeces:
            for raindex in raindeces:
                nsne = self._nsnePerSlice(universe,tindex,raindex)

                if nsne != 0:
                    coords = self._coordsPerSlice_universe(universe,tindex,raindex,nsne=nsne)
                    t0 = self._t0PerSlice(universe,tindex,raindex,nsne=nsne)

                    inrange = t0 >= trange[0]
                    inrange = np.logical_and(inrange,t0< trange[1])
                    inrange = np.logical_and(inrange,coords[0] < rarange[0])
                    inrange = np.logical_and(inrange,coords[0] < rarange[1])
                    inrange = np.logical_and(inrange,coords[1] < decrange[0])
                    inrange = np.logical_and(inrange,coords[1] < decrange[1])

                    zs = self._zPerSlice_universe(universe,tindex,raindex,nsne=nsne,subset=inrange)
                    mp = self._modelParametersPerSlice(universe,tindex,raindex,nsne=nsne, subset=inrange)

                    allzs.append(zs[inrange])
                    allcoords.append((coords[0][inrange], coords[1][inrange]))
                    allt0.append(t0[inrange])
                    allmp.append(mp[inrange])


        return allzs, allcoords, allt0, allmp

def getSNe_galaxy(self, galaxy, trange, rarange, decrange):
        tindeces, = self._getSliceIndeces(trange)
        allzs=[]
        allcoords=[]
        allt0=[]
        allmp=[]

        #this can be parallelized
        for tindex in tindeces:
            for raindex in raindeces:
                nsne = self._nsnePerSlice(universe,tindex,raindex)
                if nsne != 0:
                    zs = self._zPerSlice_universe(universe,tindex,raindex,nsne=nsne)
                    coords = self._coordsPerSlice_universe(universe,tindex,raindex,nsne=nsne)

                    t0 = self._t0PerSlice(universe,tindex,raindex,nsne=nsne)

                    mp = self._modelParametersPerSlice(universe,tindex,raindex,nsne=nsne)


                    inrange = t0 >= trange[0]
                    inrange = np.logical_and(inrange,t0< trange[1])
                    inrange = np.logical_and(inrange,coords[0] < rarange[0])
                    inrange = np.logical_and(inrange,coords[0] < rarange[1])
                    inrange = np.logical_and(inrange,coords[1] < decrange[0])
                    inrange = np.logical_and(inrange,coords[1] < decrange[1])

                    allzs.append(zs[inrange])
                    allcoords.append((coords[0][inrange], coords[1][inrange]))
                    allt0.append(t0[inrange])
                    allmp.append(mp[inrange])


        return allzs, allcoords, allt0, allmp

    # def t0(self, universe):
    #     return (np.random.uniform, {'low':self.start_mjd, 'high':self.end_mjd})



    # def z(self,universe):
    #     ans = np.random.uniform(low = 0, high=1)
    #     ans = scipy.optimize.newton(lambda zp: scipy.integrate.quad(lambda z: 
    #         universe.cosmology.differential_comoving_volume(z).value/(1+z),self.zmin,zp)[0]/self.znorm(universe.cosmology) - ans,0.5)
    #     return ans

    # def coords(self):
    #     phi = np.random.uniform(low = self.decrange[0], high=self.decrange[1])*180/np.pi
    #     theta = np.arccos(np.random.uniform(low = self.rarange[0], high=self.rarange[1]))*180/np.pi
    #     return (theta,phi)




"""
Given the rules for how supernovae go off un a universe, generate supernovae
"""

def modelsInUniverse():
    ## In this simple example all sne have the same properties
    sneInAUniverse = []
    types= ['SN IIL', 'SN IIP','SN Ib', 'SN Ib/c', 'SN Ic']
    for t in types:
        for model in Sources.registry_sources_as_models(t, subclassName="sncosmo.TimeSeriesSource"):
            sneInAUniverse.append(ModelRealizer(model))
    return sneInAUniverse

models = modelsInUniverse()
universe = Universe(astropy.cosmology.WMAP9, models,0)
ans = models[0].getSNe_universe(universe,np.array([100*365.25,101.5*365.25]),np.array([10,12])*np.pi/180,
    np.array([10,20])*np.pi/180)



wefwe
