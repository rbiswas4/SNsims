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
    def __init__(self, z, ra,dec,size, snrules, seed):
        self.seed = seed
        self.z = z
        self.ra=ra
        self.dec=dec
        self.size = size
        self.snrules=snrules

class Universe(object):
    def __init__(self, cosmology, snrules, seed):
        self.cosmology = cosmology
        self.seed = seed
        self.snrules = snrules

"""
The objective of this code is generate a set of supernovae that occur in a simulation.
The set is dercribed by a list of ntuples containing the model and its model parameters.
Reproducibility is a requirement: given a universe instance, we want to retrieve the same supernovae
independent of how the the specfic dates and solid angle used in a query. 

The model is a description of the astronomical source, its location, underlying source, and
foregrounds (dust).  sncosmo has different sources for supernovae.  Each one can have its own rules
for how they populate the universe, i.e. the pdf of their model parameters.
In addition, the realization of supernovae depends on if the simulation is galaxy- or
universe-based:

1. Given a galaxy what is the occurance of supernovae.
2. Given a universe what is the occurance of supernovae.

The generation of supernovae through galaxies is based on rates in SNU (per galaxy) units and coordinates
based relative to the galaxy.  The generation of the supernovae is based on dN/dz/dOmega (e.g. per
Mpc^3) and is not tied to galaxies.  I do not see how to make generated supernovae in a 
single universe instantiation be identical in galaxy- and universe-based realizations.  Therefore,
these two modes of generation are distict.

Each galaxy/universe must know how supernovae occur within it; the different models and for each
model the realization of objects of that model type.  Galaxy.snrules and Universe.snrules are the
object variables that contain this information, and is a list of realization rules for 
different model/generation conditions.  Each realization rule is an object.

Each realization rule is an object.  In the instantiation presented here, it is an object of
type ModelRealizer.  It is not properly designed but represents a sketch of a potential direction
for a design that is extensible.  The class has some salient features:

    - It contains a number of methods that realize specific kinds of model parameters.  Many models
    share common parameters and parameter pdfs.  This class therefore offers ready-implemented methods
    that can be of use in making rules for a particular model.

    - It contains two methods getSNe_universe and getSNe_galaxies. Generically, there are model
    parameters and source parameters. The model parameters for universe and galaxy-based simulations are
    different and hence get their own method.  However, the pattern is independent of source.

        - Realize the number of supernovae
        - Realize the angular coordinates
        - Realize the date (of explosion/maximum)
        - Realize the redshift
        - Realize the foregrounds (not implemented)
        - Realize the source parameters

    - It splits the universe into slices in space and time.  The random realization seeds are
    specified for each slice.  Queries generate SNe only for the relevent slices.  This is how
    consistent SNe are realized independent of query.  It is similar to Rahul's design except that
    the number of objects per slice is realized from a Poisson distribution, and as such the
    slices can be made small.  The implemented spatial slices are via ra, however I think eventually
    dividing the sky by healpix would be better.

    - getters for commonly used models can be provided.  For example one for SALT2.

    - units conventions
        - angles (radians)
        - time (days)

"""


class ModelRealizer(object):
    """docstring for ModelRealizerForUniverse"""
    def __init__(self, model):
        super(ModelRealizer, self).__init__()
        self.model = model
        self.rate_vol = 0.25e-4/20/365.25 # per Mpc/yr
        self.rate_snu = 0.05/365.25/20 # 2 per century distributed over all self.types
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

    def _nsnePerSlice_universe(self,universe,tindex,raindex):
        random_index = 0
        #the reproducible seed that controls generation
        self._rs.seed(np.array([universe.seed,tindex,raindex,random_index]))
        return np.random.poisson(self.rate_vol*self._rabin/2/np.pi* self._tbin*
            self.znorm(universe.cosmology))

    def _nsnePerSlice_galaxy(self,galaxy,tindex):
        random_index = 0
        #the reproducible seed that controls generation
        self._rs.seed(np.array([galaxy.seed,tindex,random_index]))
        return np.random.poisson(self.rate_snu*self._tbin)

    def _t0PerSlice(self, host,tindex, raindex=None, nsne=1):
        random_index=1
        if raindex is None:
            self._rs.seed(np.array([host.seed,tindex,random_index]))
        else:
            self._rs.seed(np.array([host.seed,tindex,raindex,random_index]))
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

    def _zPerSlice_galaxy(self, galaxy,tindex, nsne=1, subset='default'):
        if subset == 'default':
            calcind = xrange(0,nsne)
        else:
            calcind = np.where(subset)[0]

        random_index=2
        return np.zeros(nsne)+galaxy.z


    def _coordsPerSlice_universe(self, universe,tindex, raindex, nsne=1):
        random_index=3
        self._rs.seed(np.array([universe.seed,tindex,raindex,random_index]))
        phi = np.random.uniform(size=nsne,low = -np.pi/2, high = np.pi/2 )
        theta = np.arccos(np.random.uniform(size=nsne,low = np.cos(raindex*self._rabin),
            high=np.cos((raindex+1)*self._rabin)))
        return (theta,phi)

    def _coordsPerSlice_galaxy(self, galaxy,tindex, nsne=1):
        random_index=3
        self._rs.seed(np.array([galaxy.seed,tindex,random_index]))
        theta_i = 90-np.random.exponential(scale=galaxy.size)
        phi_i = np.random.uniform(low = 0, high=1)*360.
        # need to work out coordinates particularly unit conventions
        ans=astropy.modeling.rotations.RotateNative2Celestial.evaluate(phi_i,theta_i,galaxy.ra*180/np.pi,galaxy.dec*180/np.pi,0)
        return (np.array([ans[0]*np.pi/180]),np.array([ans[1]*np.pi/180]))

    def _modelParametersPerSlice(self, host,tindex, raindex=None, nsne=1, subset='default'):
        random_index=4
        if raindex is None:
            self._rs.seed(np.array([host.seed,tindex,random_index]))
        else:
            self._rs.seed(np.array([host.seed,tindex,raindex,random_index]))

        if subset == 'default':
            calcind = xrange(0,nsne)
        else:
            calcind = np.where(subset)[0]
        ans=np.zeros(nsne)
        ans[calcind] = np.random.normal(size=nsne,loc=self.M0,scale=self.M0_sig)[calcind]
        return ans

    def _dustPerSlice(self,universe,tindex,rindex):
        random_index=5
        self._rs.seed(np.array([universe.seed,tindex,raindex,random_index]))
        raise NotImplementedError()

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
        allcoords=[[],[]]
        allt0=[]
        allmp=[]

        #this series of loops can be parallelized
        for tindex in tindeces:
            for raindex in raindeces:
                nsne = self._nsnePerSlice_universe(universe,tindex,raindex)

                if nsne != 0:
                    coords = self._coordsPerSlice_universe(universe,tindex,raindex,nsne=nsne)
                    t0 = self._t0PerSlice(universe,tindex,raindex,nsne=nsne)

                    inrange = t0 >= trange[0]
                    inrange = np.logical_and(inrange,t0< trange[1])
                    inrange = np.logical_and(inrange,coords[0] < rarange[0])
                    inrange = np.logical_and(inrange,coords[0] < rarange[1])
                    inrange = np.logical_and(inrange,coords[1] < decrange[0])
                    inrange = np.logical_and(inrange,coords[1] < decrange[1])

                    if inrange.sum() > 0:
                        zs = self._zPerSlice_universe(universe,tindex,raindex,nsne=nsne,subset=inrange)
                        mp = self._modelParametersPerSlice(universe,tindex,raindex,nsne=nsne, subset=inrange)

                        allzs.extend(zs[inrange])
                        allcoords[0].extend(coords[0][inrange])
                        allcoords[1].extend(coords[1][inrange])
                        allt0.extend(t0[inrange])
                        allmp.extend(mp[inrange])

        return np.array(allzs), np.transpose(allcoords), np.array(allt0), np.array(allmp)

    def getSNe_galaxy(self, galaxy, trange):
        tindeces, raindeces = self._getSliceIndeces(trange)
        allzs=[]
        allcoords=[[],[]]
        allt0=[]
        allmp=[]

        #this series of loops can be parallelized

        for tindex in tindeces:
            nsne = self._nsnePerSlice_galaxy(galaxy,tindex)

            if nsne != 0:
                coords = self._coordsPerSlice_galaxy(galaxy,tindex,nsne=nsne)
                t0 = self._t0PerSlice(galaxy,tindex,nsne=nsne)

                inrange = t0 >= trange[0]
                inrange = np.logical_and(inrange,t0< trange[1])

                if inrange.sum() > 0:
                    zs = self._zPerSlice_galaxy(galaxy,tindex,nsne=nsne,subset=inrange)
                    mp = self._modelParametersPerSlice(galaxy,tindex,nsne=nsne, subset=inrange)

                    allzs.extend(zs[inrange])
                    allcoords[0].extend(coords[0][inrange])
                    allcoords[1].extend(coords[1][inrange])
                    allt0.extend(t0[inrange])
                    allmp.extend(mp[inrange])

        return np.array(allzs), np.transpose(allcoords), np.array(allt0), np.array(allmp)


"""
Now we can make a list of these rules for each of the models.  In the example below, I
take all the sncosmo.TimeSeriesSource of type SN IIL, SN IIP, SN Ib, SN Ib/c, SN Ic
that are in the sncosmo registry and create a ModelRealizer for each.  In this example
the behavior of each is identical but in general the internal variables for each could be
different.
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

"""
Make a universe and print out the types and the parameters for those types
"""

universe = Universe(astropy.cosmology.WMAP9, models,0)
trange = np.array([100*365.25,100.5*365.25])
rarange = np.array([10,11])*np.pi/180
decrange = np.array([10,15])*np.pi/180
for model in universe.snrules:
    print model.model.source.name
    zs,coords,t0,mp = models[0].getSNe_universe(universe,trange, rarange,decrange)
    for a,b,c,d in zip(zs,coords,t0,mp):
        print a,b,c,d

"""
Make some galaxies and print out the type and parameters for each
"""
gals = []
for i in xrange(1000):
    gals.append(Galaxy(0.5, 0, 0, 2/3600./180.*np.pi,models,i))

for gal in gals:
    for model in universe.snrules:
        zs,coords,t0,mp = models[0].getSNe_galaxy(gal, trange)
        if len(zs) !=0:
            print gal.seed, model.model.source.name
            for a,b,c,d in zip(zs,coords,t0,mp):
                print a,b,c,d

