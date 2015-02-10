#!/usr/bin/env python

"""
This is a wrapper class for `~sncosmo.Model` adding attributes and
functionality, perhaps at the cost of flexibility of use, if appropriate for
a particular project
"""
import sncosmo
import numpy as np
from cStringIO import StringIO
import sys

class SNObject(object) :
    '''
    Implementation of a SN with certain properties
    
    ra :

    dec :

    snmodelsoure

    dusttype

    parameter distribution

    rate



    '''
    def __init__(self, id, ra, dec, sourcemodel='salt2', hostdust, mwdust):

        self.id = id
        _snid = self.id
        self.ra = ra
        self.dec = dec 
        self.sourcemodel = sourcemodel

        Model.__init__(self, source=sourcemodel,
                       effects=[hostdust, mwdust],
                       effect_names=['host', 'mw'],
                       effect_frames=['rest', 'obs'])

        
    @property 
    def snid(self) :
        return self.id 

    @property 
    def snra (self) :
        return self.ra 


    @property 
    def sndec (self) :
        return self.dec 

    @property 
    def z (self) :
        return self._z 

    mabs = np.random.normal(-19.3, 0.3)
    if self.modelname == "SALT2exact":
        model = sncosmo.Model(source = 'salt2')
    else:
        raise ValueError("Don't know how to treat this supernova modelname")

    #model = sncosmo.model(source='SALT2')
    model.set(z=z)
    model.set_source_peakabsmag(mabs, 'bessellb', 'ab')
    x0 = model.get('x0')
    p = {'z':z, 't0':uniform(tmin, tmax), 'x0':x0, 'x1': normal(0., 1.), 'c': normal(0., 0.1)}
    params.append(p)


    def dictifySALTparams():
        return 0 

    for chunk in catalogIterator:
        if dtype is None:
            dtype = chunk.dtype
        for record in chunk:
            id = np.int(record['id'])
            mass = np.float(record['mass_stellar'])*10.0**10.0
            ra = np.float(record['ra'])
            dec = np.float(record['dec'])
            redshift = np.float(record['redshift'])
            sn = SNIa( id=id, ra=ra, dec=dec, z=redshift, mass=mass, modelname="SALT2exact") 

            print sn.snid, sn.z , sn.snra , sn.sndec



