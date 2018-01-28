"""
Module with concrete implementations of functions important for sampling
parameters
"""
from __future__ import absolute_import, print_function
from copy import deepcopy
import numpy as np
import pandas as pd
import healpy as hp
import sncosmo
from opsimsummary import HealpixTiles
from .populationParamSamples import (RateDistributions,
                                     SALT2Parameters,
                                     PositionSamples)
from astropy.cosmology import Planck15
from .samplingGalaxies import SersicSamples

__all__ = ['PowerLawRates', 'SimpleSALTDist', 'CoordSamples', 'TwinklesRates',
           'CatSimPositionSampling', 'TwinklesSim', 'SALT2_MMDist',
           'GMM_SALT2Params']

def double_gauss(mu, sigp, sigm, size):
    """Double Gaussian distribution. Note: mu is the mode and not the mean."""
    
    sigm = abs(sigm) # Just in case

    p = np.array([sigp, sigm], dtype=np.float64) # probability of one side is proportional to sigma on that side
    p /= sum(p)


    sig = np.random.choice([sigp, -sigm], size = size, replace = True, p = p)
    
    return abs(np.random.normal(size = size))*sig + mu

def SALT2_MMDist(numSN,
                cm=-0.0474801042369, cs1=0.0965032273527, cs2=0.042844366359,
                x1m=0.872727291354, x1s1=0.358731835038, x1s2=1.42806797468,
                mm=10.701690617, ms1=0.334359086569, ms2=1.0750402101,
                mBm=-19.0199168813, mc=-0.0838387899933, mt=10.,
                cc=3.20907949118, cx1=-0.137042055737):
    """
    Generates "noise"-free mB, x1, c. Trained on JLA SALT2-4 SDSS z < 0.2 and
    SNLS z < 0.5. Very simple model with linear alpha/beta and same
    distribution irrspective of host-mass. mB needs h=0.7 distance modulus
    added to it.

    From D. Rubin
    """
    color = double_gauss(cm, cs1, cs2, size=numSN)
    x1 = double_gauss(x1m, x1s1, x1s2, size=numSN)
    mass = double_gauss(mm, ms1, ms2, size=numSN)

    mB = mBm + mc * (mass > 10.) + cc * color + cx1 * x1

    return mB, x1, color, mass

class SimpleSALTDist(SALT2Parameters):
    """
    Concrete Implementation of `SALT2Parameters`
    """
    def __init__(self, numSN, zSamples, snids=None, alpha=0.11, beta=3.14,
                 cSigma=0.1, x1Sigma=1.0, meanM=-19.3, Mdisp=0.15, rng=None,
                 cosmo=Planck15, mjdmin=0., surveyDuration=10.):
        """
        """
        self._snids = snids
        self.alpha = alpha
        self.beta = beta
        self._numSN = numSN
        self.zSamples = zSamples
        self.x1Sigma = x1Sigma
        self.cSigma = cSigma
        self._rng = rng
        self.centralMabs = meanM
        self.Mdisp = Mdisp
        self._paramSamples = None
        self.cosmo = cosmo
        self.mjdmin = mjdmin
        self.surveyDuration = surveyDuration

    @property
    def mjdmax(self):
        return self.mjdmin + self.surveyDuration * 365.0

    @property
    def snids(self):
        if self._snids is None:
            self._snids = np.arange(self.numSN)
        return self._snids

    @property
    def numSN(self):
        if self._numSN is None:
            if self.zSamples is None:
                raise ValueError('Both zSamples and numSN cannot be None')
            self._numSN = len(self.zSamples)
        return self._numSN

    @property
    def randomState(self):
        if self._rng is None:
            raise NotImplementedError('rng must be provided')
        return self._rng

    @property
    def paramSamples(self):
        if self._paramSamples is None:
            timescale = self.mjdmax - self.mjdmin
            T0Vals = self.randomState.uniform(size=self.numSN) * timescale \
                    + self.mjdmin
            cvals = self.randomState.normal(loc=0., scale=self.cSigma,
                                            size=self.numSN)
            x1vals = self.randomState.normal(loc=0., scale=self.x1Sigma,
                                            size=self.numSN)
            M = - self.alpha * x1vals - self.beta * cvals
            Mabs = self.centralMabs + M + self.randomState.normal(loc=0., scale=self.Mdisp,
                                               size=self.numSN)
            x0 = np.zeros(self.numSN)
            mB = np.zeros(self.numSN)
            model = sncosmo.Model(source='SALT2')
            for i, z in enumerate(self.zSamples):
                model.set(z=z, x1=x1vals[i], c=cvals[i])
                model.set_source_peakabsmag(Mabs[i], 'bessellB', 'ab',
                                            cosmo=self.cosmo)
                x0[i] = model.get('x0')
                mB[i] = model.source.peakmag('bessellB', 'ab')
            df = pd.DataFrame(dict(x0=x0, mB=mB, x1=x1vals, c=cvals, M=M, Mabs=Mabs,
                                   t0=T0Vals, z=self.zSamples, snid=self.snids))
            self._paramSamples = df
        return self._paramSamples


class GMM_SALT2Params(SimpleSALTDist):
    """
    """
    def __init__(self, numSN, zSamples, snids=None, alpha=0.11, beta=3.14,
                 Mdisp=0.15, rng=None, cosmo=Planck15, mjdmin=0., surveyDuration=10.):
        super(self.__class__, self).__init__(numSN, zSamples, snids=snids, alpha=alpha,
                              beta=beta, rng=rng, cosmo=cosmo, mjdmin=mjdmin,
                              surveyDuration=SurveyDuration, Mdisp=Mdisp)

    @property
    def paramSamples(self):
        if self._paramSamples is not None:
            return self._paramSamples
        timescale = self.mjdmax - self.mjdmin
        T0Vals = self.randomState.uniform(size=self.numSN) * timescale \
            + self.mjdmin
        mB, x1, c, m = SALT2_MMDist(self.numSN)
        x0 = np.zeros(len(mB))
        mB += self.randomState.normal(loc=0., scale=self.Mdisp,
                                      size=self.numSN)
        model = sncosmo.Model(source='SALT2')
        for i, z in enumerate(self.zSamples):
            model.set(z=z, x1=x1[i], c=c[i])
            model.source.set_peakmag(mB[i], 'bessellB', 'ab')
            x0[i] = model.get('x0')
        df = pd.DataFrame(dict(x0=x0, x1=x1, c=c, mB=mB, z=self.zSamples))
        self._paramSamples = df
        return self._paramSamples

class CoordSamples(PositionSamples, HealpixTiles):
    def __init__(self, nside, hpOpSim, rng):
        self.nside = nside
        self.nside = nside
        super(self.__class__, self).__init__(nside=nside, healpixelizedOpSim=hpOpSim)
        self._rng = rng
    @property
    def randomState(self):
        if self._rng is None:
            raise ValueError('self._rng should not be None')
        return self._rng
    def _angularSamples(self, phi_c, theta_c, radius, numSamples, tileID):
        phi, theta = super(self.__class__, self).samplePatchOnSphere(phi=phi_c,
								     theta=theta_c, delta=radius, 
                                                                     size=numSamples, rng=self.randomState)
        tileIds = hp.ang2pix(nside=self.nside, theta=np.radians(theta),
			     phi=np.radians(phi), nest=True)
        inTile = tileIds == tileID
        return phi[inTile], theta[inTile]
        
    def positions(self, tileID, numSamples):
        res_phi = np.zeros(numSamples)
        res_theta = np.zeros(numSamples)
        ang = hp.pix2ang(nside=self.nside, ipix=tileID, nest=True)
        radius = np.degrees(np.sqrt(hp.nside2pixarea(self.nside) / np.pi))
        phi_c, theta_c = np.degrees(ang[::-1])
        num_already = 0
        while numSamples > 0:
            phi, theta = self._angularSamples(phi_c, theta_c, radius=2*radius,
					      numSamples=numSamples,
					      tileID=tileID)
            num_obtained = len(phi)
            res_phi[num_already:num_obtained + num_already] = phi
            res_theta[num_already:num_obtained + num_already] = theta
            num_already += num_obtained
            numSamples -= num_obtained
        return res_phi, res_theta


class PowerLawRates(RateDistributions):
    """
    This class is a concrete implementation of `RateDistributions` with the
    following properties:
    - The SN rate : The SN rate is a single power law with numerical
        coefficients (alpha, beta)  passed into the instantiation. The rate is
        the number of SN at redshift z per comoving volume per unit observer
        time over the entire sky expressed in units of numbers/Mpc^3/year 
    - A binning in redshift is used to perform the calculation of numbers of SN.
        This is assumed
    - The expected number of SN in each of these redshift bins is computed using
        the rate above, and a cosmology to compute the comoving volume for the
        redshift bin
    - The numbers of SN are determined by a Poisson Distribution about the
        expected number in each redshift bin,  determined with a random state
        passed in as an argument. This number must be integral.
    - It is assumed that the change of rates and volume within a redshift bin
        is negligible enough that samples to the true distribution may be drawn
        by obtaining number of SN samples of z from a uniform distribution
        within the z bin.

    """

    def __init__(self,
                 rng,
                 alpha=2.6e-5, beta=1.5,
                 zbinEdges=None,
                 zlower=1.0e-8,
                 zhigher=1.4,
                 numBins=20,
                 surveyDuration=10., # Unit of  years
                 fieldArea=None, # Unit of degree square
                 skyFraction=None,
                 cosmo=Planck15):
        """
        Parameters
        ----------
        rng : instance of `np.random.RandomState`
        cosmo : Instance of `astropy.cosmology` class, optional, defaults to Planck15
            data structure specifying the cosmological parameters 
        alpha : float, optional, defaults to 2.6e-5
            constant amplitude in SN rate powerlaw
        beta : float, optional, defaults to 1.5
            constant exponent in SN rate powerlaw
        surveyDuration : float, units of years, defaults to 10.
            duration of the survey in units of years
        fieldArea : float, units of degree sq, defaults to None
            area of the field over which supernovae are being simulated in
            units of square degrees.
        skyFraction : float, optional, defaults to None
            Alternative way of specifying fieldArea by supplying the unitless
            ratio between sky area where the supernovae are being simulated and
            accepting the default `None` for `fieldArea`
        """
        self.alpha = alpha
        self.beta = beta
        self.cosmo = cosmo
        self.zlower = zlower
        self.zhigher = zhigher
        self.numBins = numBins
        self._zbinEdges = zbinEdges
        self._rng = rng
        self.DeltaT = surveyDuration
        self.fieldArea = fieldArea
        self._skyFraction = skyFraction
        # not input
        self._numSN = None
        self._zSamples = None

    @property
    def skyFraction(self):
        if self._skyFraction is None:
            if self.fieldArea is None:
                raise ValueError('both fieldArea and skyFraction cannot be given')
            self._skyFraction = self.fieldArea * np.radians(1.)**2.0  / 4.0 / np.pi
        return self._skyFraction

    @property
    def randomState(self):
        if self._rng is None:
            raise NotImplemented('rng must be provided')
        return self._rng

    @property
    def zbinEdges(self):
        if self._zbinEdges is None:
            if any(x is None for x in (self.zlower, self.zhigher, self.numBins)):
                raise ValueError('Both zbinEdges, and'
                                 '(zlower, zhigher, numBins) cannot be None')
            if self.zlower >= self.zhigher:
                raise ValueError('zlower must be less than zhigher')
            self._zbinEdges = np.linspace(self.zlower, self.zhigher, self.numBins + 1)
        return self._zbinEdges




    def snRate(self, z):
        """
        The rate of SN at a redshift z in units of number of SN/ comoving
        volume in Mpc^3/yr in earth years according to the commonly used
        power-law expression 

        .. math:: rate(z) = \alpha (h/0.7)^3 (1.0 + z)^\beta
        
        Parameters
        ----------
        z : array-like, mandatory 
            redshifts at which the rate is evaluated

        Examples
        --------
        """
        res = self.alpha * (1.0 + z)**self.beta 
        res *= ((self.cosmo.h / 0.7) **3.)  
        return res

    def zSampleSize(self): 
        #, zbinEdges=self.zbinEdges, DeltaT=self.DeltaT,
        #            skyFraction=self.skyFraction,
        #            zlower=None, zhigher=None, numBins=None):
        """
        Parameters
        ----------
        zbinEdges : `nunpy.ndarray` of edges of zbins, defaults to None
            Should be of the form np.array([z0, z1, z2]) which will have
            zbins (z0, z1) and (z1, z2)
        skyFraction : np.float, optional, 
        fieldArea : optional, units of degrees squared
            area of sky considered.
        zlower : float, optional, defaults to None
            lower edge of z range
        zhigher : float, optional, defaults to None
            higher edge of z range
        numBins : int, optional, defaults to None
           if not None, overrides zbinEdges
        """
        DeltaT = self.DeltaT
        skyFraction = self.skyFraction
        zbinEdges = self.zbinEdges
        z_mids = 0.5 * (zbinEdges[1:] + zbinEdges[:-1])
        snpervolume = self.snRate(z_mids) 

        # Comoving volume of the univere in between zlower and zhigher
        vols = self.cosmo.comoving_volume(zbinEdges)

        # Comoving volume in each bin
        vol = vols[1:] - vols[:-1]
        vol *= skyFraction
        
        numSN = vol * snpervolume * DeltaT / (1.0 + z_mids)
        return numSN.value

    def numSN(self):
        """
        Return the number of expected supernovae in time DeltaT, with a rate snrate
        in a redshift range zlower, zhigher divided into numBins equal redshift
        bins. The variation of the rate within a bin is ignored.
    
        Parameters
        ----------
        zlower : mandatory, float
            lower limit on redshift range
        zhigher : mandatory, float
            upper limit on redshift range
        numBins : mandatory, integer
            number of bins
        cosmo : `astropy.cosmology` instance, mandatory
            cosmological parameters
        fieldArea : mandatory, units of radian square
            sky area considered
        """
        if self._numSN is None:
            lam = self.zSampleSize()
            self._numSN = self.randomState.poisson(lam=lam)
        return self._numSN


    @property
    def zSamples(self):
        # Calculate only once
        if self._zSamples is None:
            numSN = self.numSN()
            zbinEdges =  self.zbinEdges
            x = zbinEdges[:-1]
            y = zbinEdges[1:]
            arr = (self.randomState.uniform(low=xx, high=yy, size=zz).tolist()
                                           for (xx, yy, zz) in zip(x,y, numSN))
            self._zSamples = np.asarray(list(__x for __lst in arr
                                             for __x in __lst))
        return self._zSamples


class TwinklesRates(PowerLawRates):
    def __init__(self, galsdf, rng, alpha=2.6e-3, beta=1.5, zbinEdges=None,
                 zlower=0.0000001, zhigher=1.2, numBins=24, agnGalids=None,
                 surveyDuration=10., fieldArea=None, skyFraction=None,
                 cosmo=Planck15):
        PowerLawRates.__init__(self, rng=rng, alpha=alpha, beta=beta,
                               zbinEdges=zbinEdges, zlower=zlower,
                               zhigher=zhigher, numBins=numBins,
                               fieldArea=fieldArea, cosmo=cosmo)
        self._galsdf = galsdf
        if agnGalids is None:
            agnGalids = []
        self.agnGalTileIds = tuple(agnGalids)
        self.binWidth = np.diff(self.zbinEdges)[0]
        #self.galsdf =None
        self.binnedGals = None
        self.numGals = None
        self.gdf = None
        self.rng = rng
        self._selectedGals = None
        
    
    @property
    def galsdf(self):
        if self.gdf is not None:
            return self.gdf
        zhigher = self.zhigher
        self.addRedshiftBins()
        
        vetoedGaltileIds = tuple(self.agnGalTileIds)
        sql_query = 'redshift <= @zhigher and galtileid not in @vetoedGaltileIds'
        galsdf = self._galsdf.query(sql_query)
        self.binnedGals = galsdf.groupby('redshiftBin')
        self.numGals = self.binnedGals.redshift.count()
        galsdf['probHist'] = galsdf.redshiftBin.apply(self.probHost)
        galsdf['hostAssignmentRandom'] = self.rng.uniform(size=len(galsdf))
        self.gdf = galsdf
        return galsdf
    
    def addRedshiftBins(self):
        self._galsdf['redshiftBin'] = (self._galsdf.redshift - self.zlower) // self.binWidth
        self._galsdf.redshiftBin = self._galsdf.redshiftBin.astype(np.int)
        
    @property    
    def selectedGals(self):
        if self._selectedGals is None:
            df = self.galsdf.query('hostAssignmentRandom < probHist')
            df.galtileid = df.galtileid.astype(int)
            df['snid'] = df.id
        else:
            df = self._selectedGals
        return df
    def probHost(self, binind):
        return np.float(self.numSN()[binind]) / self.numGals[binind]
    
    @property
    def zSamples(self):
        return self.selectedGals.redshift.values

class CatSimPositionSampling(object):
    def __init__(self, rng, galdf, snAngularUnits='degrees'):
        self.galdf = galdf.copy()
        self.rng = rng
        self.ss = SersicSamples(rng=self.rng)
        self.radianOverArcSec = np.pi / 180.0 / 3600.
        self.snAngularUnits=snAngularUnits

    
    def f1(self, x):
        return self.ss.sampleRadius(x)[0]
    def f4(self, x):
        return self.ss.sampleRadius(x, sersicIndex=4)[0]
    
    def SampleDiskAngles(self, x):
        return self.ss.sampleAngles(x.a_d, x.b_d)[0]
    def SampleBulgeAngles(self, x):
        return self.ss.sampleAngles(x.a_b, x.b_b)[0]
    @staticmethod
    def theta(df, angle='diskAngle', PositionAngle='pa_disk'):
        return np.radians(df[angle] - df[PositionAngle] + 90.)

    @staticmethod
    def snInDisk(x, value=" None"):
        if x.sedFilenameDisk == value:
            return 0
        elif x.sedFilenameBulge == value:
            return 1
        else:
            return np.random.choice([0, 1], p=[0.5, 0.5])
    def addPostions(self):
        self.galdf['isinDisk'] = self.galdf.apply(self.snInDisk, axis=1)
        self.galdf['bulgeRadialPos'] = self.galdf.BulgeHalfLightRadius.apply(self.f4)
        self.galdf['diskRadialPos'] = self.galdf.DiskHalfLightRadius.apply(self.f1)
        self.galdf['bulgeAngle'] = self.galdf.apply(self.SampleBulgeAngles, axis=1)
        self.galdf['diskAngle'] = self.galdf.apply(self.SampleDiskAngles, axis=1)
        self.galdf['DeltaRaDisk'] = self.galdf.diskRadialPos * np.cos(self.theta(self.galdf)) * self.galdf.isinDisk
        self.galdf['DeltaRaBulge'] = self.galdf.bulgeRadialPos * self.theta(self.galdf, angle='bulgeAngle', 
                                                                            PositionAngle='pa_bulge')\
                                     * (1 - self.galdf.isinDisk)
        self.galdf['DeltaDecDisk'] = self.galdf.diskRadialPos * np.sin(self.theta(self.galdf)) * self.galdf.isinDisk
        self.galdf['DeltaDecBulge'] = self.galdf.bulgeRadialPos * np.sin(self.theta(self.galdf, angle='bulgeAngle',
                                                                                 PositionAngle='pa_bulge')) \
                                                                                 * (1 - self.galdf.isinDisk)
        self.galdf['snra'] = self.radianOverArcSec *\
            self.galdf[['DeltaRaDisk', 'DeltaRaBulge']].apply(np.nansum, axis=1)\
            + self.galdf.raJ2000

        self.galdf['sndec'] = self.radianOverArcSec *\
            self.galdf[['DeltaDecDisk', 'DeltaDecBulge']].apply(np.nansum, axis=1)\
            + self.galdf.decJ2000

        if self.snAngularUnits == 'degrees':
            self.galdf[['snra', 'sndec']] = \
                    self.galdf[['snra', 'sndec']].apply(np.degrees)
        elif self.snAngularUnits =='radians':
            pass
        else :
            raise NotImplementedError('conversion to snAngularUnits {} not implemented', self.snAngularUnits)

class TwinklesSim(TwinklesRates):

    def __init__(self, catsimgaldf, rng, fieldArea, cosmo, agnGalids=None, numBins=24,
                 rate_alpha=0.0026, rate_beta=1.5, zlower=1.0e-7, zhigher=1.2,
                 zbinEdges=None, tripp_alpha=0.11, tripp_beta=3.14, mjdmin=0.):
        super(TwinklesSim, self).__init__(catsimgaldf, rng=rng,
                                          cosmo=cosmo, fieldArea=fieldArea,
                                          agnGalids=agnGalids,
                                          alpha=rate_alpha, beta=rate_beta,
                                          zlower=zlower, zhigher=1.2,
                                          numBins=numBins, zbinEdges=None,
                                          skyFraction=None)
        self.cosmo = cosmo
        self.beta_rate = deepcopy(self.beta)
        self.alpha_rate = deepcopy(self.alpha)
        self.numSN = len(self.zSamples)
        self.tripp_alpha = tripp_alpha
        self.tripp_beta = tripp_beta
        self.catsimpos = CatSimPositionSampling(rng=self.rng, 
                                                galdf=self.selectedGals)
        self.catsimpos.addPostions()
        self.mjdmin = mjdmin
        self.salt2params = SimpleSALTDist(numSN=self.numSN, 
                                          zSamples=self.zSamples, 
                                          rng=self.rng,
                                          cosmo=self.cosmo,
                                          snids=self.catsimpos.galdf.snid,
                                          alpha=self.tripp_alpha,
                                          beta=self.tripp_beta,
                                          surveyDuration=10.,
                                          mjdmin=self.mjdmin)
        self._snparamdf = None

    @property
    def snparamdf(self):
        """
        dataFrame with information about the object
        """
        if self._snparamdf is None:
            self.salt2params.paramSamples.set_index('snid', inplace=True)
            self.catsimpos.galdf.set_index('snid', inplace=True)
            self._snparamdf = self.salt2params.paramSamples.join(self.catsimpos.galdf)
        return self._snparamdf
