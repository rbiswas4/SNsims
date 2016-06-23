import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catUtils.mixins import CosmologyMixin
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.photUtils import BandpassDict



class SNTile(object):
    def __init__(self, tileID):
        self.tileID = tileID

    @property
    def obsMetaDataSpatial(self):
        pass

lsst_bp = BandpassDict.loadTotalBandpassesFromFiles()

degConv = np.array([1., 1./60., 1./3600.])
raConv = degConv / 24.0 * 360.
centralRA = np.dot(np.array([3., 32., 30]), raConv) #03h 32m 30s
centralDec = np.dot(np.array([-28, 6., 0.]), degConv)
patchRadius = 0.4 * np.sqrt(2) #np.dot(np.array([0.0, 10.0, 0.]), degConv)



area = np.pi * (0.4 * np.sqrt(2.))**2
factorLarger = area / 0.16 / 0.16; print(factorLarger)
NumHighSNdesired = factorLarger * 100; print (NumHighSNdesired)


# In[13]:

print(centralRA, centralDec, patchRadius)


# In[14]:

TwinklesObsMetaData = ObservationMetaData(boundType='circle',pointingRA=centralRA,pointingDec=centralDec,
                                          boundLength=patchRadius, mjd=49540.0)


# In[15]:

TwinklesObsMetaDataSmall = ObservationMetaData(boundType='box',pointingRA=centralRA,pointingDec=centralDec,
                                          boundLength=0.167, mjd=49540.0)


# In[16]:

#The following is to get the object ids in the registry
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catUtils.baseCatalogModels.GalaxyModels import GalaxyTileObj


# In[17]:

galaxyTiled  = GalaxyTileObj()


# In[18]:

class galCopy(InstanceCatalog):
    column_outputs = ['galtileid', 'raJ2000', 'decJ2000', 'redshift', 'a_d', 'b_d', 'pa_disk']
    override_formats = {'raJ2000': '%8e', 'decJ2000': '%8e', 'a_d': '%8e', 'b_d': '%8e', 'pa_disk': '%8e'}


# In[19]:

TwinklesSmall = galCopy(galaxyTiled, obs_metadata=TwinklesObsMetaDataSmall)


# In[20]:

TwinklesSmall.write_catalog('twinklesSmall.dat')


# In[21]:

TwinklesGalaxies = galCopy(galaxyTiled, obs_metadata=TwinklesObsMetaData)


# In[22]:

TwinkSmallGalsdf = pd.read_csv('TwinklesSmall.dat', sep=',\s+', engine='python')
TwinkSmallGalsdf.rename(columns={'#galtileid':'galtileid'}, inplace=True)


# In[23]:

TwinkSmallGalsdf.head()


# In[24]:

TwinkSmallGalsdf['zbin']= TwinkSmallGalsdf['redshift'] // 0.1


# In[25]:

TwinklesGalaxies.write_catalog('TwinklesGalaxies.dat')


# In[26]:

TwinkGalsdf = pd.read_csv('TwinklesGalaxies.dat', sep=',\s+', engine='python', index_col=0)


# In[27]:

len(TwinkGalsdf)


# In[28]:

fig, ax = plt.subplots()
ax.plot(np.degrees(TwinkGalsdf.raJ2000.values), np.degrees(TwinkGalsdf.decJ2000.values), 'o')
ax.set_aspect('equal')
ax.set_ylabel('dec')
ax.set_xlabel('ra')
ax.set_title('Sprinkled Region')
fig.savefig('Twinkles_Area.png')


# In[29]:

TwinkGalsdf['zbin']= TwinkGalsdf['redshift'] // 0.1
zmids = np.arange(0.05, (TwinkGalsdf.zbin.max()+1.)* 0.1, 0.1)
print(zmids)
zbinnedGals = TwinkGalsdf.groupby('zbin')
binnedTwinks = pd.DataFrame({'zmids': zmids})
binnedTwinks['counts'] = zbinnedGals['redshift'].count()


# In[30]:

fig, ax = plt.subplots()
ax.errorbar(binnedTwinks.zmids, binnedTwinks.counts, np.sqrt(binnedTwinks.counts), fmt='o' )
ax.set_xlim(0., 1.4)
ax.set_xlabel('redshift')
ax.set_ylabel('Number of galaxies in bins')


# ## SN light Curves

# In[31]:

import sncosmo
import gedankenLSST


# In[32]:

lsstchar = gedankenLSST.LSSTReq
lsstchar['meanNumVisits'] = pd.Series(np.repeat(3650.,6), index=['u','g','r','i','z','y'])
lsstchar['meanNumVisits']
sn = gedankenLSST.GSN_Obs(mjd_center=49530., lsstrequirements=lsstchar)
sndf = sn.summary
sndf[sndf['filter'] == 'u'].hist('night',bins=80)


# In[33]:

lsstchar['medianSVD']


# In[34]:

s = gedankenLSST.SNObs(summarydf=sndf, t0=49530, lsst_bp=lsst_bp, ra=centralRA, dec=centralDec)


# In[35]:

a = []
for z in binnedTwinks.zmids.values[:12]:
    s.snState = {'z': z}
    lc = s.lightcurve
    totalEpochs = len(lc)
    highSNRlc = lc.query('SNR > 5')
    highSNREpochs = len(highSNRlc)
    highSNREpochs_u = len(highSNRlc.query("filter == 'u'"))
    highSNREpochs_g = len(highSNRlc.query("filter == 'g'"))
    highSNREpochs_r = len(highSNRlc.query("filter == 'r'"))
    highSNREpochs_i = len(highSNRlc.query("filter == 'i'"))
    highSNREpochs_z = len(highSNRlc.query("filter == 'z'"))
    highSNREpochs_y = len(highSNRlc.query("filter == 'y'"))
    
    a.append([z, highSNREpochs, highSNREpochs_u, highSNREpochs_g, highSNREpochs_r, highSNREpochs_i, highSNREpochs_z,
              highSNREpochs_y, totalEpochs, -2.5 * np.log10(s.SN.get('x0'))])


# In[36]:

FlatzSummary  = pd.DataFrame(a, columns=['redshift', 'highSNREpochs', 'u', 'g', 'r', 'i', 'z', 'y', 'totalEpochs', 'mB'])


# In[37]:

FlatzSummary['frac'] = FlatzSummary.highSNREpochs / FlatzSummary.totalEpochs


# In[38]:

numSNperZBinDesired = NumHighSNdesired /12.
FlatzSummary['NumSNperzBin'] = numSNperZBinDesired * 3650. / 80. / FlatzSummary['frac']


# In[39]:

_nsn = FlatzSummary.NumSNperzBin.replace([-np.inf, np.inf], 0.) * FlatzSummary.frac
print(_nsn.sum() * 80. / 3650.) 


# In[40]:

# Increase the numbers since some of the bins are empty


# In[41]:

FlatzSummary['NumSNperzBin'] =  FlatzSummary['NumSNperzBin'] * (12./9.)


# In[42]:

_nsn = FlatzSummary.NumSNperzBin.replace([-np.inf, np.inf], 0.) * FlatzSummary.frac
print(_nsn.sum() * 80. / 3650.) 


# In[43]:

FlatzSummary['numGalsperzBin'] = binnedTwinks['counts'].head(12)


# In[44]:

FlatzSummary['numSNperGal'] = FlatzSummary['NumSNperzBin'] / FlatzSummary['numGalsperzBin'] 


# In[45]:

FlatzSummary


# In[46]:

plt.plot(FlatzSummary.redshift, FlatzSummary['NumSNperzBin'].replace(np.inf,0), 'o')


# #  SN Table

# In[47]:

model = sncosmo.Model(source='salt2')


# In[48]:

from astropy.cosmology import FlatLambdaCDM


# ### Simulation Parameters

# In[49]:

# Astropy cosmology object for CatSim Cosmology
CatSimCosmo = FlatLambdaCDM(Om0=0.25, H0=73.)

alphaTwinkles = 0.11
betaTwinkles = -3.14
cdistTwinkles = [0., 0.1]
x1distTwinkles = [0, 1.]
MTwinkles = [-19.3, 0.15]


# In[50]:

zbinnedGals = TwinkGalsdf.groupby('zbin')


# In[51]:

def assignIds(snwithHosts, maxval=100000000 * 10000 * 100):
    snwithHosts['offset'] = 0
    sngroups = snwithHosts.groupby('galtileid')
    for host in (sngroups.count() > 0).index.values:
        sn = sngroups.get_group(host)
        idx  = sn.index
        snwithHosts.loc[idx, 'offset'] = np.arange(len(sn))
    return None


# In[52]:

def assignSNHosts(galdf, numSN, seed):
    if seed is not None:
        np.random.seed(seed)
    sngalids = np.random.choice(galdf.index.values, numSN, replace=True)
    zvals = galdf.ix[sngalids,'redshift']
    df = pd.DataFrame({'galtileid': sngalids, 
                      'redshift' : zvals.values})
    return df


# In[53]:

# Slow step: Takes about 20 mins
def assignSN(zbinnedGals, SNzSummary, binList=[0, 1], maxval=100000000 * 10000 * 100, seed=42):
    
    dfs = []
    for idx in binList:
        galdf = zbinnedGals.get_group(idx)
        numSN = SNzSummary.NumSNperzBin[idx]
        if idx == 0 :
            snWithHosts = assignSNHosts(galdf, numSN, seed)
        else:
            snWithHosts = assignSNHosts(galdf, numSN, seed=None)
        assignIds(snWithHosts, maxval=maxval)
        dfs.append(snWithHosts)
    snvals = pd.concat(dfs)
    snvals['snid'] = snvals['galtileid'] *100 + snvals['offset']
    return snvals


# In[54]:

# Slow step ~ 20 mins
snvals = assignSN(zbinnedGals, FlatzSummary, binList=[0, 1, 2, 3, 4, 5, 6, 7, 8])


# In[55]:

import time 
print time.time()


# In[56]:

snvals.set_index(snvals['snid'], drop=True, verify_integrity=True, inplace=True)


# In[57]:

def assigSNParams(sntable, seed=42, cosmo=None, T0Min=0., T0Max=3650., 
                  MabsScatter= [-19.3, 0.15], cScatter=[0., 0.1], x1Scatter=[0., 1.], alpha=0.11, beta=-3.14 ):
    if seed is not None:
        np.random.seed(seed)
        
    model = sncosmo.Model(source='salt2')
    if cosmo is None:
        cosmo = FlatLambdaCDM(Om0=0.25, H0=73.)
        
    numSN = len(sntable)
    zvals = sntable.redshift.values
    cvals = np.random.normal(cScatter[0], cScatter[1], size=numSN)
    x1vals = np.random.normal(x1Scatter[0], x1Scatter[1], size=numSN)
    M  = np.random.normal(MabsScatter[0], MabsScatter[1], size=numSN)
    
    M += -alpha * x1vals - beta * cvals 
    t0 = np.random.uniform(T0Min, T0Max, size=numSN)
    x0 = np.zeros(numSN)
    mB = np.zeros(numSN)
    # Slow Step
    for i, Mabs in enumerate(M):
        model.set(z=zvals[i], c=cvals[i], x1=x1vals[i])
        model.set_source_peakabsmag(Mabs, 'bessellB', 'ab', cosmo=cosmo)
        x0[i] = model.get('x0')
        mB[i] = model.source.peakmag('bessellB', 'ab')
    sntable['t0'] = t0
    sntable['c'] = cvals
    sntable['x1'] = x1vals
    sntable['x0'] = x0
    sntable['mB'] = mB
    sntable['M'] = M
    
    print (alpha, beta, cScatter, x1Scatter, MabsScatter)
    


# In[58]:

starttime = time.time()
assigSNParams(sntable=snvals, cosmo=CatSimCosmo, alpha=alphaTwinkles, beta=betaTwinkles, MabsScatter=MTwinkles, 
              seed=24)
endtime = time.time()
print("Time taken", endtime - starttime)


# In[ ]:

def assignPositions(sntable, Galsdf, seed=42):

    radiansPerArcSec = (np.pi / 180.)* (1./60.)**2
    if seed is not None:
        np.random.seed(seed)
    
    r1 = np.random.normal(0., 1., sntable.snid.size)
    r2 = np.random.normal(0., 1., sntable.snid.size)
    
    sntable['raJ2000'] = Galsdf.ix[sntable.galtileid, 'raJ2000'].values
    sntable['decJ2000'] = Galsdf.ix[sntable.galtileid, 'decJ2000'].values
    sntable['a_d'] = Galsdf.ix[sntable.galtileid, 'a_d'].values * radiansPerArcSec
    sntable['b_d'] = Galsdf.ix[sntable.galtileid, 'b_d'].values * radiansPerArcSec

    # convert from degrees to radians
    sntable['theta'] = np.radians(Galsdf.ix[sntable.galtileid, 'pa_disk'].values)
    
    sntable['sndec'] = np.cos(-sntable['theta']) * sntable['a_d']* r1 + np.sin(-sntable['theta'])*sntable['b_d'] * r2
    sntable['snra'] =  - np.sin(-sntable['theta']) * sntable['a_d']*r1 + np.cos(-sntable['theta'])* sntable['b_d'] * r2
    sntable['snra'] += Galsdf.ix[sntable.galtileid, 'raJ2000'].values
    sntable['sndec'] += Galsdf.ix[sntable.galtileid, 'decJ2000'].values
    
    sntable['sndec'] = np.degrees(sntable['sndec'])
    sntable['snra'] = np.degrees(sntable['snra'])
    sntable['raJ2000'] = np.degrees(sntable['raJ2000'])
    sntable['decJ2000'] = np.degrees(sntable['decJ2000'])



# In[ ]:

assignPositions(snvals, TwinkGalsdf, seed=4 )


# In[ ]:

snvals.raJ2000.hist(histtype='step', alpha=1, lw=4, color='k')
snvals.snra.hist(histtype='step', lw=1, color='r')


# In[ ]:

(snvals.snra - snvals.raJ2000).hist(histtype='step', lw=1, color='r',bins=50, **{'normed':True, 'log':True})
(snvals.sndec - snvals.decJ2000).hist(histtype='step', lw=1, color='b',bins=50, **{'normed':True, 'log':True})


# In[ ]:

plt.hist(snvals.sndec - snvals.decJ2000,histtype='step', lw=1, color='b',bins=50,
        log=True)


# In[ ]:

snvals.head()


# In[ ]:

snvals.columns


# In[ ]:

(snvals.raJ2000.iloc[0] - snvals.snra.iloc[0]) * 180. / np.pi * 3600.


# In[ ]:

snvals.head()


# In[ ]:

(3600* (snvals.snra - snvals.raJ2000).apply(np.degrees)).hist(bins=20)


# In[ ]:

snvals.head()


# In[ ]:

snvals.to_csv('TwinklesSN_new.csv')


# In[ ]:

#snvals.to_csv('TwinklesSN.csv')


# In[ ]:

len(np.unique(snvals.index.values)) == len(snvals.index.values)

