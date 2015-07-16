import numpy as np
import haloutils
#import readsnapshots.readsnapHDF5_greg as rsg
from caterpillaranalysis import *
from scipy import interpolate
#from scipy.integrate import quad
import matplotlib.pyplot as plt
import grifflib as glib
#from caterpillaranalysis import ProfilePlugin


#hpath = '/bigbang/data/AnnaGroup/caterpillar/halos/H1599988/H1599988_EX_Z127_P7_LN7_LX14_O4_NV4'

def distance(posA, posB,boxsize=100*1000.):
    dist = abs(posA-posB)
    tmp = dist > boxsize/2.0
    dist[tmp] = boxsize-dist[tmp]
    if dist.shape == (3,):
        return np.sqrt(np.sum(dist**2))
    else:
        return np.sqrt(np.sum(dist**2,axis=1))

class SubstructurePlugin(PluginBase):
    def __init__(self):
        super(SubstructurePlugin,self).__init__()
        self.filename='substructure_data.dat'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = False
        self.xlabel='' ; self.ylabel='' 

    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")

        # what to quantify substructure
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID = haloutils.load_zoomid(hpath)
        hosthalo = cat.ix[hostID]
        hostpos = np.array(hosthalo[['posX','posY','posZ']])
        subs = cat.get_subhalos_within_halo(hostID) # just 1 level deep
        subs_all = cat.get_all_subhalos_within_halo(hostID)
        
        # num halos above some vmax thresh currently
        maskv = subs_all['vmax']> 25
        nsubs_vmax = np.sum(maskv)

        # num halos that ever entered above some vmax thresh
        import MTanalysis2 as MT
        AE = MT.AllExtantData()
        AD = MT.AllDestroyedData()
        dataE = AE.read(hpath)
        dataD = AD.read(hpath)
        maskE = dataE['infall_vmax'] > 25
        maskD = dataD['infall_vmax'] > 25
        nsubs_vmax_ever = np.sum(maskE)+np.sum(maskD)
        
        # SHMF normalization (divided by halo mass)
        # need to extract from something else
        shmf = SHMFPlugin()
        x,y,sx,sy = shmf.read(hpath)
        # x in Msun, y in dN/dM
        
        # total subhalo mass?
        subhalomass = np.sum(subs['mgrav']/cat.h0)

        # subhalo mass fraction
        submass_frac = subhalomass/(hosthalo['mgrav']/cat.h0)
        
        # num subhalos within 50% of rvir vs outside of it
        subposns = np.array(subs_all[['posX','posY','posZ']])
        dist = distance(hostpos,subposns)*cat.scale/cat.h0*1000
        subs_inner_frac = np.sum(dist<.5*hosthalo['rvir']/cat.h0)/float(len(subs_all))

        nsubs = len(subs_all)

        g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,'wb')
        np.array([nsubs,nsubs_vmax, nsubs_vmax_ever,subhalomass,submass_frac,subs_inner_frac]).tofile(g)
        g.close()        

    def _read(self,hpath):
        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename)
        return data
        

substruct = SubstructurePlugin()
#nsubs_vmax, nsubs_vmax_ever,subhalomass,submass_frac,subs_inner_frac = substruct.read(hpath)
lx = 14
figdir = '/bigbang/data/AnnaGroup/GregFigs/lx'+str(lx)+'/'
#halo_paths = haloutils.find_halo_paths(levellist=[lx],require_mergertree=True,require_subfind=False,verbose=False) 

halo_paths = ['/bigbang/data/AnnaGroup/caterpillar/halos/H95289/H95289_BB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1195448/H1195448_EC_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1354437/H1354437_EA_Z127_P7_LN7_LX14_O4_NV5',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H447649/H447649_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H5320/H5320_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H581141/H581141_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1631506/H1631506_EA_Z127_P7_LN7_LX14_O4_NV5',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1130025/H1130025_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1725139/H1725139_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H581180/H581180_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1422331/H1422331_EX_Z127_P7_LN7_LX14_O4_NV5',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1232164/H1232164_EX_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1387186/H1387186_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1725272/H1725272_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1725372/H1725372_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H264569/H264569_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1292085/H1292085_EX_Z127_P7_LN7_LX14_O4_NV5',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H94687/H94687_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1268839/H1268839_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H196589/H196589_EX_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1599988/H1599988_EX_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1079897/H1079897_EX_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H94638/H94638_EX_Z127_P7_LN7_LX14_O4_NV5',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H796175/H796175_EB_Z127_P7_LN7_LX14_O4_NV4']


for hpath in halo_paths:
    nsubs, nsubs_vmax, nsubs_vmax_ever,subhalomass,submass_frac,subs_inner_frac = substruct.read(hpath)
    print subs_inner_frac, 'inner fraction'
    print 'done with', hpath
