import numpy as np
import pylab as plt
from caterpillaranalysis import *
import readsnapshots.readsnapHDF5_greg as rsg
import haloutils
import time
from scipy import interpolate
from scipy.integrate import quad
import os, subprocess
import pandas


# compute distance from posA to posB.       
 # posA can be an array. boxsize must be in same units as positions. 
def distance(posA, posB,boxsize=100.):
    dist = abs(posA-posB)
    tmp = dist > boxsize/2.0
    dist[tmp] = boxsize-dist[tmp]
    if dist.shape == (3,):
        return np.sqrt(np.sum(dist**2))
    else:
        return np.sqrt(np.sum(dist**2,axis=1))

# functions for pericenter
def get_pericenter(sub_mb, iSnap, host_mb,cat):
    """
    rsid: value of rockstar id at z=0                    
    """
    orbit_pos = trackPosition(sub_mb,host_mb,iSnap)*1000
    scale_list = sub_mb['scale'][0:len(orbit_pos)][::-1]
    orbital_dist = distance(orbit_pos, np.array([0,0,0]),boxsize=cat.boxsize*1000)*scale_list/cat.h0
    perisnap, min_peridist, mean_peridist = getPericenterInfo(orbital_dist,sub_mb)
    #print perisnap, min_peridist, mean_peridist, 'snap, peri, mean peri'                 
    return min_peridist, mean_peridist

def trackPosition(sub_mb,host_mb,iSnap,verbose=False):
    """                                                                                                              
    @param sub: MTCatalogueTree of the main branch of the subhalo 
    @param host: MTCatalogueTree of the main branch of the host halo                    
    @return: position of subhalo relative to host over time. 
    starting with position at infall, ending with z=0 
    """
    # make both masks start at 0 end at min(len(sub_mb), len(host_mb))
    mask_sub = sub_mb['snap']>=iSnap
    mask_host = host_mb['snap']>=iSnap
    ll = min(np.sum(mask_sub), np.sum(mask_host))
    assert(host_mb[0]['snap']==sub_mb[0]['snap']), "ERROR: Mismatch in alignment of host_mb and sub_mb"

    sub_pos = np.array(map(list,sub_mb[mask_sub][0:ll][['posX','posY','posZ']]))[::-1]
    host_pos = np.array(map(list,host_mb[mask_host][0:ll][['posX','posY','posZ']]))[::-1]
    return sub_pos-host_pos

    """
    if len(sub_mb) < len(host_mb):
        mask = sub_mb['snap']>=iSnap
    else:
        mask = host_mb['snap']>=iSnap
    # need to change list of tuples into list of lists. use map              
    sub_pos = np.array(map(list,sub_mb[mask][['posX','posY','posZ']]))[::-1]
    host_pos = np.array(map(list,host_mb[mask][['posX','posY','posZ']]))[::-1]
    return sub_pos-host_pos
    """

def getPericenterInfo(orbital_dist,sub_mb,nsnaps=320):
    localMins = np.r_[True, orbital_dist[1:] < orbital_dist[:-1]] & np.r_[orbital_dist[:-1] < orbital_dist[1:], True]
    snap_list = sub_mb['snap'][0:len(orbital_dist)][::-1]
    perisnaps = snap_list[localMins]
    peridists = orbital_dist[localMins]

    perisnap = perisnaps[-1]
    peridist = peridists[-1]
    # if perisnap is z=0 snapshot, instead choose the one true minimum before it   
    # but what if z=0 min is less than previous pericenter?    
    if perisnap == nsnaps-1:
        if len(perisnaps)>1:
            if peridists[-2]<peridist:
                perisnap = perisnaps[-2]
                peridist = peridists[-2]
    return perisnap, np.min(peridists), np.mean(peridists)



# 1232164  this halo failed for some reason on analyze. this is halo 20, hpaths[19]

def getInfall(sub_mb, host_mb, maxmass=''):
    assert(host_mb[0]['snap']==sub_mb[0]['snap']), "ERROR: Mismatch in alignment of host_mb and sub_mb"

    host_ids = host_mb['id']
    sub_upids = sub_mb['upid']
    if len(host_ids) < len(sub_upids):
        still_sub = np.where(host_ids == sub_upids[0:len(host_ids)])[0]
    else:
        still_sub = np.where(host_ids[0:len(sub_upids)] == sub_upids)[0]
    if len(still_sub) ==0:
        #print 'ERROR: "subhalo" never actually a subhalo. Mass is '+str(maxmass)
        return None, None
    if still_sub[-1] == len(sub_upids)-1:
        #print 'subhalo began as a subhalo. Mass is '+str(maxmass)
        return None, None # "subhalo" began as a subhalo
    else:
        loc = still_sub[-1]+1 #position of infall in array
        # tagging it right before it enters virial radius of host
        iSnap = sub_mb[loc]['snap']
    if sub_mb[loc]['phantom']!=0:
        #phantom halo in merger tree. Go forward one snapshot
        phantom = sub_mb[loc]['phantom']
        loc-=phantom
        iSnap+=phantom
        if loc<0 or sub_mb[loc]['phantom']!=0:
            #print 'subhalo is phantom too much. Mass is '+str(maxmass)
            return None, None
        else:
            dummy=1
            #print 'encountered phantom, but ok'
    return loc, iSnap


# for just getting extant data
class ExtantDataReionization(PluginBase):
    def __init__(self):
        super(ExtantDataReionization,self).__init__()
        self.filename='ExtantDataReionizationVmax.dat'

    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
        mto = MT_Object(hpath)
        subs = mto.cat.get_all_subhalos_within_halo(mto.hostID)
        host_mb = mto.hosttree.getMainBranch(0, mto.mmp_map)
        good = 0; toosmall=0
        for subRSID, subrank in zip(np.array(subs['id']),np.arange(1,len(subs)+1)):
            if mto.mtc.Trees.has_key(subRSID):
                mto.searchtree = mto.mtc.Trees[subRSID] # can it be none?
            else:
                #print 'subhalo',subRSID,'not in merger tree'
                continue
            mto.mmp_map = mto.searchtree.get_mmp_map()
            mto.non_mmp_map = mto.searchtree.get_non_mmp_map()
            sub_mb = mto.searchtree.getMainBranch(0,mto.mmp_map)      
            if sub_mb is None:
                print 'subhalo', sub_rank, 'main branch not found in MT. Skipping it. Z=0 Mass: %.4e, Vmax: %.4f' %(mto.cat.ix[subRSID]['mgrav'], mto.cat.ix[subRSID]['vmax'])
                sys.stdout.flush()
                continue # skip to next subhalo

            # Get maximal mass. Use this to ignore the small halos.
            max_mass = np.max(sub_mb['mvir'])
            if max_mass/mto.cat.h0 < mto.min_mass:
                toosmall+=1
                continue

            # get infall time, if possible
            iLoc, iSnap = getInfall(sub_mb, host_mb, max_mass) 
            add_data(mto,sub_mb,iLoc,subrank,depth=0, host_mb=host_mb, iSnap=iSnap)
            # now get all things that merged into this extant sub
            #print 'halo',subrank,'in host level 0, added of depth 0'
            if iLoc is None:
                iLoc=len(sub_mb)-1 # only for deciding new end_level
            auxiliary_add(mto,host_mb[1:],cur_line=0,level=0, end_level=iLoc, subrank=-abs(subrank), depth=1)

            if subrank%50==0 or subrank==1:
                print subrank, '/', len(subs), 'finished. Time = ', (time.time()-mto.start_time)/60., 'minutes'
            sys.stdout.flush()
            good+=1
        g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,'wb')
        np.array(mto.otherdata).tofile(g)
        g.close()
        # to test, use haloutils.get_hpath_lx(hid,lx)
        print good, 'halos good out of', len(subs)
        print toosmall, 'num halos too small'
        print len(subs)-good-toosmall, 'number of subhalo failures'
        print 'All finished. Time = ', (time.time()-mto.start_time)/60., 'minutes'

    # try for reionization redshifts at 9.31,10.71,11.97,13.50  (Sawala did 11.5)
    # closest snapshots are: 9.33- 59,   10.73 - 50,  12.1 - 43,  13.57 - 37
        
    # now I want 14.25, 13.25, 12.25, 11.25, 10.25, 9.25
    # closest snapshots are: 

    # add_data is where these snapshot numbers are specified.
    # to get a list of snapshots and redshifts: 
    # hpaths = dm.get_hpaths(field=False)
    # snaps = np.arange(20,95)
    # a = zip(htils.get_z_snap(hpaths[0], snaps), snaps)
    # try for new snapshots
    def _read(self,hpath):
        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename)
        #pdtype = ['sub_rank_reion','rsid_reion','depth_reion','vmax_9','vmax_11','vmax_12','vmax_13']
        #pdtype = ['sub_rank_reion','rsid_reion','depth_reion','vmax_8', 'vmax_9','vmax_10','vmax_11','vmax_12','vmax_13','vmax_14','m200_8', 'm200_9','m200_10', 'm200_11', 'm200_12', 'm200_13','m200_14','tvir_8', 'tvir_9','tvir_10', 'tvir_11', 'tvir_12', 'tvir_13','tvir_14', 'max_mass_full', 'max_mass_half', 'max_mass_third', 'max_mass_fourth', 'min_peri', 'mean_peri']
        pdtype = ['sub_rank_reion','rsid_reion','depth_reion','vmax_6','vmax_7','vmax_8', 'vmax_9','vmax_10','vmax_11','vmax_12','vmax_13','vmax_14','m200_6','m200_7','m200_8', 'm200_9','m200_10', 'm200_11', 'm200_12', 'm200_13','m200_14','tvir_6','tvir_7','tvir_8', 'tvir_9','tvir_10', 'tvir_11', 'tvir_12', 'tvir_13','tvir_14', 'max_mass_full', 'max_mass_half', 'max_mass_third', 'max_mass_fourth', 'min_peri', 'mean_peri']
        n = len(pdtype)
        import pandas
        return pandas.DataFrame(data.reshape(len(data)/n,n), columns=pdtype)

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return

class MT_Object(object):
    def __init__(self,hpath):
        print 'creating MT object'
        if hpath=="/bigbang/data/AnnaGroup/caterpillar/halos/H1387186/H1387186_EB_Z127_P7_LN7_LX15_O4_NV4":
            self.cat = haloutils.load_rscat(hpath,165,rmaxcut=False)
            print int(self.cat['id'][0:1]), 'what biggest halo is?'
            print self.cat.ix[int(self.cat['id'][0:1])]['mgrav']/1.e12, 'mass of this halo e12'
            self.hostID = 648070
            print '648070 what brendan says host halo is'
            print self.cat.ix[648070]['mgrav']/1.e12, 'mass of brendans selection e12'
            
        else:
            print 'hpath not the level 15 run'
            snap_z0 = haloutils.get_numsnaps(hpath)-1
            self.cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
            self.hostID = haloutils.load_zoomid(hpath)
        
        self.hosthalo = self.cat.ix[self.hostID]
        self.mtc = haloutils.load_mtc(hpath,haloids=[self.hostID],indexbyrsid=True)
        print 'loaded mtc'
        self.hosttree = self.mtc.Trees[self.hostID]
        self.searchtree = self.hosttree
        self.mmp_map = self.searchtree.get_mmp_map()
        self.non_mmp_map = self.searchtree.get_non_mmp_map()
        self.min_mass = 10**7.5
        self.otherdata = []
        self.start_time = time.time()
        self.hpath = hpath

def Tvir(Mvir,z, delta=200):
    # try to use m200 for this
    h0 = 0.6711
    omega_m = 0.3125
    omega_lambda = 1-omega_m
    def f(x):
        return omega_m * (1+x)**3 + omega_lambda
    Tpower = 4.4 * 10**4 * (delta/200 )**(1/3.) * (Mvir*h0 /(10**10))**(2.0/3) * f(z)
    Tgriff = ( Mvir / (18*3.14*3.14/178/omega_m)**0.5 / (0.6*10/1.22/1.98e4)**1.5/ (1.0e8/h0) )**(2./3) *(1.0+z)
    return Tgriff        # Tpower, Tgriff
    # M = 1e8/h * (T/(1.+z))**1.5 * (0.6*10/1.22/1.98e4)**1.5 * (18*3.14*3.14/178/omega_m)**0.5 #in solar masses
    # mass without h0



def add_data(mto,sub_mb,iLoc,subrank,depth, host_mb, iSnap):
    # can use sub_mb and host_mb to get pericenter here
    min_peri, mean_peri = get_pericenter(sub_mb, iSnap, host_mb,mto.cat)
    #print min_peri, mean_peri, 'min and mean peri'

    # need to get index loc of the snapshots  59, 50, 43, 37
    # do I care about phantom values?
    # z = 9.33, 10.73, 12.1, 13.57 at these snapshots

    # max mass with every other position. test snapshot frequency convergence
    # phantom halo in merger tree. Find peak of non phantom values   
    mask = np.where(sub_mb['phantom']==0)[0]
    tmploc = np.argmax(sub_mb[mask]['mvir'])
    tmploc_2 = 2*np.argmax(sub_mb[mask]['mvir'][::2])
    tmploc_3 = 3*np.argmax(sub_mb[mask]['mvir'][::3])
    tmploc_4 = 4*np.argmax(sub_mb[mask]['mvir'][::4])
    max_mass_loc = mask[tmploc]
    max_mass_loc_2 = mask[tmploc_2]
    max_mass_loc_3 = mask[tmploc_3]
    max_mass_loc_4 = mask[tmploc_4]
    max_mass_mvir = sub_mb[max_mass_loc]['mvir']
    max_mass_mvir_2 = sub_mb[max_mass_loc_2]['mvir']
    max_mass_mvir_3 = sub_mb[max_mass_loc_3]['mvir']
    max_mass_mvir_4 = sub_mb[max_mass_loc_4]['mvir']


    h0 = 0.6711
    snaps = sub_mb['snap']  # sub mb counts from snap = 200 to snap = 0
    zsnaps = [90,78,67,59,53,47,42,38,34]  # corresponds to z = 6.33, 7.26, 8.346, 9.33, 10.22, 11.28, 12.33, 13.31, 14.44
    vmaxes = [0]*len(zsnaps); tvirs = [0]*len(zsnaps); m200s = [0]*len(zsnaps)
    
    for i in range(len(zsnaps)):
        cursnap = zsnaps[i]
        if snaps[-1] > cursnap or snaps[0]<cursnap:
            vmaxes[i] = 0; tvirs[i] = 0; m200s[i] = 0
        else:
            loc = np.where(snaps==cursnap)[0][0]
            vmaxes[i] = np.max(sub_mb[loc:]['vmax'])
            m200 = sub_mb[loc:]['m200c_all']    #['altm2']
            mvir = sub_mb[loc:]['mvir']  
            snapshots = sub_mb[loc:]['snap']
            scales = haloutils.get_scale_snap(mto.hpath,snapshots)
            z = 1.0/scales - 1
            tvirs[i] = np.max(Tvir(mvir/h0, z))  # need to test that this functions and gets reasonable values. So far it does not
            m200s[i] = np.max(m200)
            
    """
    if snaps[-1] > 37 or snaps[0]<37:
        vmax_13 = 0; tvir_13 = 0; m200_13 = 0
    else:
        loc37 = np.where(snaps==37)[0][0]
        vmax_13 = np.max(sub_mb[loc37:]['vmax'])
        m200 = sub_mb[loc37:]['m200c_all']    #['altm2']
        mvir = sub_mb[loc37:]['mvir']  
        snapshots = sub_mb[loc37:]['snap']
        scales = haloutils.get_scale_snap(mto.hpath,snapshots)
        z = 1.0/scales - 1
        tvir_13 = np.max(Tvir(mvir/h0, z))  # need to test that this functions and gets reasonable values. So far it does not
        m200_13 = np.max(m200)
    """ 
    #pdtype = ['sub_rank_reion','rsid_reion','depth_reion','vmax_9','vmax_11','vmax_12','vmax_13', 'm200_9', 'm200_11', 'm200_12', 'm200_13', 'tvir_9', 'tvir_11', 'tvir_12', 'tvir_13']        
    # all values should be the mx on the merger tree before the fixed redshift
    # 'vmax_9','vmax_11','vmax_12','vmax_13'
    # 'm200_9', 'm200_11', 'm200_12', 'm200_13'     # from madau, m200 = 10^6 is about H2 cooling level.
    # 'Tvir_9', 'Tvir_11', 'Tvir_12','Tvir_13' 10^3 Kelvin is about H2 cooling level for the first halos
    
    #I can just looop over vmax, m, t arrays. vmax[0], vmax[1], etc.

    mto.otherdata=np.r_[mto.otherdata,subrank,sub_mb['origid'][0],depth, vmaxes[0],vmaxes[1],vmaxes[2],vmaxes[3],vmaxes[4],vmaxes[5],vmaxes[6],vmaxes[7],vmaxes[8], m200s[0], m200s[1], m200s[2], m200s[3], m200s[4], m200s[5],m200s[6],m200s[7],m200s[8], tvirs[0],tvirs[1],tvirs[2],tvirs[3],tvirs[4],tvirs[5],tvirs[6],tvirs[7],tvirs[8],max_mass_mvir, max_mass_mvir_2, max_mass_mvir_3, max_mass_mvir_4, min_peri, mean_peri]
    return


def auxiliary_add(mto, host_mb, cur_line, level, end_level, subrank,depth,destr=False):
    while level!=end_level:
        merged_subs = mto.searchtree.getNonMMPprogenitors(cur_line, mto.non_mmp_map) # merged_subs are one step up
        for subline in merged_subs:
            # now get main branch of this halo
            sub_mb = mto.searchtree.getMainBranch(subline, mto.mmp_map)
            max_mass = np.max(sub_mb['mvir'])
            if max_mass/mto.cat.h0 < mto.min_mass: 
                continue # skip if too small
            #print 'adding sub with mass', np.log10(max_mass/mto.cat.h0)
            
            # get infall time, if possible
            iLoc, iSnap = getInfall(sub_mb,host_mb, max_mass)# if None, still write it

            add_data(mto,sub_mb, iLoc,subrank,depth, host_mb, iSnap)
            #print 'halo',subrank,', in host level', level, ', added of depth', depth
            sys.stdout.flush()

            if iLoc is None:
                iLoc=len(sub_mb)-1 # only for deciding new end_level
            # go recursively deep
            auxiliary_add(mto,host_mb[1:],cur_line=subline,level=level, end_level=level+iLoc, subrank=-abs(subrank), depth=depth+1) # subrank < 0 always if from a merged sub
            if subrank > 0 and destr:
                subrank+=1
        level+=1
        cur_line = mto.searchtree.getMMP(cur_line,mto.mmp_map)
        host_mb=host_mb[1:]
        if subrank > 0 and destr:
            print 'finished level', level, 'Time = ',(time.time()-mto.start_time)/60., 'minutes'
    return
