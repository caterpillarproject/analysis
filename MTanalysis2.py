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


"""
hid = 'H1631506'
hpath = haloutils.get_hpath_lx(hid,14)
subRSID = 148923 # this is a sub-sub halo!! And its parent is not in the merger tree....

mass = 2.13e8
sub_rank = 46
flags = 9

its hostID is:  162183, which is a subhalo. So it is a sub-sub halo.
162183 is also not found in the merger tree
1.06 kpc away from center of host, mgrav is 3.032336896e9. mvir is 4.8e9

162185 is the main host

hid = 'H581141'
hpath = haloutils.get_hpath_lx(hid,14)
subRSID = 118959
its hostID is 119299
main hostID is 155162
--> it is a sub-sub again
its first host is found in the merger tree. at a distance of 130.7 kpc

hid = 'H1387186'
hpath = haloutils.get_hpath_lx(hid,14)
subRSID = 119955
its host id is: 150793
real main halo id: 150793
This halo is in the merger tree as a destroyed halo
it is dataD[682]
infall snap is 175.
backsnap is 19
infall scale is 0.4473


hid = 'H581180'
hpath = haloutils.get_hpath_lx(hid,14)
subRSID = 166360
its host is 166371
main host is 166371
flags = 9. same as the subhalos that do work
"""

# for quick testing:
#hpath = '/bigbang/data/AnnaGroup/caterpillar/halos/H1599988/H1599988_EX_Z127_P7_LN7_LX14_O4_NV4'

#/bigbang/data/AnnaGroup/caterpillar/halos/H95289/H95289_BB_Z127_P7_LN7_LX14_O4_NV4
# subhalo ids:
#44661. Mass in RSCatalog way too high. Mass in merger tree very reasonsble.
#its host is 42818, which is much smaller than 44661

# trying to push
#  git config --global push.default matching
#  git config --global push.default simple

######################### EXAMPLE TO LOAD STAR IDS AND MASS
#AE = AllExtantData()
#AD = AllDestroyedData()
#TM = TagMass()
# dataE = AE.read(hpath)
# dataD = AD.read(hpath)
#idsE,massE,idsD,massD = TM.read(hpath)
# get stars from specific halo:
# halostars = getStars(dataE, idsE, row)
# halostarmass - getStarMass(dataE, massE,row)
########################################


# tagging 1% - because start_pos is not sequential,
# might want to 
# quick solution - use TagMass to get ids, mass for getStars_1
# but use TagMass_1 to get ids1, etc.

# better - re-run AllExtantData
# use 'nstars_1', 'start_pos_1'
# instead of writing to file, use generate ids, mass

# make it possible to get any fraction of 3%



def getStars(data, ids, row):
    sp = data['start_pos'][row]
    nstars = data['nstars'][row]
    return ids[sp:sp+nstars]

def getStars_x(data,ids,row, fmb=1):
    if fmb  > 3:
        print "ERROR fmb cannot be > 3"
        return None
    sp = data['start_pos'][row]
    nstars = np.round(data['nstars'][row]/(3./fmb))
    if nstars==0 and data['nstars'][row]>0:
        nstars=1
    return ids[sp:sp+nstars]

def getStarMass(data, mass, row):
    sp = data['start_pos'][row]
    nstars = data['nstars'][row]
    return mass[sp:sp+nstars]

def getStarMass_x(data, mass, row, fmb=1):
    if fmb  > 3:
        print "ERROR fmb cannot be > 3"
        return None
    sp = data['start_pos'][row]
    nstars = np.round(data['nstars'][row]/(3./fmb))
    if nstars==0 and data['nstars'][row]>0:
        nstars=1
    return mass[sp:sp+nstars]


def getIdsMass_x(hpath, fmb=1):
    TM = TagMass()
    idsE,massE,idsD,massD = TM.read(hpath)
    AE = AllExtantData()
    AD = AllDestroyedData()
    dataE = AE.read(hpath)
    dataD = AD.read(hpath)
    
    idsE_1 = []; mper_arrE = []
    for row in xrange(len(dataE)):
        idsE_1 = np.r_[idsE_1, getStars_x(dataE,idsE,row,fmb)]
        mper_arrE = np.r_[mper_arrE, getStarMass_x(dataE,massE,row,fmb)]

    idsD_1 = []; mper_arrD = []
    for row in xrange(len(dataD)):
        idsD_1 = np.r_[idsD_1, getStars_x(dataD,idsD,row,fmb)]
        mper_arrD = np.r_[mper_arrD, getStarMass_x(dataD,massD,row,fmb)]
    return np.array(idsE_1), np.array(mper_arrE), np.array(idsD_1), np.array(mper_arrD)


def distance(posA, posB,boxsize=100.):
    dist = abs(posA-posB)
    tmp = dist > boxsize/2.0
    dist[tmp] = boxsize-dist[tmp]
    if dist.shape == (3,):
        return np.sqrt(np.sum(dist**2))
    else:
        return np.sqrt(np.sum(dist**2,axis=1))

def getInfall(sub_mb, host_mb, maxmass=''):
    host_ids = host_mb['id']
    sub_upids = sub_mb['upid']
    if len(host_ids) < len(sub_upids):
        still_sub = np.where(host_ids == sub_upids[0:len(host_ids)])[0]
    else:
        still_sub = np.where(host_ids[0:len(sub_upids)] == sub_upids)[0]
    if len(still_sub) ==0:
        print 'ERROR: "subhalo" never actually a subhalo. Mass is '+str(maxmass)
        return None, None
    if still_sub[-1] == len(sub_upids)-1:
        print 'subhalo began as a subhalo. Mass is '+str(maxmass)
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
            print 'subhalo is phantom too much. Mass is '+str(maxmass)
            return None, None
        else:
            print 'encountered phantom, but ok'
    return loc, iSnap

def getSubTree(mtc,rsid, hostrow=0):
    """ find tree of subhalo with id = rsid 
    subhalo is a member of host belonging to row=hostrow """
    for i in mtc.getSubTrees(hostrow):
        if mtc.Trees[i].rockstar_id == rsid:
            return mtc.Trees[i]
    print 'halo with RSID =', rsid, 'not found in mt' 
    return None


# for just getting extant data
class ExtantDataFirstPass(PluginBase):
    def __init__(self):
        super(ExtantDataFirstPass,self).__init__()
        self.filename='ExtantDataFirstPass.dat'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = True
        self.xlabel='scale factor' ; self.ylabel='Mass Accreted'    # want these to be adjustable
        self.autofigname='MergerHistory'
        self.min_mass = 10**6.0

    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
        # copy tagExtant code here
        start = 0;
        start_time = time.time()
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID1 = int(cat['id'][0:1])
        hostID = haloutils.load_zoomid(hpath)
        if hostID != hostID1:
            print 'host IDs do not match!!'
        hosthalo = cat.ix[hostID]
        subs = cat.get_all_subhalos_within_halo(hostID)

        otherdata=[]
        print 'loading mtc'
        
        sys.stdout.flush()
        mtc = haloutils.load_mtc(hpath,haloids=[hostID])
        print 'loaded mtc'
        sys.stdout.flush()
        host = mtc.Trees[0]
        host_mb = host.getMainBranch(0)
        
        good = 0; toosmall=0; sub_rank=start-1
        for subRSID in np.array(subs['id']):
            sub_rank+=1
            sub = getSubTree(mtc,subRSID)
            if sub==None:
                print sub_rank, 'subhalo not found in MTCatalogue. Z=0 Mass: %.4e, Vmax: %.4f' %(cat.ix[subRSID]['mgrav'], cat.ix[subRSID]['vmax']), 'Time = ', (time.time()-start_time)/60., 'minutes'
                sys.stdout.flush()
                continue
            sub_mb = sub.getMainBranch(0)
            if sub_mb == None:
                print 'subhalo', sub_rank, 'main branch not found in MT. Skipping it. Z=0 Mass: %.4e, Vmax: %.4f' %(cat.ix[subRSID]['mgrav'], cat.ix[subRSID]['vmax'])
                sys.stdout.flush()
                continue # skip to next subhalo

            # Get maximal mass. Use this to ignore the small halos.
            max_mass = np.max(sub_mb['mvir'])
            if max_mass/cat.h0 < self.min_mass: 
                #print sub_rank, 'subhalo too small'
                sys.stdout.flush()
                toosmall+=1
                continue


            # get infall time, if possible
            iLoc, iSnap = getInfall(sub_mb, host_mb, max_mass)
            if iLoc==None:
                # within getInfall, print how it errored
                continue

            # get max_mass values.
            max_mass_loc = np.argmax(sub_mb['mvir'])
            if sub_mb[max_mass_loc]['phantom']!=0:
                # phantom halo in merger tree. Find peak of non phantom values
                mask = np.where(sub_mb['phantom']==0)[0]
                tmploc = np.argmax(sub_mb[mask]['mvir'])
                max_mass_loc = mask[tmploc]

            max_mass_vmax = sub_mb[max_mass_loc]['vmax']
            max_mass_snap = sub_mb[max_mass_loc]['snap']
            max_mass_rsid = sub_mb[max_mass_loc]['origid']
            max_mass_mvir = sub_mb[max_mass_loc]['mvir']
            max_mass_posx = sub_mb[max_mass_loc]['posX']
            max_mass_posy = sub_mb[max_mass_loc]['posY']
            max_mass_posz = sub_mb[max_mass_loc]['posZ']
            max_mass_pecvx = sub_mb[max_mass_loc]['pecVX']
            max_mass_pecvy = sub_mb[max_mass_loc]['pecVY']
            max_mass_pecvz = sub_mb[max_mass_loc]['pecVZ']
            max_mass_virialratio = sub_mb[max_mass_loc]['T/|U|']
            max_mass_hostid_MT = sub_mb[max_mass_loc]['pid'] # merger tree ID of host, one level up
            max_mass_rvir = sub_mb[max_mass_loc]['rvir']
            max_mass_spinbullock = sub_mb[max_mass_loc]['spin_bullock']
            max_mass_rs = sub_mb[max_mass_loc]['rs']
            max_mass_scale_of_last_MM = sub_mb[max_mass_loc]['scale_of_last_MM']
            max_mass_Jx = sub_mb[max_mass_loc]['Jx']
            max_mass_Jy = sub_mb[max_mass_loc]['Jy']
            max_mass_Jz = sub_mb[max_mass_loc]['Jz']
            max_mass_xoff = sub_mb[max_mass_loc]['xoff']
            

            # get all peak values. Peak values all based on when vmax reaches its peak.
            peak_loc = np.argmax(sub_mb['vmax'])
            if sub_mb[peak_loc]['phantom']!=0:
                # phantom halo in merger tree. Find peak of non phantom values
                mask = np.where(sub_mb['phantom']==0)[0]
                tmploc = np.argmax(sub_mb[mask]['vmax'])
                peak_loc = mask[tmploc]

            peak_vmax = sub_mb[peak_loc]['vmax']
            peak_snap = sub_mb[peak_loc]['snap']
            peak_rsid = sub_mb[peak_loc]['origid']
            peak_mvir = sub_mb[peak_loc]['mvir']
            peak_posx = sub_mb[peak_loc]['posX']
            peak_posy = sub_mb[peak_loc]['posY']
            peak_posz = sub_mb[peak_loc]['posZ']
            peak_pecvx = sub_mb[peak_loc]['pecVX']
            peak_pecvy = sub_mb[peak_loc]['pecVY']
            peak_pecvz = sub_mb[peak_loc]['pecVZ']
            peak_virialratio = sub_mb[peak_loc]['T/|U|']
            peak_hostid_MT = sub_mb[peak_loc]['pid'] # merger tree ID of host, one level up
            peak_rvir = sub_mb[peak_loc]['rvir']
            peak_spinbullock = sub_mb[peak_loc]['spin_bullock']
            peak_rs = sub_mb[peak_loc]['rs']
            peak_scale_of_last_MM = sub_mb[peak_loc]['scale_of_last_MM']
            peak_Jx = sub_mb[peak_loc]['Jx']
            peak_Jy = sub_mb[peak_loc]['Jy']
            peak_Jz = sub_mb[peak_loc]['Jz']
            peak_xoff = sub_mb[peak_loc]['xoff']
            
            # Get infall parameters
            infall_snap = sub_mb[iLoc]['snap']
            #infall_scale = sub_mb[iLoc]['scale']
            infall_rsid = sub_mb[iLoc]['origid']
            infall_vmax = sub_mb[iLoc]['vmax']
            infall_mvir = sub_mb[iLoc]['mvir']
            infall_posx = sub_mb[iLoc]['posX']
            infall_posy = sub_mb[iLoc]['posY']
            infall_posz = sub_mb[iLoc]['posZ']
            infall_pecvx = sub_mb[iLoc]['pecVX']
            infall_pecvy = sub_mb[iLoc]['pecVY']
            infall_pecvz = sub_mb[iLoc]['pecVZ']
            infall_virialratio = sub_mb[iLoc]['T/|U|']
            infall_hostid_MT = sub_mb[iLoc]['pid']
            infall_rvir = sub_mb[iLoc]['rvir']
            infall_spinbullock = sub_mb[iLoc]['spin_bullock']
            infall_rs = sub_mb[iLoc]['rs']
            infall_scale_of_last_MM = sub_mb[iLoc]['scale_of_last_MM']
            infall_Jx = sub_mb[iLoc]['Jx']
            infall_Jy = sub_mb[iLoc]['Jy']
            infall_Jz = sub_mb[iLoc]['Jz']
            infall_xoff = sub_mb[iLoc]['xoff']

            otherdata=np.r_[otherdata,sub_rank,subRSID, max_mass_rsid, max_mass_snap, max_mass_vmax,max_mass_mvir,max_mass_posx,max_mass_posy,max_mass_posz,max_mass_pecvx,max_mass_pecvy,max_mass_pecvz,max_mass_virialratio,max_mass_hostid_MT,max_mass_rvir,max_mass_spinbullock,max_mass_rs,max_mass_scale_of_last_MM,max_mass_Jx,max_mass_Jy,max_mass_Jz,max_mass_xoff,peak_rsid, peak_snap, peak_vmax,peak_mvir,peak_posx,peak_posy,peak_posz,peak_pecvx,peak_pecvy,peak_pecvz,peak_virialratio,peak_hostid_MT,peak_rvir,peak_spinbullock,peak_rs,peak_scale_of_last_MM,peak_Jx,peak_Jy,peak_Jz,peak_xoff,infall_rsid,infall_snap,infall_vmax,infall_mvir,infall_posx,infall_posy,infall_posz,infall_pecvx,infall_pecvy,infall_pecvz,infall_virialratio,infall_hostid_MT,infall_rvir,infall_spinbullock,infall_rs,infall_scale_of_last_MM,infall_Jx,infall_Jy,infall_Jz,infall_xoff]
            if sub_rank%100==0:
                print sub_rank, '/', len(subs), 'finished. Time = ', (time.time()-start_time)/60., 'minutes'
            sys.stdout.flush()
            good+=1
        g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,'wb')
        np.array(otherdata).tofile(g)
        g.close()
# to test, use haloutils.get_hpath_lx(hid,lx)
        print good, 'halos good out of', len(subs)
        print toosmall, 'num halos too small'
        print len(subs)-good-toosmall, 'number of subhalo failures'

    def _read(self,hpath):
        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename)
        #dt = "float64"
        #dtype = [('sub_rank',dt),('rsid',dt),('max_mass',dt),('max_mass_snap',dt), ('peak_rsid',dt), ('peak_snap',dt), ('peak_vmax',dt),('peak_mvir',dt),('peak_posx',dt),('peak_posy',dt),('peak_posz',dt),('peak_pecvx',dt),('peak_pecvy',dt),('peak_pecvz',dt),('peak_virialratio',dt),('peak_hostid_MT',dt),('peak_rvir',dt),('peak_spinbullock',dt),('infall_rsid',dt),('infall_snap',dt),('infall_vmax',dt),('infall_mvir',dt),('infall_posx',dt),('infall_posy',dt),('infall_posz',dt),('infall_pecvx',dt),('infall_pecvy',dt),('infall_pecvz',dt),('infall_virialratio',dt),('infall_hostid_MT',dt),('infall_rvir',dt),('infall_spinbullock',dt)]
        #n = len(dtype)
        #holder = np.ndarray( (len(data)/n,n), dtype=dtype )
        #data2 = data.reshape(len(data)/n,n)
        #for i in range(data2.shape[0]):
        #    holder[i]=data2[i]
        #return holder
        pdtype = ['sub_rank','rsid','max_mass_rsid','max_mass_snap','max_mass_vmax','max_mass','max_mass_posx','max_mass_posy','max_mass_posz','max_mass_pecvx','max_mass_pecvy','max_mass_pecvz','max_mass_virialratio','max_mass_hostid_MT','max_mass_rvir','max_mass_spinbullock','max_mass_rs','max_mass_scale_of_last_MM','max_mass_Jx','max_mass_Jy','max_mass_Jz','max_mass_xoff','peak_rsid','peak_snap','peak_vmax','peak_mvir','peak_posx','peak_posy','peak_posz','peak_pecvx','peak_pecvy','peak_pecvz','peak_virialratio','peak_hostid_MT','peak_rvir','peak_spinbullock','peak_rs','peak_scale_of_last_MM','peak_Jx','peak_Jy','peak_Jz','peak_xoff','infall_rsid','infall_snap','infall_vmax','infall_mvir','infall_posx','infall_posy','infall_posz','infall_pecvx','infall_pecvy','infall_pecvz','infall_virialratio','infall_hostid_MT','infall_rvir','infall_spinbullock','infall_rs','infall_scale_of_last_MM','infall_Jx','infall_Jy','infall_Jz','infall_xoff']
        n = len(pdtype)
        import pandas
        return pandas.DataFrame(data.reshape(len(data)/n,n), columns=pdtype)

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return


class AllExtantData(PluginBase):
    def __init__(self):
        super(AllExtantData,self).__init__()
        self.filename='AllExtantData.dat'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = True
        self.xlabel='' ; self.ylabel=''  # want these to be adjustable
        self.autofigname=''
        self.min_mass = 10**6.0

    def _analyze(self,hpath):
        ED = ExtantDataFirstPass()
        dataE = ED.read(hpath)
        dtype = ['peak_mgrav','infall_mgrav','peak_hostid_RS','infall_hostid_RS','peak_rvmax','infall_rvmax','peak_corevelx','peak_corevely','peak_corevelz','infall_corevelx','infall_corevely','infall_corevelz', 'nstars', 'start_pos']
        data_newE = pandas.DataFrame(np.zeros((len(dataE),len(dtype)))-1,columns=dtype)
        peak_dataE = {}
        for peaksnap,line in zip(dataE['peak_snap'],dataE.index):
            peak_dataE.setdefault(peaksnap, []).append(line)

        infall_dataE = {}
        for infallsnap,line in zip(dataE['infall_snap'],dataE.index):
            infall_dataE.setdefault(infallsnap, []).append(line)       

        # initialize arrays for tagging
        allstars=[]; start_pos=0     
       
        for snap in range(haloutils.get_numsnaps(hpath)):
            print snap, 'snap in get extra parameters Extant'
            sys.stdout.flush()
            if peak_dataE.has_key(snap) or infall_dataE.has_key(snap):
                cat=haloutils.load_rscat(hpath,snap,rmaxcut=False)

                if peak_dataE.has_key(snap):
                    for line in peak_dataE[snap]:
                        peak_rsid = int(dataE.ix[line]['peak_rsid'])
                        data_newE.ix[line]['peak_mgrav'] = cat.ix[peak_rsid]['mgrav']
                        data_newE.ix[line]['peak_hostid_RS'] = cat.ix[peak_rsid]['hostID']
                        data_newE.ix[line]['peak_rvmax'] = cat.ix[peak_rsid]['rvmax']
                        data_newE.ix[line]['peak_corevelx'] = cat.ix[peak_rsid]['corevelx']
                        data_newE.ix[line]['peak_corevely'] = cat.ix[peak_rsid]['corevely']
                        data_newE.ix[line]['peak_corevelz'] = cat.ix[peak_rsid]['corevelz']
                        
                if infall_dataE.has_key(snap):
                    for line in infall_dataE[snap]:
                        infall_rsid = int(dataE.ix[line]['infall_rsid'])
                        data_newE.ix[line]['infall_mgrav'] = cat.ix[infall_rsid]['mgrav']
                        data_newE.ix[line]['infall_hostid_RS'] = cat.ix[infall_rsid]['hostID']
                        data_newE.ix[line]['infall_rvmax'] = cat.ix[infall_rsid]['rvmax']
                        data_newE.ix[line]['infall_corevelx'] = cat.ix[infall_rsid]['corevelx']
                        data_newE.ix[line]['infall_corevely'] = cat.ix[infall_rsid]['corevely']
                        data_newE.ix[line]['infall_corevelz'] = cat.ix[infall_rsid]['corevelz']
                        
                        iPids = cat.get_all_particles_from_halo(infall_rsid)
                        star_pids = iPids[0:int(np.round(len(iPids)*.03))]
                        data_newE.ix[line]['nstars'] = len(star_pids)
                        data_newE.ix[line]['start_pos'] = start_pos
                        allstars=np.r_[allstars,star_pids]
                        start_pos+=len(star_pids)
                
        fulldataE = pandas.concat((dataE,data_newE),axis=1)
        fulldataE.to_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t')
        f = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'extantPIDs.dat', 'wb')
        np.array(allstars).tofile(f)
        f.close()

    def _read(self,hpath):
        return pandas.read_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t')

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return








# This currently does not identify sub-subs in merger tree
# need to port over code from auxiliary_add in tagDestroyed.py on Odyssey TagCode4.0 folder
# Scenario 1: sub-sub enters main host on its own, then falls into sub, then merges with sub
# Secenario 2: sub-sub falls into sub, then sub and sub-sub system fall into main host together,
# then sub-sub merges with sub
# In both cases, I am currently tagging sub-sub when it first enters the main host.
class DestroyedDataFirstPass(PluginBase):
    def __init__(self):
        super(DestroyedDataFirstPass,self).__init__()
        self.filename='DestroyedDataFirstPass.dat'
        self.filestring='DestroyedDataFirstPass'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = True
        self.xlabel='scale factor' ; self.ylabel='Mass Accreted'    # want these to be adjustable
        self.autofigname='MergerHistory'
        self.min_mass = 10**6.0
        # corresponds to 10**7.776 Msun
       
    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
        start = 0
        start_time = time.time()
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID1 = int(cat['id'][0:1])
        hostID = haloutils.load_zoomid(hpath)
        if hostID != hostID1:
            print 'host IDs do not match!!'
        hosthalo = cat.ix[hostID]
        mtc = haloutils.load_mtc(hpath,haloids=[hostID])
        print 'loaded MTC'
        host = mtc.Trees[0]
        #### determine end ahead of time.
        host_mb = host.getMainBranch(0)
        end = len(host_mb)-1
        ####
        cur_host_line = 0
        i=start # skip to ith level in MT to start
        for k in range(i):
            cur_host_line = host.getMMP(cur_host_line)
        
        while i!=end:
            merged_subs = host.getNonMMPprogenitors(cur_host_line)
            j=-1; start_pos=0; good=0
            otherdata=[]
            host_mb = host.getMainBranch(host.getMMP(cur_host_line)) 
            for subline in merged_subs:
                j+=1
                sub_mb = host.getMainBranch(subline)

                # Get maximal mass. Use this to ignore the small halos.
                max_mass = np.max(sub_mb['mvir'])
                if max_mass/cat.h0 < self.min_mass: 
                    sys.stdout.flush()
                    continue


                # get infall time, if possible
                iLoc, iSnap = getInfall(sub_mb,host_mb, max_mass)
                if iLoc == None:
                    print 'subhalo', j, 'is bad in MT. Reason to follow.'
                    sys.stdout.flush()
                    continue


                # get all max_mass values
                max_mass_loc = np.argmax(sub_mb['mvir'])
                if sub_mb[max_mass_loc]['phantom']!=0:
                    # phantom halo in merger tree. Find peak of non phantom values
                    mask = np.where(sub_mb['phantom']==0)[0]
                    tmploc = np.argmax(sub_mb[mask]['mvir'])
                    max_mass_loc = mask[tmploc]

                max_mass_vmax = sub_mb[max_mass_loc]['vmax']
                max_mass_snap = sub_mb[max_mass_loc]['snap']
                max_mass_rsid = sub_mb[max_mass_loc]['origid']
                max_mass_mvir = sub_mb[max_mass_loc]['mvir']
                max_mass_posx = sub_mb[max_mass_loc]['posX']
                max_mass_posy = sub_mb[max_mass_loc]['posY']
                max_mass_posz = sub_mb[max_mass_loc]['posZ']
                max_mass_pecvx = sub_mb[max_mass_loc]['pecVX']
                max_mass_pecvy = sub_mb[max_mass_loc]['pecVY']
                max_mass_pecvz = sub_mb[max_mass_loc]['pecVZ']
                max_mass_virialratio = sub_mb[max_mass_loc]['T/|U|']
                max_mass_hostid_MT = sub_mb[max_mass_loc]['pid'] # merger tree ID of host, one level up
                max_mass_rvir = sub_mb[max_mass_loc]['rvir']
                max_mass_spinbullock = sub_mb[max_mass_loc]['spin_bullock']
                max_mass_rs = sub_mb[max_mass_loc]['rs']
                max_mass_scale_of_last_MM = sub_mb[max_mass_loc]['scale_of_last_MM']
                max_mass_Jx = sub_mb[max_mass_loc]['Jx']
                max_mass_Jy = sub_mb[max_mass_loc]['Jy']
                max_mass_Jz = sub_mb[max_mass_loc]['Jz']
                max_mass_xoff = sub_mb[max_mass_loc]['xoff']


                # get all peak values. Peak values all based on when vmax reaches its peak.
                peak_loc = np.argmax(sub_mb['vmax'])
                if sub_mb[peak_loc]['phantom']!=0:
                    # phantom halo in merger tree. Find peak of non phantom values
                    mask = np.where(sub_mb['phantom']==0)[0]
                    tmploc = np.argmax(sub_mb[mask]['vmax'])
                    peak_loc = mask[tmploc]
                peak_vmax = sub_mb[peak_loc]['vmax']
                peak_snap = sub_mb[peak_loc]['snap']
                peak_rsid = sub_mb[peak_loc]['origid']
                peak_mvir = sub_mb[peak_loc]['mvir']
                peak_posx = sub_mb[peak_loc]['posX']
                peak_posy = sub_mb[peak_loc]['posY']
                peak_posz = sub_mb[peak_loc]['posZ']
                peak_pecvx = sub_mb[peak_loc]['pecVX']
                peak_pecvy = sub_mb[peak_loc]['pecVY']
                peak_pecvz = sub_mb[peak_loc]['pecVZ']
                peak_virialratio = sub_mb[peak_loc]['T/|U|']
                peak_hostid_MT = sub_mb[peak_loc]['pid'] # merger tree ID of host, one level up
                peak_rvir = sub_mb[peak_loc]['rvir']
                peak_spinbullock = sub_mb[peak_loc]['spin_bullock']
                peak_rs = sub_mb[peak_loc]['rs']
                peak_scale_of_last_MM = sub_mb[peak_loc]['scale_of_last_MM']
                peak_Jx = sub_mb[peak_loc]['Jx']
                peak_Jy = sub_mb[peak_loc]['Jy']
                peak_Jz = sub_mb[peak_loc]['Jz']
                peak_xoff = sub_mb[peak_loc]['xoff']

                # Get infall parameters
                infall_snap = sub_mb[iLoc]['snap']
                #infall_scale = sub_mb[iLoc]['scale']
                infall_rsid = sub_mb[iLoc]['origid']
                infall_vmax = sub_mb[iLoc]['vmax']
                infall_mvir = sub_mb[iLoc]['mvir']
                infall_posx = sub_mb[iLoc]['posX']
                infall_posy = sub_mb[iLoc]['posY']
                infall_posz = sub_mb[iLoc]['posZ']
                infall_pecvx = sub_mb[iLoc]['pecVX']
                infall_pecvy = sub_mb[iLoc]['pecVY']
                infall_pecvz = sub_mb[iLoc]['pecVZ']
                infall_virialratio = sub_mb[iLoc]['T/|U|']
                infall_hostid_MT = sub_mb[iLoc]['pid']
                infall_rvir = sub_mb[iLoc]['rvir']
                infall_spinbullock = sub_mb[iLoc]['spin_bullock']
                infall_rs = sub_mb[iLoc]['rs']
                infall_scale_of_last_MM = sub_mb[iLoc]['scale_of_last_MM']
                infall_Jx = sub_mb[iLoc]['Jx']
                infall_Jy = sub_mb[iLoc]['Jy']
                infall_Jz = sub_mb[iLoc]['Jz']
                infall_xoff = sub_mb[iLoc]['xoff']


                otherdata=np.r_[otherdata,j,sub_mb['origid'][0],max_mass_rsid, max_mass_snap, max_mass_vmax,max_mass_mvir,max_mass_posx,max_mass_posy,max_mass_posz,max_mass_pecvx,max_mass_pecvy,max_mass_pecvz,max_mass_virialratio,max_mass_hostid_MT,max_mass_rvir,max_mass_spinbullock,max_mass_rs,max_mass_scale_of_last_MM,max_mass_Jx,max_mass_Jy,max_mass_Jz,max_mass_xoff, peak_rsid, peak_snap, peak_vmax,peak_mvir,peak_posx,peak_posy,peak_posz,peak_pecvx,peak_pecvy,peak_pecvz,peak_virialratio,peak_hostid_MT,peak_rvir,peak_spinbullock,peak_rs,peak_scale_of_last_MM,peak_Jx,peak_Jy,peak_Jz,peak_xoff,infall_rsid,infall_snap,infall_vmax,infall_mvir,infall_posx,infall_posy,infall_posz,infall_pecvx,infall_pecvy,infall_pecvz,infall_virialratio,infall_hostid_MT,infall_rvir,infall_spinbullock,infall_rs,infall_scale_of_last_MM,infall_Jx,infall_Jy,infall_Jz,infall_xoff,i]
                #print j, 'halo in host level', i
                good+=1
                sys.stdout.flush()
            print i, 'host level finished. Time = ', (time.time()-start_time)/60., 'minutes'
            print good,'/',j+1,'were tagged'
            sys.stdout.flush()
            if not os.path.exists(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed'):
                subprocess.call("mkdir -p "+hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed',shell=True)
            g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/'+self.filestring+'_'+str(i)+'.dat','wb')
            np.array(otherdata).tofile(g)
            g.close()    
            cur_host_line = host.getMMP(cur_host_line)
            i+=1
        print 'wrote final set of data'
        self.combinefiles(hpath)
        print 'combined files'

## convert all data into one file
    def combinefiles(self,hpath):
            i = 0; data=[]
            while os.path.exists(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/'+self.filestring+'_'+str(i)+'.dat'):
                tmp = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/'+self.filestring+'_'+str(i)+'.dat')
                data = np.r_[data,tmp]
                i+=1
            dt = "float64"
            dtype = [('sub_rank',dt),('rsid',dt), ('max_mass_rsid',dt), ('max_mass_snap',dt), ('max_mass_vmax',dt),('max_mass',dt),('max_mass_posx',dt),('max_mass_posy',dt),('max_mass_posz',dt),('max_mass_pecvx',dt),('max_mass_pecvy',dt),('max_mass_pecvz',dt),('max_mass_virialratio',dt),('max_mass_hostid_MT',dt),('max_mass_rvir',dt),('max_mass_spinbullock',dt),('max_mass_rs',dt),('max_mass_scale_of_last_MM',dt),('max_mass_Jx',dt),('max_mass_Jy',dt),('max_mass_Jz',dt),('max_mass_xoff',dt), ('peak_rsid',dt), ('peak_snap',dt), ('peak_vmax',dt),('peak_mvir',dt),('peak_posx',dt),('peak_posy',dt),('peak_posz',dt),('peak_pecvx',dt),('peak_pecvy',dt),('peak_pecvz',dt),('peak_virialratio',dt),('peak_hostid_MT',dt),('peak_rvir',dt),('peak_spinbullock',dt),('peak_rs',dt),('peak_scale_of_last_MM',dt),('peak_Jx',dt),('peak_Jy',dt),('peak_Jz',dt),('peak_xoff',dt), ('infall_rsid',dt),('infall_snap',dt),('infall_vmax',dt),('infall_mvir',dt),('infall_posx',dt),('infall_posy',dt),('infall_posz',dt),('infall_pecvx',dt),('infall_pecvy',dt),('infall_pecvz',dt),('infall_virialratio',dt),('infall_hostid_MT',dt),('infall_rvir',dt),('infall_spinbullock',dt),('infall_rs',dt),('infall_scale_of_last_MM',dt),('infall_Jx',dt),('infall_Jy',dt),('infall_Jz',dt),('infall_xoff',dt),('backsnap',dt)]
            n = len(dtype)
            holder = np.ndarray( (len(data)/n,), dtype=dtype )
            data2 = data.reshape(len(data)/n,n)
            for j in range(data2.shape[0]):
                holder[j]=data2[j]
            g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,'wb')
            np.array(holder).tofile(g)
            g.close()    

            #np.array(bound,dtype=np.float32)
   
    def _read(self,hpath):
        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename)
        pdtype = ['sub_rank','rsid','max_mass_rsid','max_mass_snap','max_mass_vmax','max_mass','max_mass_posx','max_mass_posy','max_mass_posz','max_mass_pecvx','max_mass_pecvy','max_mass_pecvz','max_mass_virialratio','max_mass_hostid_MT','max_mass_rvir','max_mass_spinbullock','max_mass_rs','max_mass_scale_of_last_MM','max_mass_Jx','max_mass_Jy','max_mass_Jz','max_mass_xoff', 'peak_rsid','peak_snap','peak_vmax','peak_mvir','peak_posx','peak_posy','peak_posz','peak_pecvx','peak_pecvy','peak_pecvz','peak_virialratio','peak_hostid_MT','peak_rvir','peak_spinbullock','peak_rs','peak_scale_of_last_MM','peak_Jx','peak_Jy','peak_Jz','peak_xoff','infall_rsid','infall_snap','infall_vmax','infall_mvir','infall_posx','infall_posy','infall_posz','infall_pecvx','infall_pecvy','infall_pecvz','infall_virialratio','infall_hostid_MT','infall_rvir','infall_spinbullock','infall_rs','infall_scale_of_last_MM','infall_Jx','infall_Jy','infall_Jz','infall_xoff','backsnap']
        n = len(pdtype)
        import pandas
        return pandas.DataFrame(data.reshape(len(data)/n,n), columns=pdtype)
    
        """
        dt = "float64"
        dtype = [('sub_rank',dt),('rsid',dt),('max_mass',dt),('max_mass_snap',dt), ('peak_rsid',dt), ('peak_snap',dt), ('peak_vmax',dt),('peak_mvir',dt),('peak_posx',dt),('peak_posy',dt),('peak_posz',dt),('peak_pecvx',dt),('peak_pecvy',dt),('peak_pecvz',dt),('peak_virialratio',dt),('peak_hostid_MT',dt),('peak_rvir',dt),('peak_spinbullock',dt),('infall_rsid',dt),('infall_snap',dt),('infall_vmax',dt),('infall_mvir',dt),('infall_posx',dt),('infall_posy',dt),('infall_posz',dt),('infall_pecvx',dt),('infall_pecvy',dt),('infall_pecvz',dt),('infall_virialratio',dt),('infall_hostid_MT',dt),('infall_rvir',dt),('infall_spinbullock',dt),('backsnap',dt)]
        n = len(dtype)
        holder = np.ndarray( (len(data)/n,n), dtype=dtype )
        data2 = data.reshape(len(data)/n,n)
        for i in range(data2.shape[0]):
            holder[i]=data2[i]
        return holder
        """
   
    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return


class AllDestroyedData(PluginBase):
    def __init__(self):
        super(AllDestroyedData,self).__init__()
        self.filename='AllDestroyedData.dat'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = True
        self.xlabel='' ; self.ylabel=''  # want these to be adjustable
        self.autofigname=''
        self.min_mass = 10**6.0

    def _analyze(self,hpath):
        DD = DestroyedDataFirstPass()
        dataD = DD.read(hpath)
        dtype = ['peak_mgrav','infall_mgrav','peak_hostid_RS','infall_hostid_RS','peak_rvmax','infall_rvmax','peak_corevelx','peak_corevely','peak_corevelz','infall_corevelx','infall_corevely','infall_corevelz','nstars','start_pos']
        data_newD = pandas.DataFrame(np.zeros((len(dataD),len(dtype)))-1,columns=dtype)
        peak_dataD = {}
        for peaksnap,line in zip(dataD['peak_snap'],dataD.index):
            peak_dataD.setdefault(peaksnap, []).append(line)

        infall_dataD = {}
        for infallsnap,line in zip(dataD['infall_snap'],dataD.index):
            infall_dataD.setdefault(infallsnap, []).append(line)       

        # initialize arrays
        allstars=[]; start_pos=0

        for snap in range(haloutils.get_numsnaps(hpath)):
            print snap, 'snap in get extra parameters Destroyed'
            sys.stdout.flush()
            if peak_dataD.has_key(snap) or infall_dataD.has_key(snap):
                cat=haloutils.load_rscat(hpath,snap,rmaxcut=False)

                if peak_dataD.has_key(snap):
                    for line in peak_dataD[snap]:
                        peak_rsid = int(dataD.ix[line]['peak_rsid'])
                        data_newD.ix[line]['peak_mgrav'] = cat.ix[peak_rsid]['mgrav']
                        data_newD.ix[line]['peak_hostid_RS'] = cat.ix[peak_rsid]['hostID']
                        data_newD.ix[line]['peak_rvmax'] = cat.ix[peak_rsid]['rvmax']
                        data_newD.ix[line]['peak_corevelx'] = cat.ix[peak_rsid]['corevelx']
                        data_newD.ix[line]['peak_corevely'] = cat.ix[peak_rsid]['corevely']
                        data_newD.ix[line]['peak_corevelz'] = cat.ix[peak_rsid]['corevelz']
                        
                if infall_dataD.has_key(snap):
                    for line in infall_dataD[snap]:
                        infall_rsid = int(dataD.ix[line]['infall_rsid'])
                        data_newD.ix[line]['infall_mgrav'] = cat.ix[infall_rsid]['mgrav']
                        data_newD.ix[line]['infall_hostid_RS'] = cat.ix[infall_rsid]['hostID']
                        data_newD.ix[line]['infall_rvmax'] = cat.ix[infall_rsid]['rvmax']
                        data_newD.ix[line]['infall_corevelx'] = cat.ix[infall_rsid]['corevelx']
                        data_newD.ix[line]['infall_corevely'] = cat.ix[infall_rsid]['corevely']
                        data_newD.ix[line]['infall_corevelz'] = cat.ix[infall_rsid]['corevelz']

                        iPids = cat.get_all_particles_from_halo(infall_rsid)
                        star_pids = iPids[0:int(np.round(len(iPids)*.03))]
                        data_newD.ix[line]['nstars'] = len(star_pids)
                        data_newD.ix[line]['start_pos'] = start_pos
                        allstars=np.r_[allstars,star_pids]
                        start_pos+=len(star_pids)
                
        fulldataD = pandas.concat((dataD,data_newD),axis=1)
        fulldataD.to_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t')
        f = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'destroyedPIDs.dat', 'wb')
        np.array(allstars).tofile(f)
        f.close()

    def _read(self,hpath):
        return pandas.read_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t')

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return



# Moster et al stellar mass to halo mass relationship
def getFraction(M, a):
    M10 = 11.590
    M11 = 1.195
    N10 = .0351
    N11 = -0.0247
    B10 = 1.376
    B11 = -0.826
    G10 = 0.608
    G11 = 0.329
    def M1(a):
        return 10**(M10+M11*(1-a))
    def N(a):
        return N10 + N11*(1-a)
    def beta(a):
        return B10 + B11*(1-a)
    def gamma(a):
        return G10 + G11*(1-a)
    return 2*N(a)*( (M/M1(a) )**-beta(a) + ( M/M1(a) )**gamma(a) )**-1
    

class TagMass(PluginBase):
    def __init__(self):
        super(TagMass,self).__init__()
        self.filename='destroyedMass_moster.dat'
        self.xmin=1;     self.xmax=400
        self.ymin=10**-5; self.ymax=10*6
        self.xlog= True; self.ylog = True
        self.xlabel='' ; self.ylabel=r''
        self.autofigname=''

    def _analyze(self,hpath):
        # RetagExtant
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
        Extant = AllExtantData()
        dataE = Extant.read(hpath)
        idsE = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'extantPIDs.dat')
        fracs = getFraction(dataE['infall_mvir']/cat.h0, haloutils.get_scale_snap(hpath, np.array(dataE['infall_snap'],dtype=np.int32)) )
        # need something better than getScale
        nstars = np.array(dataE['nstars'],dtype=np.int32)
        mper=(dataE['infall_mvir']/cat.h0*fracs)/nstars
        sp = np.array(dataE['start_pos'],dtype=np.int32)

        # must initialize mper_arr as same length as ids
        # then must index it as mper_arr[start_pos:start_pos+nstars]=mper[i]*nstars[i]. So that they properly align.
        mper_arr = np.zeros(len(idsE))
        for i in range(len(dataE)):
            mper_arr[sp[i]:sp[i]+nstars[i]] = [mper[i]]*nstars[i]
        # rewrite data properly here
        g=open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'extantMass_moster.dat','wb')
        mper_arr.tofile(g)
        g.close()

        # Now retag destroyed data
        Destroyed = AllDestroyedData()
        dataD = Destroyed.read(hpath)
        idsD = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'destroyedPIDs.dat')
        mper_arr=np.zeros(len(idsD))
        fracs = getFraction(dataD['infall_mvir']/cat.h0, haloutils.get_scale_snap(hpath, np.array(dataD['infall_snap'],dtype=np.int32)) )
        nstars = np.array(dataD['nstars'],dtype=np.int32)
        mper=(dataD['infall_mvir']/cat.h0*fracs)/nstars
        sp = np.array(dataD['start_pos'],dtype=np.int32)
        for i in range(len(dataD)):
            mper_arr[sp[i]:sp[i]+nstars[i]] = [mper[i]]*nstars[i]
        g=open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'destroyedMass_moster.dat','wb')
        mper_arr.tofile(g)
        g.close()

    def _read(self,hpath):
        # extant data
        idsE = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'extantPIDs.dat')
        massE = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'extantMass_moster.dat')
        # destroyed data
        idsD = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'destroyedPIDs.dat')
        massD = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'destroyedMass_moster.dat')
        return np.array(idsE,dtype=np.int64), massE, np.array(idsD,dtype=np.int64), massD

"""
class TagMass_1(PluginBase):
    def __init__(self):
        super(TagMass_1,self).__init__()
        self.filename='destroyedMass_moster_1.dat'
        self.xmin=1;     self.xmax=400
        self.ymin=10**-5; self.ymax=10*6
        self.xlog= True; self.ylog = True
        self.xlabel='' ; self.ylabel=r''
        self.autofigname=''

    def _analyze(self,hpath):
        # RetagExtant
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
        Extant = AllExtantData()
        dataE = Extant.read(hpath)
        idsE = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'extantPIDs.dat')
        # re-write idsE
        idsE_1 = []; mper_arr = []
        for row in xrange(len(dataE)):
            idsE_1 = np.r_[idsE_1, getStars_1(dataE,idsE,row)]
            mper_arr = get
        g=open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'extantPIDs_1.dat','wb')
        np.array(idsE_1).tofile(g)
        g.close()
        # done rewriting idsE

        fracs = getFraction(dataE['infall_mvir']/cat.h0, haloutils.get_scale_snap(hpath, np.array(dataE['infall_snap'],dtype=np.int32)) )
        nstars = np.array(np.round(dataE['nstars']/3.),dtype=np.int32)
        mask0 = np.array(dataE['nstars'])==0
        mask = nstars==0
        nstars[mask]=1
        nstars[mask0]=0
        mper=(dataE['infall_mvir']/cat.h0*fracs)/nstars
        sp = np.array(dataE['start_pos'],dtype=np.int32) # THIS IS WRONG

        # must initialize mper_arr as same length as ids
        mper_arr = np.zeros(len(idsE_1))
        for i in xrange(len(dataE)):
            mper_arr[sp[i]:sp[i]+nstars[i]] = [mper[i]]*nstars[i]
        # rewrite data properly here
        g=open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'extantMass_moster_1.dat','wb')
        mper_arr.tofile(g)
        g.close()

        # Now retag destroyed data
        Destroyed = AllDestroyedData()
        dataD = Destroyed.read(hpath)
        idsD = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'destroyedPIDs.dat')
        # re-write idsD
        idsD_1 = []
        for row in xrange(len(dataD)):
            idsD_1 = np.r_[idsD_1, getStars_1(dataD,idsD,row)]
        g=open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'destroyedPIDs_1.dat','wb')
        np.array(idsD_1).tofile(g)
        g.close()
        # done rewriting idsD
        mper_arr=np.zeros(len(idsD))
        fracs = getFraction(dataD['infall_mvir']/cat.h0, haloutils.get_scale_snap(hpath, np.array(dataD['infall_snap'],dtype=np.int32)) )
        nstars = np.array(np.round(dataD['nstars']/3.),dtype=np.int32)
        mask0 = np.array(dataD['nstars'])==0
        mask = nstars==0
        nstars[mask]=1
        nstars[mask0]=0

        mper=(dataD['infall_mvir']/cat.h0*fracs)/nstars
        sp = np.array(dataD['start_pos'],dtype=np.int32) # THIS IS WRONG
        for i in xrange(len(dataD)):
            mper_arr[sp[i]:sp[i]+nstars[i]] = [mper[i]]*nstars[i]
        g=open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'destroyedMass_moster_1.dat','wb')
        mper_arr.tofile(g)
        g.close()

    def _read(self,hpath):
        # extant data
        idsE = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'extantPIDs_1.dat')
        massE = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'extantMass_moster_1.dat')
        # destroyed data
        idsD = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'destroyedPIDs_1.dat')
        massD = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'destroyedMass_moster_1.dat')
        return np.array(idsE,dtype=np.int64), massE, np.array(idsD,dtype=np.int64), massD
"""
