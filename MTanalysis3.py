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

# make it possible to get any fraction of 5%
h0 = .6711

def load_idsE(hpath):
    return np.fromfile(hpath+'/'+'analysis'+'/'+'extantPIDs.dat') 

def load_idsD(hpath):
    return np.fromfile(hpath+'/'+'analysis'+'/'+'destroyedPIDs.dat') 

def mass_per_particle(hpath,data, row):
    frac = getFraction(data[row:row+1]['infall_mvir']/h0, haloutils.get_scale_snap(hpath, int(data[row:row+1]['infall_snap'])))
    nstars = int(data['nstars'][row:row+1])
    return float((data['infall_mvir'][row:row+1]/h0*frac)/nstars)


def getStars(data, ids, row):
    sp = data['start_pos'][row]
    nstars = data['nstars'][row]
    return np.array(ids[sp:sp+nstars],dtype=np.int64)

def getStars_x(data,ids,row, fmb=1):
    if fmb  > 5:
        print "ERROR fmb cannot be > 5"
        return None
    sp = data['start_pos'][row]
    nstars = np.round(data['nstars'][row]/(5./fmb))
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


def Hubble(cat):
    to_inv_sec = 3.241*10**-20
    a = cat.scale
    Ok = 1-cat.Om-cat.Ol
    return to_inv_sec * 100*cat.h0*(cat.Om*(1/a)**3 + Ok*(1/a)**2 + cat.Ol)**.5 # value in 1/sec

def m_enclNFW(r,rs,rho_0):
    return rho_0*4*np.pi*rs**3 *( np.log((rs+r)/rs) - r/(r+rs) )

def rho_enlcNWF(r,rs,rho_0):
    mencl = rho_0*4*np.pi*rs**3 *( np.log((rs+r)/rs) - r/(r+rs) )
    return mencl / (4/3. * np.pi * r**3)

# this will likely be a slight underestimate of what rockstar would determine
# this is due to not including the unbound particles
def getM_xcrit(hpath,pids,cat,rsid,delta=350):
    # delta*rho_crit
    G = 4.157e-39  # kpc^3/Msun/s^2
    H = Hubble(cat) # in 1/s
    rho_crit = (3*H**2)/(8*np.pi*G) # in Msun/kpc^3

    # get mltr interpolating function
    pids = np.sort(pids)
    halo = cat.ix[rsid]
    halo_center = np.array(halo[['posX','posY','posZ']])
    part_pos = haloutils.load_partblock(hpath,cat.snap_num,'POS ',parttype=1,ids=pids)
    dr = distance(part_pos,halo_center,cat.boxsize)*cat.scale/cat.h0*1000  # in kpc
    maxr = 1.1*float(halo['rvir'])*cat.scale/cat.h0
    minr = 0.1
    binwidth = 0.04
    nbins = np.ceil((np.log10(maxr)-np.log10(minr))/binwidth)
    rarr = 10**np.linspace(np.log10(minr), np.log10(minr)+nbins*binwidth,nbins+1)
    h_r, x_r = np.histogram(dr, bins=np.concatenate(([0],rarr)))
    m_lt_r = np.cumsum(h_r)*cat.particle_mass/cat.h0
    rho = m_lt_r / (4/3. * np.pi * rarr**3)
    tck_mltr = interpolate.splrep(rarr,m_lt_r)
    tck_rho = interpolate.splrep(rarr,rho)

    # also solve assuming NFW with r_scale
    rs_NFW = halo['rs']*cat.scale/cat.h0 # in kpc
    m200 = halo['altm2']/cat.h0 # in Msun
    r200 = (m200*3/(4*np.pi*200*rho_crit))**(1/3.) # in kpc
    rho_0 = m200/ (4*np.pi*rs_NFW**3 *( np.log((rs_NFW+r200)/rs_NFW) - r200/(r200+rs_NFW) )) # msun/kpc^3
    
    from scipy.optimize import fsolve
    def func(x):
        return interpolate.splev(x,tck_rho)/(rho_crit) - delta
    rvir = fsolve(func,1)
    mvir = interpolate.splev(rvir,tck_mltr)
    
    def funcNFW(x):
        return rho_enlcNWF(x,rs_NFW,rho_0)/(rho_crit) - delta
    rvirNFW = fsolve(funcNFW,1)
    mvirNFW = 4./3 * np.pi * rvirNFW**3 * rho_crit * 350
    #print rvir, rvirNFW, 'rvir350 and rvir350 NFW'
    return mvir[0]*cat.h0, mvirNFW[0]*cat.h0
    







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
class ExtantDataFirstPass(PluginBase):
    def __init__(self):
        super(ExtantDataFirstPass,self).__init__()
        self.filename='ExtantDataFirstPass.dat'

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
            add_data(mto,sub_mb,iLoc,subrank,depth=0)
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

    def _read(self,hpath):
        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename)
        pdtype = ['sub_rank','rsid','pid', 'backsnap','depth','max_mass_rsid','max_mass_snap','max_mass_vmax','max_mass','max_mass_posx','max_mass_posy','max_mass_posz','max_mass_pecvx','max_mass_pecvy','max_mass_pecvz','max_mass_virialratio','max_mass_hostid_MT','max_mass_rvir','max_mass_spinbullock','max_mass_rs','max_mass_scale_of_last_MM','max_mass_Jx','max_mass_Jy','max_mass_Jz','max_mass_xoff', 'peak_rsid','peak_snap','peak_vmax','peak_mvir','peak_posx','peak_posy','peak_posz','peak_pecvx','peak_pecvy','peak_pecvz','peak_virialratio','peak_hostid_MT','peak_rvir','peak_spinbullock','peak_rs','peak_scale_of_last_MM','peak_Jx','peak_Jy','peak_Jz','peak_xoff','infall_rsid','infall_snap','infall_vmax','infall_mvir','infall_posx','infall_posy','infall_posz','infall_pecvx','infall_pecvy','infall_pecvz','infall_virialratio','infall_hostid_MT','infall_rvir','infall_spinbullock','infall_rs','infall_scale_of_last_MM','infall_Jx','infall_Jy','infall_Jz','infall_xoff']
        n = len(pdtype)
        import pandas
        return pandas.DataFrame(data.reshape(len(data)/n,n), columns=pdtype)

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return


class AllExtantData(PluginBase):
    def __init__(self):
        super(AllExtantData,self).__init__()
        #self.filename='AllExtantData.dat'
        self.filename='ExtantDataSecondPass.dat'

    # tag 5% of particles
    def _analyze(self,hpath):
        frac_to_tag = .05
        ED = ExtantDataFirstPass()
        dataE = ED.read(hpath)
        # adding alternate mass definitions at the end
        dtype = ['peak_mgrav','infall_mgrav','peak_hostid_RS','infall_hostid_RS','peak_rvmax','infall_rvmax','peak_corevelx','peak_corevely','peak_corevelz','infall_corevelx','infall_corevely','infall_corevelz', 'nstars', 'start_pos', 'max_mass350','max_mass350NFW','infall_mass200','max_mass200']
        data_newE = pandas.DataFrame(np.zeros((len(dataE),len(dtype)))-1,columns=dtype)
        peak_dataE = {}
        for peaksnap,line in zip(dataE['peak_snap'],dataE.index):
            peak_dataE.setdefault(peaksnap, []).append(line)

        infall_dataE = {}
        for infallsnap,line in zip(dataE['infall_snap'],dataE.index):
            infall_dataE.setdefault(infallsnap, []).append(line)       

        maxmass_dataE = {}
        for maxmass_snap,line in zip(dataE['max_mass_snap'],dataE.index):
            maxmass_dataE.setdefault(maxmass_snap, []).append(line)       


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

                if maxmass_dataE.has_key(snap):
                    for line in maxmass_dataE[snap]:
                        maxmass_rsid = int(dataE.ix[line]['max_mass_rsid'])
                        data_newE.ix[line]['max_mass200'] = cat.ix[maxmass_rsid]['altm2']
                        pids = cat.get_all_particles_from_halo(maxmass_rsid)
                        m350,m350NFW = getM_xcrit(hpath,pids,cat,maxmass_rsid,delta=350)
                        data_newE.ix[line]['max_mass350'] = m350
                        data_newE.ix[line]['max_mass350NFW'] = m350NFW

                if infall_dataE.has_key(snap):
                    for line in infall_dataE[snap]:
                        infall_rsid = int(dataE.ix[line]['infall_rsid'])
                        data_newE.ix[line]['infall_mgrav'] = cat.ix[infall_rsid]['mgrav']
                        data_newE.ix[line]['infall_hostid_RS'] = cat.ix[infall_rsid]['hostID']
                        data_newE.ix[line]['infall_rvmax'] = cat.ix[infall_rsid]['rvmax']
                        data_newE.ix[line]['infall_corevelx'] = cat.ix[infall_rsid]['corevelx']
                        data_newE.ix[line]['infall_corevely'] = cat.ix[infall_rsid]['corevely']
                        data_newE.ix[line]['infall_corevelz'] = cat.ix[infall_rsid]['corevelz']
                        
                        data_newE.ix[line]['infall_mass200'] = cat.ix[infall_rsid]['altm2']

                        iPids = cat.get_all_particles_from_halo(infall_rsid)
                        star_pids = iPids[0:int(np.round(len(iPids)*frac_to_tag))]
                        data_newE.ix[line]['nstars'] = len(star_pids)
                        data_newE.ix[line]['start_pos'] = start_pos
                        allstars=np.r_[allstars,star_pids]
                        start_pos+=len(star_pids)
                
        #fulldataE = pandas.concat((dataE,data_newE),axis=1)
        #fulldataE.to_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t')
        data_newE.to_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t')
        f = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'extantPIDs.dat', 'wb')
        np.array(allstars).tofile(f)
        f.close()

    def _read(self,hpath):
        ED = ExtantDataFirstPass()
        dataE = ED.read(hpath)
        data_newE = pandas.read_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t', index_col=0)
        fulldataE = pandas.concat((dataE,data_newE),axis=1)
        return fulldataE        
        #return pandas.read_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t')

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return


# only valid with AllExtantData in original form
class extract_dataE(PluginBase):
    def __init__(self):
        super(extract_dataE,self).__init__()
        self.filename='ExtantDataSecondPass.dat'

    def _analyze(self,hpath):
        # load up the full catalog, split it off, and re-write just the new data
        AE = AllExtantData()
        data = AE.read(hpath)
        data_newE = data.ix[:,-18:]
        data_newE.to_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t')

    def _read(self,hpath):
        ED = ExtantDataFirstPass()
        dataE = ED.read(hpath)
        data_newE = pandas.read_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t', index_col=0)
        print data_newE.columns
        fulldataE = pandas.concat((dataE,data_newE),axis=1)
        return fulldataE        
        #return pandas.read_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t')

        # In AllExtantData, write out EDSP to file, and combine first pass and second pass in the read function



#### will need to add in this function to get vmax at z=11,10,9,8 etc.
## can it be done without having to re-run AllExtantData??
# I could copy the code, and write out an entirely new dataset with just the info needed for reionization
# make it the same order as the original data set so that they can be aligned for more info
# reionization values needed: vmax at all reionization redshifts. vmax at peak, snap at peak. mass at peak.
# mass at infall, m200 at infall, m350 at peak. those are slow, don't want to re-compute, so I have to load
# old catalog anyway. make a new catalog with just the 4 new values
# then do the pandas.concat((dataE,data_reion),axis=1) in dwarf methods, get extant data.
# make one overlapping column to be sure its correct. like depth, sub_rank, rsid. but name them
# depth_reion, sub_rank_reion, rsid_reion.
# first pass should be fast. wake up early have it run before I meet Anna.

def add_data(mto,sub_mb,iLoc,subrank,depth):
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
    if iLoc !=None:
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
    else:
        infall_snap, infall_rsid, infall_vmax, infall_mvir, infall_posx, infall_posy = [-1]*6
        infall_posz, infall_pecvx, infall_pecvy, infall_pecvz, infall_virialratio = [-1]*5
        infall_hostid_MT, infall_rvir, infall_spinbullock, infall_rs = [-1]*4
        infall_scale_of_last_MM, infall_Jx, infall_Jy, infall_Jz, infall_xoff = [-1]*5
        
    if depth ==0:
        RSID = int(sub_mb['origid'][0])
        PID = int(mto.cat.ix[RSID]['hostID'])
    else:
        PID = -100

    mto.otherdata=np.r_[mto.otherdata,subrank,sub_mb['origid'][0], PID, sub_mb['snap'][0],depth,max_mass_rsid, max_mass_snap, max_mass_vmax,max_mass_mvir,max_mass_posx,max_mass_posy,max_mass_posz,max_mass_pecvx,max_mass_pecvy,max_mass_pecvz,max_mass_virialratio,max_mass_hostid_MT,max_mass_rvir,max_mass_spinbullock,max_mass_rs,max_mass_scale_of_last_MM,max_mass_Jx,max_mass_Jy,max_mass_Jz,max_mass_xoff, peak_rsid, peak_snap, peak_vmax,peak_mvir,peak_posx,peak_posy,peak_posz,peak_pecvx,peak_pecvy,peak_pecvz,peak_virialratio,peak_hostid_MT,peak_rvir,peak_spinbullock,peak_rs,peak_scale_of_last_MM,peak_Jx,peak_Jy,peak_Jz,peak_xoff,infall_rsid,infall_snap,infall_vmax,infall_mvir,infall_posx,infall_posy,infall_posz,infall_pecvx,infall_pecvy,infall_pecvz,infall_virialratio,infall_hostid_MT,infall_rvir,infall_spinbullock,infall_rs,infall_scale_of_last_MM,infall_Jx,infall_Jy,infall_Jz,infall_xoff]
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
            add_data(mto,sub_mb, iLoc,subrank,depth)
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


class DestroyedDataFirstPass(PluginBase):
    def __init__(self):
        super(DestroyedDataFirstPass,self).__init__()
        self.filename='DestroyedDataFirstPass.dat'
        self.filestring='DestroyedDataFirstPass'
       
    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
        mto = MT_Object(hpath)
        host_mb = mto.hosttree.getMainBranch(0, mto.mmp_map)
        print len(host_mb), 'lengh of host main branch'
        auxiliary_add(mto,host_mb[1:],cur_line=0,level=0,end_level=len(host_mb), subrank=0.5,depth=1, destr=True)
        print 'about to write data to file'
        g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,'wb')
        print 'opened file'
        np.array(mto.otherdata).tofile(g)
        print 'wrote data, file still open'
        g.close()    
        print 'Finished Destroyed First Pass'
        
    def _read(self,hpath):
        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename)
        pdtype = ['sub_rank','rsid','backsnap','depth','max_mass_rsid','max_mass_snap','max_mass_vmax','max_mass','max_mass_posx','max_mass_posy','max_mass_posz','max_mass_pecvx','max_mass_pecvy','max_mass_pecvz','max_mass_virialratio','max_mass_hostid_MT','max_mass_rvir','max_mass_spinbullock','max_mass_rs','max_mass_scale_of_last_MM','max_mass_Jx','max_mass_Jy','max_mass_Jz','max_mass_xoff', 'peak_rsid','peak_snap','peak_vmax','peak_mvir','peak_posx','peak_posy','peak_posz','peak_pecvx','peak_pecvy','peak_pecvz','peak_virialratio','peak_hostid_MT','peak_rvir','peak_spinbullock','peak_rs','peak_scale_of_last_MM','peak_Jx','peak_Jy','peak_Jz','peak_xoff','infall_rsid','infall_snap','infall_vmax','infall_mvir','infall_posx','infall_posy','infall_posz','infall_pecvx','infall_pecvy','infall_pecvz','infall_virialratio','infall_hostid_MT','infall_rvir','infall_spinbullock','infall_rs','infall_scale_of_last_MM','infall_Jx','infall_Jy','infall_Jz','infall_xoff']
        n = len(pdtype)
        import pandas
        return pandas.DataFrame(data.reshape(len(data)/n,n), columns=pdtype)
       
    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return


class AllDestroyedData(PluginBase):
    def __init__(self):
        super(AllDestroyedData,self).__init__()
        self.filename='AllDestroyedData.dat'

    def _analyze(self,hpath):
        frac_to_tag = .05
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
                        star_pids = iPids[0:int(np.round(len(iPids)*frac_to_tag))]
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
