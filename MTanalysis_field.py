import numpy as np
import readsnapshots.readsnapHDF5_greg as rsg
import haloutils
import pandas
import time
from caterpillaranalysis import *
import DwarfMethods as dm

# Load all field halos outside host rvir, within contamination radius
# get all subhalos in those field halos, and record peak mass, peak time, infall mass, infall time


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
    dr = dm.distance(part_pos,halo_center,cat.boxsize)*cat.scale/cat.h0*1000  # in kpc
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
    mvirNFW = 4./3 * np.pi * rvirNFW**3 * rho_crit * delta
    #print rvir, rvirNFW, 'rvir350 and rvir350 NFW'
    return mvir[0]*cat.h0, mvirNFW[0]*cat.h0
    

def get_field_halos(cat,hpath,hostID,mlow=10,mhigh=11.5):
    # get halos beyond host halo virial radius, but less than
    # contamination radius                     
    hosts = cat.get_hosts()
    dists = dm.distance(cat.ix[hostID][['posX','posY','posZ']], hosts[['posX','posY','posZ']])
    contam_dist = dm.get_contam_dist(hpath)
    mask = (dists < contam_dist)&(dists > .001)
    mass_mask = (10**mlow < hosts['mgrav']/cat.h0) &(10**mhigh > hosts['mgrav']/cat.h0)
    return hosts[mask & mass_mask]

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




# want to write out z=0 mgrav of the subhalo
# want to compute and write out mgrav, m350, and m350NFW of the z=0 field halo hosts

# sub_mb['origid'][0] is the subhalo z=0 rsid
# mto.fieldid is the z=0 rsid of the host halo


#### will need to add in this function to get vmax at z=11,10,9,8 etc.
### also get peak vmax

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
    max_mass_scale = sub_mb[max_mass_loc]['scale']
    max_mass_rsid = sub_mb[max_mass_loc]['origid']
    max_mass_mvir = sub_mb[max_mass_loc]['mvir']

    # Get infall parameters 
    if iLoc !=None:
        infall_snap = sub_mb[iLoc]['snap']
        infall_scale = sub_mb[iLoc]['scale']
        infall_rsid = sub_mb[iLoc]['origid']
        infall_vmax = sub_mb[iLoc]['vmax']
        infall_mvir = sub_mb[iLoc]['mvir']
    else:
        infall_snap, infall_scale, infall_rsid, infall_vmax, infall_mvir = [-1]*5
    mto.otherdata=np.r_[mto.otherdata, mto.fieldid,subrank,sub_mb['origid'][0], sub_mb['snap'][0],depth,max_mass_rsid, max_mass_snap, max_mass_scale, max_mass_vmax,max_mass_mvir,infall_rsid,infall_snap,infall_scale,infall_vmax,infall_mvir, sub_mb['mvir'][0], mto.cat.ix[mto.fieldid]['mgrav'], mto.m350_host, mto.m350NFW_host] 
    return


class MT_Object(object):
    def __init__(self,hpath):
        print 'creating MT object'
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        self.cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
        self.hostID = haloutils.load_zoomid(hpath)
        self.field_halos = get_field_halos(self.cat,hpath,self.hostID,mlow=10,mhigh=12)
        self.field_rsids = np.array(self.field_halos['id'])
        self.hosthalo = self.cat.ix[self.hostID]
        self.mtc = haloutils.load_mtc(hpath,haloids=self.field_rsids,indexbyrsid=True)
        print 'loaded mtc'
        self.min_mass = 10**7.0
        self.otherdata = []
        self.start_time = time.time()

        # not sure about these yet
        #self.hosttree 
        #self.searchtree
        #self.mmp_map = 
        #self.non_mmp_map = 



class HostHaloM350(PluginBase):
    def __init__(self):
        super(HostHaloM350,self).__init__()
        self.filename='HostM350.dat'
    def _analyze(self,hpath):
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID = haloutils.load_zoomid(hpath)                
        pids = cat.get_all_particles_from_halo(hostID)
        m350_host, m350NFW_host = getM_xcrit(hpath,pids,cat,hostID,delta=350)

        g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,'wb')
        np.array([m350_host, m350NFW_host]).tofile(g)
        g.close()


    def _read(self,hpath):
        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename)
        pdtype = ['m350_host','m350NFW_host']
        n = len(pdtype)
        import pandas
        return pandas.DataFrame(data.reshape(len(data)/n,n), columns=pdtype)



# need to add peak vmax values.
# need to add vmax at the time of         


# start by just getting peak mass.
# infall mass harder to compute
class FieldHaloSubstructureFirstPass(PluginBase):
    def __init__(self):
        super(FieldHaloSubstructureFirstPass,self).__init__()
        self.filename='FieldHaloSubstructure.dat'

    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
                raise IOError("No rockstar")
        mto = MT_Object(hpath)
        
        for field_rsid, num in zip(mto.field_rsids, np.arange(1,len(mto.field_rsids)+1)):
            mto.fieldid = field_rsid

            ### slows things down, may not need in future versions
            pids = mto.cat.get_all_particles_from_halo(field_rsid)
            mto.m350_host, mto.m350NFW_host = getM_xcrit(hpath,pids,mto.cat,field_rsid,delta=350)
            print  mto.m350_host, mto.m350NFW_host, mto.cat.ix[mto.fieldid]['mgrav'], 'm350 and m350 NFW and mgrav'
            ###
            subs = mto.cat.get_subhalos_within_halo(field_rsid)
            hosttree = mto.mtc.Trees[field_rsid] # host is the field halo
            host_mmp_map = hosttree.get_mmp_map()
            host_mb = hosttree.getMainBranch(0, host_mmp_map)
            good=0; toosmall=0
            for subRSID, subrank in zip(np.array(subs['id']),np.arange(1,len(subs)+1)):
                if mto.mtc.Trees.has_key(subRSID):
                    mto.searchtree = mto.mtc.Trees[subRSID]
                else:
                    if np.log10(mto.cat.ix[subRSID]['mgrav']/mto.cat.h0) > 7.0:
                        print 'subhalo',subRSID,'not in merger tree. mass = ', np.log10(mto.cat.ix[subRSID]['mgrav']/mto.cat.h0)
                    continue
                mto.mmp_map = mto.searchtree.get_mmp_map()
                mto.non_mmp_map = mto.searchtree.get_non_mmp_map()
                sub_mb = mto.searchtree.getMainBranch(0,mto.mmp_map)
                if sub_mb is None:
                    print 'subhalo', sub_rank, 'main branch not found in MT. Skipping it. Z=0 Mass: %.4e, Vmax: %.4f' %(mto.cat.ix[subRSID]['mgrav'], mto.cat.ix[subRSID]['vmax'])
                    sys.stdout.flush()
                    continue # skip to next subhalo 
                max_mass = np.max(sub_mb['mvir'])
                if max_mass/mto.cat.h0 < mto.min_mass:
                    toosmall+=1
                    continue
                iLoc, iSnap = getInfall(sub_mb, host_mb, max_mass)
                add_data(mto,sub_mb,iLoc,subrank,depth=0)
                if iLoc is None:
                    iLoc=len(sub_mb)-1 # only for deciding new end_level
                auxiliary_add(mto,host_mb[1:],cur_line=0,level=0, end_level=iLoc, subrank=-abs(subrank), depth=1)
                if subrank%50==0 or subrank==1:
                    print subrank, '/', len(subs), 'finished. Time = ', (time.time()-mto.start_time)/60., 'minutes'
                sys.stdout.flush()
                good+=1
            print 'Done with Field Halo', num,'/',len(mto.field_rsids), 'Time = ', (time.time()-mto.start_time)/60., 'minutes'
            print good, 'halos good out of', len(subs)
            print toosmall, 'num halos too small'
            print len(subs)-good-toosmall, 'number of subhalo failures'
            sys.stdout.flush()
                
        g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,'wb')
        np.array(mto.otherdata).tofile(g)
        g.close()
        print 'All finished. Time = ', (time.time()-mto.start_time)/60., 'minutes'

    def _read(self,hpath):
        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename)
        pdtype = ['field_rsid','sub_rank','rsid','backsnap','depth','max_mass_rsid','max_mass_snap','max_mass_scale','max_mass_vmax','max_mass','infall_rsid','infall_snap','infall_scale','infall_vmax','infall_mvir', 'mgravz0', 'host_mgravz0','host_m350z0','host_m350NFWz0']
        n = len(pdtype)
        import pandas
        return pandas.DataFrame(data.reshape(len(data)/n,n), columns=pdtype)


class AllExtantFieldData(PluginBase):
    def __init__(self):
        super(AllExtantFieldData,self).__init__()
        self.filename='AllExtantFieldData.dat'

    def _analyze(self,hpath):
        #frac_to_tag = .05
        ED = FieldHaloSubstructureFirstPass()
        dataE = ED.read(hpath)
        dtype = ['max_mass_mgrav','infall_mgrav','max_mass_hostid_RS','infall_hostid_RS','max_mass200','infall_mass200','max_mass350', 'max_mass350NFW']
        data_newE = pandas.DataFrame(np.zeros((len(dataE),len(dtype)))-1,columns=dtype)
        maxmass_dataE = {}
        for maxmass_snap,line in zip(dataE['max_mass_snap'],dataE.index):
            maxmass_dataE.setdefault(maxmass_snap, []).append(line)

        infall_dataE = {}
        for infallsnap,line in zip(dataE['infall_snap'],dataE.index):
            infall_dataE.setdefault(infallsnap, []).append(line)

        # initialize arrays for tagging                                                                                           
        allstars=[]; start_pos=0

        for snap in range(haloutils.get_numsnaps(hpath)):
            sys.stdout.flush()
            if maxmass_dataE.has_key(snap) or infall_dataE.has_key(snap):
                cat=haloutils.load_rscat(hpath,snap,rmaxcut=False)
                print snap, 'snap in get extra parameters Extant'
                if maxmass_dataE.has_key(snap):
                    for line in maxmass_dataE[snap]:
                        maxmass_rsid = int(dataE.ix[line]['max_mass_rsid'])
                        data_newE.ix[line]['max_mass_mgrav'] = cat.ix[maxmass_rsid]['mgrav']
                        data_newE.ix[line]['max_mass_hostid_RS'] = cat.ix[maxmass_rsid]['hostID']
                        data_newE.ix[line]['max_mass200'] = cat.ix[maxmass_rsid]['altm2']
                        
                        pids = cat.get_all_particles_from_halo(maxmass_rsid)
                        m350,m350NFW = getM_xcrit(hpath,pids,cat,maxmass_rsid,delta=350)
                        data_newE.ix[line]['max_mass350'] = m350
                        data_newE.ix[line]['max_mass350NFW'] = m350NFW
                        #print cat.ix[maxmass_rsid]['mgrav']/cat.h0, 'max mass mgrav'
                        #print m350/cat.h0, 'max mass m350'
                        #print m350NFW/cat.h0, 'max mass m350 NFW'
                        #print cat.ix[maxmass_rsid]['altm2']/cat.h0, 'max mass 200'

                if infall_dataE.has_key(snap):
                    for line in infall_dataE[snap]:
                        infall_rsid = int(dataE.ix[line]['infall_rsid'])
                        data_newE.ix[line]['infall_mgrav'] = cat.ix[infall_rsid]['mgrav']
                        data_newE.ix[line]['infall_hostid_RS'] = cat.ix[infall_rsid]['hostID']
                        data_newE.ix[line]['infall_mass200'] = cat.ix[infall_rsid]['altm2']

                        #iPids = cat.get_all_particles_from_halo(infall_rsid)
                        #star_pids = iPids[0:int(np.round(len(iPids)*frac_to_tag))]
                        #data_newE.ix[line]['nstars'] = len(star_pids)
                        #data_newE.ix[line]['start_pos'] = start_pos
                        #allstars=np.r_[allstars,star_pids]
                        #start_pos+=len(star_pids)

        fulldataE = pandas.concat((dataE,data_newE),axis=1)
        fulldataE.to_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t')
        #f = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'extantPIDs.dat', 'wb')
        #np.array(allstars).tofile(f)
        #f.close()

    def _read(self,hpath):
        return pandas.read_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t')

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return
