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
 

Tvir = 2000
append_time = 0.0
main_branch_time = 0.0
mini_halo_time = 0.0
 
 
def mcrit(T,z):
    h = 0.6711
    omega_m = 0.3125
    M = 1e8/h * (T/(1.+z))**1.5 * (0.6*10/1.22/1.98e4)**1.5 * (18*3.14*3.14/178/omega_m)**0.5 #in solar masses
    return M
 
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
        self.filename='ExtantMiniHalosFirstPass.dat'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = True
        self.xlabel='scale factor' ; self.ylabel='Mass Accreted'    # want these to be adjustable
        self.autofigname='MergerHistory'
        self.min_mass = 10**5.5
 
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
         
        #sys.stdout.flush()
        mtc = haloutils.load_mtc(hpath,haloids=[hostID])
        print 'loaded mtc'
        #sys.stdout.flush()
        host = mtc.Trees[0]
        host_mb = host.getMainBranch(0)
         
        good = 0; too_small=0; sub_rank=start-1; start_big = 0
        for subRSID in np.array(subs['id']):
            sub_rank+=1; flag=0
            sub = getSubTree(mtc,subRSID)
            mmp_map = sub.get_mmp_map()
            non_mmp_map = sub.get_non_mmp_map()

            if sub==None:
                #print sub_rank, 'subhalo not found in MTCatalogue. Z=0 Mass: %.4e, Vmax: %.4f' %(cat.ix[subRSID]['mgrav'], cat.ix[subRSID]['vmax']), 'Time = ', (time.time()-start_time)/60., 'minutes'
                #sys.stdout.flush()
                continue
            sub_mb = sub.getMainBranch(0,mmp_map)
            if sub_mb == None:
                #print 'subhalo', sub_rank, 'main branch not found in MT. Skipping it. Z=0 Mass: %.4e, Vmax: %.4f' %(cat.ix[subRSID]['mgrav'], cat.ix[subRSID]['vmax'])
                #sys.stdout.flush()
                continue # skip to next subhalo

            # minihalo_loc is index of sub_mb where this first happens
            sub_redshifts = (1./sub_mb['scale']) - 1
            mcut = mcrit(Tvir,sub_redshifts)
            # check for phantoms
            index = np.where(np.logical_and(sub_mb['mvir']/cat.h0 > mcut, sub_mb['phantom']==0))[0]
            if len(index)==0:
                 minihalo_loc = None
                 too_small+=1
            else:
                 minihalo_loc = index[0]
            
            # special case: if main branch starts above 10**6. what to do?
            mcut_last = mcrit(Tvir,sub_redshifts[-1])
            if sub_mb['mvir'][-1]/cat.h0 > mcut_last:
                start_big += 1
                flag = 10

            otherdata, start_big = auxiliary_add(cat,host_mb,otherdata,host=sub,subline=0,ii=0, j=sub_rank, end=len(sub_mb)-1,min_mass=self.min_mass,start_big=start_big,flag=-1,mmp_map=mmp_map, non_mmp_map=non_mmp_map)
            if minihalo_loc == None:
                continue
            otherdata = add_data(otherdata,sub_mb, minihalo_loc,sub_rank,flag)

            if sub_rank%100==0:
                print sub_rank, '/', len(subs), 'finished. Time = ', (time.time()-start_time)/60., 'minutes'
            sys.stdout.flush()
            good+=1
        g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,'wb')
        np.array(otherdata).tofile(g)
        g.close()
# to test, use haloutils.get_hpath_lx(hid,lx)
        print good, 'halos good out of', len(subs)
        print too_small,'/',len(subs),'always too small'
        print start_big, "number of halos above the threshold to begin with"
        sys.stdout.flush()
 
    def _read(self,hpath):
        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename)
        pdtype = ['sub_rank','rsid','minihalo_rsid','minihalo_snap','minihalo_vmax','max_mass','minihalo_posx','minihalo_posy','minihalo_posz','minihalo_pecvx','minihalo_pecvy','minihalo_pecvz','minihalo_virialratio','minihalo_hostid_MT','minihalo_rvir','minihalo_spinbullock','minihalo_rs','minihalo_scale_of_last_MM','minihalo_Jx','minihalo_Jy','minihalo_Jz','minihalo_xoff','root_snap','flags']
        n = len(pdtype)
        import pandas
        return pandas.DataFrame(data.reshape(len(data)/n,n), columns=pdtype)
 
    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return


class AllExtantData(PluginBase):
    def __init__(self):
        super(AllExtantData,self).__init__()
        self.filename='MiniHaloAllExtantData.dat'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = True
        self.xlabel='' ; self.ylabel=''  # want these to be adjustable
        self.autofigname=''
        self.min_mass = 10**5.5
 
    def _analyze(self,hpath):
        ED = ExtantDataFirstPass()
        dataE = ED.read(hpath)
        # any values in RSCatalog you want that are not in MT
        dtype = ['minihalo_mgrav','minihalo_hostid_RS','minihalo_rvmax','minihalo_corevelx','minihalo_corevely','minihalo_corevelz', 'nstars', 'start_pos'] # nstars and start_pos should stay
        data_newE = pandas.DataFrame(np.zeros((len(dataE),len(dtype)))-1,columns=dtype)
 
        minihalo_dataE = {}
        for infallsnap,line in zip(dataE['minihalo_snap'],dataE.index):
            minihalo_dataE.setdefault(infallsnap, []).append(line)       
 
        # initialize arrays for tagging
        allstars=[]; start_pos=0    
        
        for snap in range(haloutils.get_numsnaps(hpath)):
            print snap, 'snap in get extra parameters Extant'
            sys.stdout.flush()
            if minihalo_dataE.has_key(snap):
                cat=haloutils.load_rscat(hpath,snap,rmaxcut=False)
                         
                if minihalo_dataE.has_key(snap):
                    for line in minihalo_dataE[snap]:
                        # extract data from rs catalog here
                        minihalo_rsid = int(dataE.ix[line]['minihalo_rsid'])
                        data_newE.ix[line]['minihalo_mgrav'] = cat.ix[minihalo_rsid]['mgrav']
                        data_newE.ix[line]['minihalo_hostid_RS'] = cat.ix[minihalo_rsid]['hostID']
                        data_newE.ix[line]['minihalo_rvmax'] = cat.ix[minihalo_rsid]['rvmax']
                        data_newE.ix[line]['minihalo_corevelx'] = cat.ix[minihalo_rsid]['corevelx']
                        data_newE.ix[line]['minihalo_corevely'] = cat.ix[minihalo_rsid]['corevely']
                        data_newE.ix[line]['minihalo_corevelz'] = cat.ix[minihalo_rsid]['corevelz']
                         
                        # collect most bound particle
                        iPids = cat.get_all_particles_from_halo(minihalo_rsid)
                        star_pids = iPids[0:1]
                        data_newE.ix[line]['nstars'] = len(star_pids)
                        data_newE.ix[line]['start_pos'] = start_pos
                        allstars=np.r_[allstars,star_pids]
                        start_pos+=len(star_pids)
                 
        fulldataE = pandas.concat((dataE,data_newE),axis=1)
        fulldataE.to_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t')
        f = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'MiniHalos_extantPIDs.dat', 'wb')
        np.array(allstars).tofile(f)
        f.close()
 
    def _read(self,hpath):
        return pandas.read_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t')
 
    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return
 
 
# get parameters and append to otherdata 
def add_data(otherdata,sub_mb, minihalo_loc,subrank,flag):
    global main_branch_time,append_time, mini_halo_time

    mb_root_snap = sub_mb['snap'][0]
    minihalo_vmax = sub_mb[minihalo_loc]['vmax']
    minihalo_snap = sub_mb[minihalo_loc]['snap']
    minihalo_rsid = sub_mb[minihalo_loc]['origid']
    minihalo_mvir = sub_mb[minihalo_loc]['mvir']
    minihalo_posx = sub_mb[minihalo_loc]['posX']
    minihalo_posy = sub_mb[minihalo_loc]['posY']
    minihalo_posz = sub_mb[minihalo_loc]['posZ']
    minihalo_pecvx = sub_mb[minihalo_loc]['pecVX']
    minihalo_pecvy = sub_mb[minihalo_loc]['pecVY']
    minihalo_pecvz = sub_mb[minihalo_loc]['pecVZ']
    minihalo_virialratio = sub_mb[minihalo_loc]['T/|U|']
    minihalo_hostid_MT = sub_mb[minihalo_loc]['pid'] # merger tree ID of host, one level up
    minihalo_rvir = sub_mb[minihalo_loc]['rvir']
    minihalo_spinbullock = sub_mb[minihalo_loc]['spin_bullock']
    minihalo_rs = sub_mb[minihalo_loc]['rs']
    minihalo_scale_of_last_MM = sub_mb[minihalo_loc]['scale_of_last_MM']
    minihalo_Jx = sub_mb[minihalo_loc]['Jx']
    minihalo_Jy = sub_mb[minihalo_loc]['Jy']
    minihalo_Jz = sub_mb[minihalo_loc]['Jz']
    minihalo_xoff = sub_mb[minihalo_loc]['xoff']    
 
    # j is the position of the subhalo i the whole list of subs. As a sub of j, we will report it as -j.
    tt=time.time()
    otherdata=np.r_[otherdata,subrank,sub_mb['origid'][0],minihalo_rsid, minihalo_snap, minihalo_vmax,minihalo_mvir,minihalo_posx,minihalo_posy,minihalo_posz,minihalo_pecvx,minihalo_pecvy,minihalo_pecvz,minihalo_virialratio,minihalo_hostid_MT,minihalo_rvir,minihalo_spinbullock,minihalo_rs,minihalo_scale_of_last_MM,minihalo_Jx,minihalo_Jy,minihalo_Jz,minihalo_xoff,mb_root_snap,flag]
    append_time+=time.time()-tt
    return otherdata
 
  
# this tags sub-subs that have merged with their host subhalo.
# Scenario 1: sub-sub enters main host on its own, then falls into sub, then merges with sub
# Secenario 2: sub-sub falls into sub, then sub and sub-sub system fall into main host together,
# then sub-sub merges with sub
# In both cases, I am currently tagging sub-sub when it first enters the main host.
# j is the position of the subhalo i the whole list of subs. As a sub of j, we will report it as -j.
def auxiliary_add(cat, host_mb, otherdata, host, subline, ii, j, end, min_mass,start_big,flag,mmp_map,non_mmp_map):
    #print end-ii, 'iterations in auxiliary_add'
    global main_branch_time,append_time, mini_halo_time

    while ii!=end:  # ii is backsnap
        tt=time.time()
        merged_subs = host.getNonMMPprogenitors(subline,non_mmp_map) # merged_subs are one step up
        main_branch_time+=time.time()-tt

        host_mb = host_mb[1:] # main branch of our main host
        for subsubline in merged_subs:
            flag_add=0
            # now get main branch of this halo
            tt=time.time()
            sub_mb = host.getMainBranch(subsubline,mmp_map)
            main_branch_time+=time.time()-tt

            # Get maximal mass. Use this to ignore the small halos.
            max_mass = np.max(sub_mb['mvir'])
            if max_mass/cat.h0 < min_mass: 
                sys.stdout.flush()
                continue

            # get infall time, if possible
            ##### ADD CODE HERE
            #  ADD the threshold minihalo identifying code again

            tt=time.time()
            sub_redshifts = (1./sub_mb['scale']) - 1
            mcut = mcrit(Tvir,sub_redshifts)
            # check for phantoms
            index = np.where(np.logical_and(sub_mb['mvir']/cat.h0 > mcut, sub_mb['phantom']==0))[0]
            if len(index)==0:
                 minihalo_loc = None
            else:
                 minihalo_loc = index[0]
            mini_halo_time+=time.time()-tt


            # go recursively deep
            otherdata, start_big = auxiliary_add(cat,host_mb, otherdata,host,subsubline,ii, 10000+j, end=ii+len(sub_mb)-1,min_mass=min_mass,start_big=start_big,flag=flag-1,mmp_map=mmp_map, non_mmp_map=non_mmp_map)
            
            #mask = np.where(sub_mb['phantom']==0)[0]
            # special case: if main branch starts above 10**6. what to do?
            tt=time.time()
            mcut_last = mcrit(Tvir,sub_redshifts[-1])
            if sub_mb['mvir'][-1]/cat.h0 > mcut_last:
                start_big += 1
                print 'encountered halo that started too large'
                flag_add=10
             
            if minihalo_loc == None:
                continue
            mini_halo_time+=time.time()-tt
            ####            
            #print 'adding sub-sub with mass', np.log10(max_mass/cat.h0)
            #print ii, 'value of ii', end, 'value of end', j, 'value of subrank'

            otherdata = add_data(otherdata,sub_mb,minihalo_loc,subrank=-j,flag=flag+flag_add)
            #sys.stdout.flush()
 
        ii+=1
        tt=time.time()
        subline = host.getMMP(subline,mmp_map)
        main_branch_time+=time.time()-tt
    return otherdata,start_big
 
 
# Scenario 1: sub-sub enters main host on its own, then falls into sub, then merges with sub
# Secenario 2: sub-sub falls into sub, then sub and sub-sub system fall into main host together,
# then sub-sub merges with sub
# In both cases, I am currently tagging sub-sub when it first enters the main host.
class DestroyedDataFirstPass(PluginBase):
    def __init__(self):
        super(DestroyedDataFirstPass,self).__init__()
        self.filename='MiniHaloDestroyedDataFirstPass.dat'
        self.filestring='MiniHaloDestroyedDataFirstPass'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = True
        self.xlabel='scale factor' ; self.ylabel='Mass Accreted'    # want these to be adjustable
        self.autofigname='MergerHistory'
        self.min_mass = 10**5.5 
        
    def _analyze(self,hpath):
        global main_branch_time,append_time, mini_halo_time 
        
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
        mtc = haloutils.load_zoom_mtc(hpath) #,haloids=[hostID])
        host = mtc.Trees[0]
        mmp_map = host.get_mmp_map()
        non_mmp_map = host.get_non_mmp_map()
        print 'loaded MTC', time.time()-start_time, 'time to load'
        #### determine end ahead of time.
        tt=time.time()
        host_mb = host.getMainBranch(0,mmp_map)
        main_branch_time+=time.time()-tt
        end = len(host_mb)-1
        ####
        cur_host_line = 0
        i=start # skip to ith level in MT to start
        for k in range(i):
            cur_host_line = host.getMMP(cur_host_line,mmp_map)
         
        while i!=end:
            tt=time.time()
            merged_subs = host.getNonMMPprogenitors(cur_host_line,non_mmp_map)
            main_branch_time+=time.time()-tt
            j=-1; start_pos=0; good=0; start_big=0; too_small=0
 
            otherdata=[]
            host_mb = host.getMainBranch(host.getMMP(cur_host_line,mmp_map),mmp_map) 
            for subline in merged_subs:
                j+=1; flag=0
                tt=time.time()
                sub_mb = host.getMainBranch(subline,mmp_map)
                main_branch_time+=time.time()-tt 
                ##### ADD CODE HERE
                #  ADD the threshold minihalo identifying code again
                
                tt=time.time()
                sub_redshifts = (1./sub_mb['scale']) - 1
                mcut = mcrit(Tvir,sub_redshifts)
                # check for phantoms
                # need mcut indexed to get a single value
                index = np.where(np.logical_and(sub_mb['mvir']/cat.h0 > mcut, sub_mb['phantom']==0))[0]
                if len(index)==0:
                     minihalo_loc = None
                     too_small+=1
                else:
                     minihalo_loc = index[0]
                mini_halo_time+=time.time()-tt

                # even if current branch doesn't start as a minihalo, need to 
                # search all of its progenitors tracing all the way back
                otherdata, start_big = auxiliary_add(cat,host_mb, otherdata,host,subline,i, j, end=i+len(sub_mb)-1,min_mass=self.min_mass,start_big=start_big,flag=-1,mmp_map=mmp_map, non_mmp_map=non_mmp_map)
            
                # special case: if main branch starts above 10**6. what to do?
                tt=time.time()
                mcut_last = mcrit(Tvir,sub_redshifts[-1])
                if sub_mb['mvir'][-1]/cat.h0 > mcut_last:
                    #print 'main branch started while already a minihalo'
                    start_big += 1
                    flag = 10
            
                if minihalo_loc == None:
                    continue
                mini_halo_time+=time.time()-tt
 
                otherdata = add_data(otherdata,sub_mb, minihalo_loc,subrank=j,flag=flag)
                #print j, 'halo in host level', i
                good+=1
                #sys.stdout.flush()
                print j, 'subrank'
                  
            print i, 'host level finished. Time = ', (time.time()-start_time)/60., 'minutes'
            print good,'/',j+1,'had mini halos'
            print too_small,'/',j+1,'always too small'
            print start_big, "number of halos above the threshold to begin with"
            print mini_halo_time/60., 'mini halo time (minutes)'
            print main_branch_time/60., 'main branch time (minutes)'
            print append_time/60., 'append time (minutes)'


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
            dtype = [('sub_rank',dt),('rsid',dt), ('minihalo_rsid',dt), ('minihalo_snap',dt), ('minihalo_vmax',dt),('max_mass',dt),('minihalo_posx',dt),('minihalo_posy',dt),('minihalo_posz',dt),('minihalo_pecvx',dt),('minihalo_pecvy',dt),('minihalo_pecvz',dt),('minihalo_virialratio',dt),('minihalo_hostid_MT',dt),('minihalo_rvir',dt),('minihalo_spinbullock',dt),('minihalo_rs',dt),('minihalo_scale_of_last_MM',dt),('minihalo_Jx',dt),('minihalo_Jy',dt),('minihalo_Jz',dt),('minihalo_xoff',dt),('root_snap',dt),('flags',dt)]
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
        pdtype = ['sub_rank','rsid','minihalo_rsid','minihalo_snap','minihalo_vmax','max_mass','minihalo_posx','minihalo_posy','minihalo_posz','minihalo_pecvx','minihalo_pecvy','minihalo_pecvz','minihalo_virialratio','minihalo_hostid_MT','minihalo_rvir','minihalo_spinbullock','minihalo_rs','minihalo_scale_of_last_MM','minihalo_Jx','minihalo_Jy','minihalo_Jz','minihalo_xoff','root_snap','flags']
        n = len(pdtype)
        import pandas
        return pandas.DataFrame(data.reshape(len(data)/n,n), columns=pdtype)
     
    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return
 
 
class AllDestroyedData(PluginBase):
    def __init__(self):
        super(AllDestroyedData,self).__init__()
        self.filename='AllMiniHaloDestroyedData.dat'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = True
        self.xlabel='' ; self.ylabel=''  # want these to be adjustable
        self.autofigname=''
        self.min_mass = 10**5.5
 
    def _analyze(self,hpath):
        DD = DestroyedDataFirstPass()
        dataD = DD.read(hpath)
        dtype = ['minihalo_mgrav','minihalo_hostid_RS','minihalo_rvmax','minihalo_corevelx','minihalo_corevely','minihalo_corevelz','nstars','start_pos']
        data_newD = pandas.DataFrame(np.zeros((len(dataD),len(dtype)))-1,columns=dtype)
         
        # replace infall data with the minihalo data
        minihalo_dataD = {}
        for minihalo_snap,line in zip(dataD['minihalo_snap'],dataD.index):
            minihalo_dataD.setdefault(minihalo_snap, []).append(line)       
 
        # initialize arrays
        allstars=[]; start_pos=0
 
        for snap in range(haloutils.get_numsnaps(hpath)):
            print snap, 'snap in get extra parameters Destroyed'
            sys.stdout.flush()
            if minihalo_dataD.has_key(snap):
                cat=haloutils.load_rscat(hpath,snap,rmaxcut=False)
                         
                if minihalo_dataD.has_key(snap):
                    for line in minihalo_dataD[snap]:
                        minihalo_rsid = int(dataD.ix[line]['minihalo_rsid'])
                        data_newD.ix[line]['minihalo_mgrav'] = cat.ix[minihalo_rsid]['mgrav']
                        data_newD.ix[line]['minihalo_hostid_RS'] = cat.ix[minihalo_rsid]['hostID']
                        data_newD.ix[line]['minihalo_rvmax'] = cat.ix[minihalo_rsid]['rvmax']
                        data_newD.ix[line]['minihalo_corevelx'] = cat.ix[minihalo_rsid]['corevelx']
                        data_newD.ix[line]['minihalo_corevely'] = cat.ix[minihalo_rsid]['corevely']
                        data_newD.ix[line]['minihalo_corevelz'] = cat.ix[minihalo_rsid]['corevelz']
 
                        iPids = cat.get_all_particles_from_halo(minihalo_rsid)
                        star_pids = iPids[0:1]
                        data_newD.ix[line]['nstars'] = len(star_pids)
                        data_newD.ix[line]['start_pos'] = start_pos
                        allstars=np.r_[allstars,star_pids]
                        start_pos+=len(star_pids)
                 
        fulldataD = pandas.concat((dataD,data_newD),axis=1)
        fulldataD.to_csv(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,sep='\t')
        f = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'MiniHalo_destroyedPIDs.dat', 'wb')
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
        fracs = getFraction(dataE['minihalo_mvir']/cat.h0, haloutils.get_scale_snap(hpath, np.array(dataE['minihalo_snap'],dtype=np.int32)) )
        # need something better than getScale
        nstars = np.array(dataE['nstars'],dtype=np.int32)
        mper=(dataE['minihalo_mvir']/cat.h0*fracs)/nstars
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
        fracs = getFraction(dataD['minihalo_mvir']/cat.h0, haloutils.get_scale_snap(hpath, np.array(dataD['minihalo_snap'],dtype=np.int32)) )
        nstars = np.array(dataD['nstars'],dtype=np.int32)
        mper=(dataD['minihalo_mvir']/cat.h0*fracs)/nstars
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
