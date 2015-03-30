import numpy as np
import pylab as plt
from caterpillaranalysis import *
import readsnapshots.readsnapHDF5_greg as rsg
import haloutils
import time
from scipy import interpolate
from scipy.integrate import quad
import os, subprocess

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
        subs = cat.get_subhalos_within_halo(hostID)
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
                print sub_rank, 'subhalo not found in MTCatalogue. Z=0 Mass: %.4e' %cat.ix[subRSID]['mvir'], 'Time = ', (time.time()-start_time)/60., 'minutes'
                sys.stdout.flush()
                continue
            sub_mb = sub.getMainBranch(0)
            if sub_mb == None:
                print 'subhalo', sub_rank, 'main branch not found in MT. Skipping it. Z=0 Mass: %.4e' %cat.ix[subRSID]['mvir']
                sys.stdout.flush()
                continue # skip to next subhalo

            # Get maximal mass. Use this to ignore the small halos.
            max_mass = np.max(sub_mb['mvir'])
            if max_mass/cat.h0 < self.min_mass: 
                #print sub_rank, 'subhalo too small'
                sys.stdout.flush()
                toosmall+=1
                continue
            max_mass_snap = sub_mb[np.argmax(sub_mb['mvir'])]['snap']

            # get infall time, if possible
            iLoc, iSnap = getInfall(sub_mb, host_mb, max_mass)
            if iLoc==None:
                # within getInfall, print how it errored
                continue

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
            #peak_mgrav = peak_sub['mgrav']
            peak_posx = sub_mb[peak_loc]['posX']
            peak_posy = sub_mb[peak_loc]['posY']
            peak_posz = sub_mb[peak_loc]['posZ']
            peak_pecvx = sub_mb[peak_loc]['pecVX']
            peak_pecvy = sub_mb[peak_loc]['pecVY']
            peak_pecvz = sub_mb[peak_loc]['pecVZ']
            peak_virialratio = sub_mb[peak_loc]['T/|U|']
            #peak_rvmax = peak_sub['rvmax']
            #peak_hostid_RS = peak_sub['hostID']  # rockstar ID of host, one level up
            peak_hostid_MT = sub_mb[peak_loc]['pid'] # merger tree ID of host, one level up
            peak_rvir = sub_mb[peak_loc]['rvir']
            peak_spinbullock = sub_mb[peak_loc]['spin_bullock']
            
            # Get infall parameters
            infall_snap = sub_mb[iLoc]['snap']
            infall_scale = sub_mb[iLoc]['scale']
            infall_rsid = sub_mb[iLoc]['origid']

            infall_vmax = sub_mb[iLoc]['vmax']
            infall_mvir = sub_mb[iLoc]['mvir']
            #infall_mgrav = infall_sub['mgrav']
            infall_posx = sub_mb[iLoc]['posX']
            infall_posy = sub_mb[iLoc]['posY']
            infall_posz = sub_mb[iLoc]['posZ']
            infall_pecvx = sub_mb[iLoc]['pecVX']
            infall_pecvy = sub_mb[iLoc]['pecVY']
            infall_pecvz = sub_mb[iLoc]['pecVZ']
            infall_virialratio = sub_mb[iLoc]['T/|U|']
            #infall_rvmax = infall_sub['rvmax']
            #infall_hostid_RS = infall_sub['hostID']
            infall_hostid_MT = sub_mb[iLoc]['pid']
            infall_rvir = sub_mb[iLoc]['rvir']
            infall_spinbullock = sub_mb[iLoc]['spin_bullock']            

            otherdata=np.r_[otherdata,sub_rank,subRSID,max_mass,max_mass_snap, peak_rsid, peak_snap, peak_vmax,peak_mvir,peak_posx,peak_posy,peak_posz,peak_pecvx,peak_pecvy,peak_pecvz,peak_virialratio,peak_hostid_MT,peak_rvir,peak_spinbullock,infall_rsid,infall_snap,infall_vmax,infall_mvir,infall_posx,infall_posy,infall_posz,infall_pecvx,infall_pecvy,infall_pecvz,infall_virialratio,infall_hostid_MT,infall_rvir,infall_spinbullock]
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
        dt = "float64"
        dtype = [('sub_rank',dt),('rsid',dt),('max_mass',dt),('max_mass_snap',dt), ('peak_rsid',dt), ('peak_snap',dt), ('peak_vmax',dt),('peak_mvir',dt),('peak_posx',dt),('peak_posy',dt),('peak_posz',dt),('peak_pecvx',dt),('peak_pecvy',dt),('peak_pecvz',dt),('peak_virialratio',dt),('peak_hostid_MT',dt),('peak_rvir',dt),('peak_spinbullock',dt),('infall_rsid',dt),('infall_snap',dt),('infall_vmax',dt),('infall_mvir',dt),('infall_posx',dt),('infall_posy',dt),('infall_posz',dt),('infall_pecvx',dt),('infall_pecvy',dt),('infall_pecvz',dt),('infall_virialratio',dt),('infall_hostid_MT',dt),('infall_rvir',dt),('infall_spinbullock',dt)]
        n = len(dtype)
        holder = np.ndarray( (len(data)/n,), dtype=dtype )
        data2 = data.reshape(len(data)/n,n)
        for i in range(data2.shape[0]):
            holder[i]=data2[i]
        return holder

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return


## need to change most everything to match the above line
class DestroyedDataFirstPass(PluginBase):
    def __init__(self):
        super(DestroyedDataFirstPass,self).__init__()
        self.filename='DestroyedDataFirstPass'
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
                max_mass_snap = sub_mb[np.argmax(sub_mb['mvir'])]['snap']

                # get infall time, if possible
                iLoc, iSnap = getInfall(sub_mb,host_mb, max_mass)
                if iLoc == None:
                    print 'subhalo', j, 'is bad in MT. Reason to follow.'
                    sys.stdout.flush()
                    continue

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
                #peak_mgrav = peak_sub['mgrav']
                peak_posx = sub_mb[peak_loc]['posX']
                peak_posy = sub_mb[peak_loc]['posY']
                peak_posz = sub_mb[peak_loc]['posZ']
                peak_pecvx = sub_mb[peak_loc]['pecVX']
                peak_pecvy = sub_mb[peak_loc]['pecVY']
                peak_pecvz = sub_mb[peak_loc]['pecVZ']
                peak_virialratio = sub_mb[peak_loc]['T/|U|']
                #peak_rvmax = sub_mb[peak_loc]['rvmax']
                #peak_hostid_RS = sub_mb[peak_loc]['hostID']  # rockstar ID of host, one level up
                peak_hostid_MT = sub_mb[peak_loc]['pid'] # merger tree ID of host, one level up
                peak_rvir = sub_mb[peak_loc]['rvir']
                peak_spinbullock = sub_mb[peak_loc]['spin_bullock']

                # Get infall parameters
                infall_snap = sub_mb[iLoc]['snap']
                infall_scale = sub_mb[iLoc]['scale']
                infall_rsid = sub_mb[iLoc]['origid']
    
                infall_vmax = sub_mb[iLoc]['vmax']
                infall_mvir = sub_mb[iLoc]['mvir']
                #infall_mgrav = infall_sub['mgrav']
                infall_posx = sub_mb[iLoc]['posX']
                infall_posy = sub_mb[iLoc]['posY']
                infall_posz = sub_mb[iLoc]['posZ']
                infall_pecvx = sub_mb[iLoc]['pecVX']
                infall_pecvy = sub_mb[iLoc]['pecVY']
                infall_pecvz = sub_mb[iLoc]['pecVZ']
                infall_virialratio = sub_mb[iLoc]['T/|U|']
                #infall_rvmax = infall_sub['rvmax']
                #infall_hostid_RS = infall_sub['hostID']
                infall_hostid_MT = sub_mb[iLoc]['pid']
                infall_rvir = sub_mb[iLoc]['rvir']
                infall_spinbullock = sub_mb[iLoc]['spin_bullock']            

                otherdata=np.r_[otherdata,j,sub_mb['origid'][0],max_mass,max_mass_snap, peak_rsid, peak_snap, peak_vmax,peak_mvir,peak_posx,peak_posy,peak_posz,peak_pecvx,peak_pecvy,peak_pecvz,peak_virialratio,peak_hostid_MT,peak_rvir,peak_spinbullock,infall_rsid,infall_snap,infall_vmax,infall_mvir,infall_posx,infall_posy,infall_posz,infall_pecvx,infall_pecvy,infall_pecvz,infall_virialratio,infall_hostid_MT,infall_rvir,infall_spinbullock,i]
                #print j, 'halo in host level', i
                good+=1
                sys.stdout.flush()
            print i, 'host level finished. Time = ', (time.time()-start_time)/60., 'minutes'
            print good,'/',j+1,'were tagged'
            sys.stdout.flush()
            if not os.path.exists(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed'):
                subprocess.call("mkdir -p "+hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed',shell=True)
            g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/'+self.filename+'_'+str(i)+'.dat','wb')
            np.array(otherdata).tofile(g)
            g.close()    
            cur_host_line = host.getMMP(cur_host_line)
            i+=1
        print 'wrote final set of data'

## convert all data into one file
    def combinefiles(self,hpath):
            i = 0; data=[]
            while os.path.exists(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/'+self.filename+'_'+str(i)+'.dat'):
                tmp = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/'+self.filename+'_'+str(i)+'.dat')              
                data = np.r_[data,tmp]
                i+=1
            dt = "float64"
            dtype = [('sub_rank',dt),('rsid',dt),('max_mass',dt),('max_mass_snap',dt), ('peak_rsid',dt), ('peak_snap',dt), ('peak_vmax',dt),('peak_mvir',dt),('peak_posx',dt),('peak_posy',dt),('peak_posz',dt),('peak_pecvx',dt),('peak_pecvy',dt),('peak_pecvz',dt),('peak_virialratio',dt),('peak_hostid_MT',dt),('peak_rvir',dt),('peak_spinbullock',dt),('infall_rsid',dt),('infall_snap',dt),('infall_vmax',dt),('infall_mvir',dt),('infall_posx',dt),('infall_posy',dt),('infall_posz',dt),('infall_pecvx',dt),('infall_pecvy',dt),('infall_pecvz',dt),('infall_virialratio',dt),('infall_hostid_MT',dt),('infall_rvir',dt),('infall_spinbullock',dt),('backsnap',dt)]
            n = len(dtype)
            holder = np.ndarray( (len(data)/n,), dtype=dtype )
            data2 = data.reshape(len(data)/n,n)
            for j in range(data2.shape[0]):
                holder[j]=data2[j]
            g = open(hpath+'/'+self.OUTPUTFOLDERNAME+self.filename+'.dat','wb')
            np.array(holder).tofile(g)
            g.close()    

            #np.array(bound,dtype=np.float32)
    def _read(self,hpath):
        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename+'.dat')
        dt = "float64"
        dtype = [('sub_rank',dt),('rsid',dt),('max_mass',dt),('max_mass_snap',dt), ('peak_rsid',dt), ('peak_snap',dt), ('peak_vmax',dt),('peak_mvir',dt),('peak_posx',dt),('peak_posy',dt),('peak_posz',dt),('peak_pecvx',dt),('peak_pecvy',dt),('peak_pecvz',dt),('peak_virialratio',dt),('peak_hostid_MT',dt),('peak_rvir',dt),('peak_spinbullock',dt),('infall_rsid',dt),('infall_snap',dt),('infall_vmax',dt),('infall_mvir',dt),('infall_posx',dt),('infall_posy',dt),('infall_posz',dt),('infall_pecvx',dt),('infall_pecvy',dt),('infall_pecvz',dt),('infall_virialratio',dt),('infall_hostid_MT',dt),('infall_rvir',dt),('infall_spinbullock',dt),('backsnap',dt)]
        n = len(dtype)
        holder = np.ndarray( (len(data)/n,), dtype=dtype )
        data2 = data.reshape(len(data)/n,n)
        for i in range(data2.shape[0]):
            holder[i]=data2[i]
        return np.array(ids,dtype=np.int64), holder

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return
