import numpy as np
import pylab as plt
from caterpillaranalysis import *
import readsnapshots.readsnapHDF5_greg as rsg
import haloutils
import time
from scipy import interpolate
import TidalMethods as TM
from scipy.integrate import quad

# make x-axis scale factor not snapshot
# based on get_outfname, hpath is the full path to a halo: i.e. 
# '/bigbang/data/AnnaGroup/caterpillar/halos/H1130025/H1130025_EB_Z127_P7_LN7_LX14_O4_NV4/halos'
#then the data is stored in path/analysis/self.filename  where self.filename is specified in the __init__ function of the plugin.


# problems - automatic analysis failed...
# I should not get this error! 

# its in the read function
# autocalc is set to true.
# fails on self._read(hpath) in function read.
# this is becauses its trying to read lx=11,12,13 also.

def H(a,Om=.25,Ol=.75,h=0.73):
    to_inv_sec = 3.241*10**-20
    Ok = 1-Om-Ol
    return 100*h*(Om*(1/a)**3 + Ok*(1/a)**2 + Ol)**.5

def TagParticles(iSub,iSnap,iPids,iCat,iScale,iMass,snap_z0,hpath):
    #print snap_z0, 'snap_0'
    halopos = np.array(iSub[['posX','posY','posZ']])
    halovel = np.array(iSub[['pecVX','pecVY','pecVZ']])
    iPos = haloutils.load_partblock(hpath,iSnap,'POS ',parttype=1,ids=iPids)
    iVel = np.sqrt(iScale)*haloutils.load_partblock(hpath,iSnap,'VEL ',parttype=1,ids=iPids)  
    dr = iScale*TM.distance(iPos,halopos,boxsize=iCat.boxsize)/iCat.h0 #in MPC physical
    peculiarVEL = iVel-halovel
    Hflow = H(iScale, h=iCat.h0)*(iPos-halopos)*iScale/iCat.h0
    physicalVEL = peculiarVEL+Hflow
    vel = np.sqrt(sum((physicalVEL**2).T))
    U = PotentialE(dr,iCat) # dr should be physical, no little h
    T = .5*vel**2
    Etot = T+U

    boundsort = np.argsort(Etot)
    n_tag = int(np.round(len(dr)*.1))
    
    # tagging code relic of SIDM tagging
    mask = np.array([False]*len(boundsort))
    totag = np.arange(n_tag)
    mask[totag] = True
    star_particles = iPids[boundsort[mask]]
    return star_particles


def PotentialE(dr,cat):
    G = 1.326*10**11 # in km^3/s^2/Msun
    mpc_to_km = 3.086*10**19

    rarr = 10**np.linspace(np.log10(min(dr))-.01, np.log10(max(dr))+.01,70) # in Mpc
    h_r, x_r = np.histogram(dr, bins=np.concatenate(([0],rarr)))
    m_lt_r = np.cumsum(h_r)*cat.particle_mass/cat.h0
    tck = interpolate.splrep(rarr,m_lt_r) # gives mass in Msun
    def Ufunc(x):
        return interpolate.splev(x,tck)/(x**2)
    
    # do it even faster by using an interpolative function
    # for computing potential energy
    # pick 60 - 100 data points
    # compute potential for all, then use an interpolation scheme
    U = np.zeros(len(rarr))
    for i in range(len(rarr)):
        r = rarr[i]
        if r > max(dr)+.05:
            print 'warning - particle outside of halo. likely inaccurate PE'
            U[i] = -G*m_lt_r[-1]/(r*mpc_to_km)
        else:
            tmp = -G*m_lt_r[-1]/(max(dr)*mpc_to_km)
            U[i] = tmp+G*quad(Ufunc,max(dr),r)[0]/mpc_to_km
    tck2 = interpolate.splrep(rarr,U)
    return interpolate.splev(dr,tck2)


def getInfall(sub_mb, host_mb):
    host_ids = host_mb['id']
    sub_upids = sub_mb['upid']
    if len(host_ids) < len(sub_upids):
        still_sub = np.where(host_ids == sub_upids[0:len(host_ids)])[0]
    else:
        still_sub = np.where(host_ids[0:len(sub_upids)] == sub_upids)[0]
    if len(still_sub) ==0:
        print 'ERROR: "subhalo" never actually a subhalo'
        return None, None
    if still_sub[-1] == len(sub_upids)-1:
        print 'subhalo began as a subhalo'
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
            print 'subhalo is phantom too much.'
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

class FullTreePlugin(MultiPlugin):
    def __init__(self):
        DestroyedTree = DestroyedTreePlugin()
        ExtantTree = ExtantTreePlugin()
        super(FullTreePlugin,self).__init__([DestroyedTree, ExtantTree])
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = True
        self.xlabel='scale factor' ; self.ylabel='Mass Accreted' # want these to be adjustable
        self.autofigname='MergerHistory'

    def _plot(self,hpath,datalist,ax,lx=None,labelon=False,**kwargs):
        destr, ext = datalist
        if destr == None:
            return
        print destr, 'destr', hpath
        nsnaps = 256
        nhalos = np.array([0]*nsnaps)
        masstot = np.array([0]*nsnaps)
        for i in range(nsnaps):
            mask = destr['isnap'] == i
            nhalos[i] = np.sum(mask)
            masstot[i] = np.sum(destr[mask]['infall_mass'])

        nhalosE = np.array([0]*nsnaps)
        masstotE = np.array([0]*nsnaps)  
        for i in range(nsnaps):
            mask = ext['isnap'] ==i
            nhalosE[i]=np.sum(mask)
            masstotE[i]=np.sum(ext[mask]['infall_mass'])        
        ax.bar(np.arange(nsnaps),masstot,label='destroyed',color='blue')
        ax.bar(np.arange(nsnaps),masstotE,label='extant',color='red',bottom=masstot)
        #ax.legend()



class DestroyedTreePlugin(PluginBase):
    def __init__(self):
        super(DestroyedTreePlugin,self).__init__()
        self.filename='DestroyedData.dat'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = True
        self.xlabel='scale factor' ; self.ylabel='Mass Accreted'    # want these to be adjustable
        self.autofigname='MergerHistory'

    # hpath has to exist already for _analyze to be not called to run.
    # suggests hpath must be path to the filename of the data written out to file
    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")

    # must change this to read in data properly
    # should combine all individual files into one single file.
    # then must move each file from my spacebase folder 
    # over to the halos/ folder to be read in later

    # i combined data to all be the same name
    # should just have to read 'this file name'
    # and then convert it to the right dtype array.
    def _read(self,hpath):
        thisfilename = self.get_filename(hpath)
        data = np.fromfile(thisfilename)
        dtype = [('rank',"float64"),('rsid',"float64"),('peakmass',"float64"),('infall_mass',"float64"),('isnap',"float64"),('peaksnap',"float64"),('backsnap',"float64")]
        holder = np.ndarray( (len(data)/7,), dtype=dtype )
        data2 = data.reshape(len(data)/7,7)
        for j in range(data2.shape[0]):
            holder[j]=data2[j]
        return holder


    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        nsnaps=256
        nhalos = np.array([0]*nsnaps)
        masstot = np.array([0]*nsnaps)
        for i in range(nsnaps):
            mask = data['isnap'] == i
            nhalos[i] = np.sum(mask)
            masstot[i] = np.sum(data[mask]['infall_mass'])
        
        ax.bar(np.arange(nsnaps),masstot,label='destroyed',**kwargs)
        # convert snapshot to scale factor

class ExtantTreePlugin(PluginBase):
    def __init__(self):
        super(ExtantTreePlugin,self).__init__()
        self.filename='ExtantData.dat'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = True
        self.xlabel='scale factor' ; self.ylabel='Mass Accreted'    # want these to be adjustable
        self.autofigname='MergerHistory'

    # hpath has to exist already for _analyze to be not called to run.
    # suggests hpath must be path to the filename of the data written out to file
    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")

    # must change this to read in data properly
    # should combine all individual files into one single file.
    # then must move each file from my spacebase folder 
    # over to the halos/ folder to be read in later

    # i combined data to all be the same name
    # should just have to read 'this file name'
    # and then convert it to the right dtype array.
    def _read(self,hpath):
        thisfilename = self.get_filename(hpath)
        data = np.fromfile(thisfilename)
        dtype = [('sub_rank',"float64"),('rsid',"float64"),('peakmass',"float64"),('infall_mass',"float64"),('mvir','float64'),('isnap',"float64"),('peaksnap',"float64")]
        holder = np.ndarray( (len(data)/7,), dtype=dtype )
        data2 = data.reshape(len(data)/7,7)
        for j in range(data2.shape[0]):
            holder[j]=data2[j]
        return holder


    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        nsnaps=256
        nhalos = np.array([0]*nsnaps)
        masstot = np.array([0]*nsnaps)
        for i in range(nsnaps):
            mask = data['isnap'] == i
            nhalos[i] = np.sum(mask)
            masstot[i] = np.sum(data[mask]['infall_mass'])
        
        ax.bar(np.arange(nsnaps),masstot,label='Extant',**kwargs)
        #ax.legend()
        #ax.gca().invert_xaxis()
        # convert snapshot to scale factor
        # combine extant data too


class TagExtantPlugin(PluginBase):
    def __init__(self):
        super(TagExtantPlugin,self).__init__()
        self.filename='ExtantPIDs.dat'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = True
        self.xlabel='scale factor' ; self.ylabel='Mass Accreted'    # want these to be adjustable
        self.autofigname='MergerHistory'
        self.min_particles = 2000 #minimum particles per halo at peak for tagging. Correspons to mvir = 7.776 Msun.

    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
        # copy tagExtant code here
        start = 0; end = 10 # for testing purposes, just tag the first 10

        start_time = time.time()
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        cat = haloutils.load_rscat(hpath,snap_z0)
        hostID1 = int(cat['id'][0:1])
        hostID = haloutils.load_zoomid(hpath)
        if hostID != hostID1:
            print 'host IDs do not match!!'
        hosthalo = cat.ix[hostID]
        subs = cat.get_subhalos_within_halo(hostID)
        
        allstars=[]
        otherdata=[]
        print 'loading mtc'
        sys.stdout.flush()
        mtc = haloutils.load_mtc(hpath)
        print 'loaded mtc'
        sys.stdout.flush()
        host = mtc.Trees[0]
        host_mb = host.getMainBranch(0)
        
        good = 0; start_pos=0; toosmall=0; sub_rank=start-1
        for subRSID in np.array(subs['id']): # EDIT should remove 0:100 after test
            sub_rank+=1
            sub = getSubTree(mtc,subRSID)
            if sub==None:
                print sub_rank, 'subhalo not found in MTCatalogue. Mass: %.4e' %cat.ix[subRSID]['mvir'], 'Time = ', (time.time()-start_time)/60., 'minutes'
                sys.stdout.flush()
                continue
            sub_mb = sub.getMainBranch(0)
            if sub_mb == None:
                print 'subhalo', sub_rank, 'is bad in MT. Skipping it'
                sys.stdout.flush()
                continue # skip to next subhalo
            peakmass = np.max(sub_mb['mvir'])
            peaksnap = sub_mb[np.argmax(sub_mb['mvir'])]['snap']
            print sub_rank, 'peakmass', peakmass
            if peakmass/cat.particle_mass < self.min_particles:
                #print sub_rank, 'subhalo too small'
                sys.stdout.flush()
                toosmall+=1
                continue
            iLoc, iSnap = getInfall(sub_mb, host_mb)
            if iLoc==None:
                continue
            iMass = sub_mb['mvir'][iLoc]
            iCat = haloutils.load_rscat(hpath,iSnap)
            iScale = sub_mb['scale'][iLoc]
            iRSID = sub_mb['origid'][iLoc]
            iSub = iCat.ix[iRSID] 
            iPids = iCat.get_all_particles_from_halo(iRSID)
            iPids = np.sort(iPids)
            star_pids = TagParticles(iSub,iSnap,iPids,iCat,iScale,iMass,snap_z0,hpath)
            print sub_rank, len(star_pids), 'num stars tagged'
            allstars=np.r_[allstars,star_pids]
            otherdata=np.r_[otherdata,sub_rank,subRSID,iRSID,start_pos,len(star_pids),peakmass,iMass,sub_mb['mvir'][0],iSnap,peaksnap]
            start_pos+=len(star_pids)
            print sub_rank, '/', len(subs), 'finished. Time = ', (time.time()-start_time)/60., 'minutes'
            sys.stdout.flush()
            good+=1
        f = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'ExtantPIDs.dat', 'wb')
        np.array(allstars).tofile(f)
        f.close()
        g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'ExtantData.dat','wb')
        np.array(otherdata).tofile(g)
        g.close()
# to test, use haloutils.get_hpath_lx(hid,lx)
        print good, 'halos good out of', len(subs)
        print toosmall, 'num halos too small'

    def _read(self,hpath):
        ids = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'ExtantPIDs.dat')
        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'ExtantData.dat')
        dtype = [('sub_rank',"float64"),('rsid',"float64"),('iRsid',"float64"),('start_pos',"float64"),('nstars',"float64"),('peakmass',"float64"),('infall_mass',"float64"),('mvir',"float64"),('isnap',"float64"),('peaksnap',"float64")]
        holder = np.ndarray( (len(data)/10,), dtype=dtype )
        data2 = data.reshape(len(data)/10,10)
        for i in range(data2.shape[0]):
            holder[i]=data2[i]
        return np.array(ids,dtype=np.int64), holder


    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return


class TagDestroyedPlugin(PluginBase):
    def __init__(self):
        super(TagDestroyedPlugin,self).__init__()
        self.filename='DestroyedPIDs.dat'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = True
        self.xlabel='scale factor' ; self.ylabel='Mass Accreted'    # want these to be adjustable
        self.autofigname='MergerHistory'
        self.min_particles = 2000 #minimum particles per halo at peak for tagging.
        # corresponds to 10**7.776 Msun

    def _analyze(self,hpath):
        print 'in tag destroyed'
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
        # copy tagExtant code here
        start = 0; end = 255 # for testing purposes, just tag the first 10

        start_time = time.time()
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        cat = haloutils.load_rscat(hpath,snap_z0)
        hostID1 = int(cat['id'][0:1])
        hostID = haloutils.load_zoomid(hpath)
        if hostID != hostID1:
            print 'host IDs do not match!!'
        hosthalo = cat.ix[hostID]
        subs = cat.get_subhalos_within_halo(hostID)
        mtc = haloutils.load_mtc(hpath)
        print 'loaded MTC'
        host = mtc.Trees[0]
        cur_host_line = 0
        i=start # skip to ith level in MT to start
        for k in range(i):
            cur_host_line = host.getMMP(cur_host_line)
        
        while i!=end:
            merged_subs = host.getNonMMPprogenitors(cur_host_line)
            j=-1; start_pos=0; good=0
            allstars=[]; otherdata=[]
            host_mb = host.getMainBranch(host.getMMP(cur_host_line))
            for subline in merged_subs:
                j+=1
                sub_mb = host.getMainBranch(subline)
                peakmass = np.max(sub_mb['mvir'])
                peaksnap = sub_mb[np.argmax(sub_mb['mvir'])]['snap']
                if peakmass/cat.particle_mass < self.min_particles:
                    sys.stdout.flush()
                    continue
                iLoc, iSnap = getInfall(sub_mb,host_mb)
                if iLoc == None:
                    print 'subhalo', j, 'is bad in MT. Skipping it'
                    sys.stdout.flush()
                    continue
                
                iMass = sub_mb['mvir'][iLoc]
                iCat = haloutils.load_rscat(hpath,iSnap)
                iScale = iCat.scale
                iRSID = sub_mb['origid'][iLoc]
                iSub = iCat.ix[iRSID]
                iPids = iCat.get_all_particles_from_halo(iRSID)
                iPids = np.sort(iPids)
            
                star_pids = TagParticles(iSub,iSnap,iPids,iCat,iScale,iMass,snap_z0,hpath)
                allstars=np.r_[allstars,star_pids]
                otherdata=np.r_[otherdata,j,sub_mb['origid'][0],iRSID,start_pos,len(star_pids),peakmass,iMass,iSnap,peaksnap,i]
                start_pos+=len(star_pids)
                print j, 'halo in host level', i
                good+=1
                sys.stdout.flush()
            print i, 'host level finished. Time = ', (time.time()-start_time)/60., 'minutes'
            print good,'/',j+1,'were tagged'
            sys.stdout.flush()
            f = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'DestroyedPIDs_'+str(i)+'.dat', 'wb')
            np.array(allstars).tofile(f)
            f.close()
            g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'DestroyedData_'+str(i)+'.dat','wb')
            np.array(otherdata).tofile(g)
            g.close()    
            cur_host_line = host.getMMP(cur_host_line)
            i+=1
        print 'wrote final set of data'

    def _read(self,hpath):
        thisfilename = self.get_filename(hpath)
        data = np.fromfile(thisfilename)

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return

