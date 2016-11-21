import numpy as np
import pylab as plt
from caterpillaranalysis import *
import readsnapshots.readsnapHDF5_greg as rsg
import haloutils
import time
from scipy import interpolate
from scipy.integrate import quad
import os, subprocess

# make x-axis scale factor not snapshot
# based on get_outfname, hpath is the full path to a halo: i.e. 
# '/bigbang/data/AnnaGroup/caterpillar/halos/H1130025/H1130025_EB_Z127_P7_LN7_LX14_O4_NV4/halos'
#then the data is stored in path/analysis/self.filename  where self.filename is specified in the __init__ function of the plugin.

# re-do tagging based on Alex's new rockstar.
# put old functions in bottom of pile.
# currently used functions on top, with big space.

def distance(posA, posB,boxsize=100.):
    dist = abs(posA-posB)
    tmp = dist > boxsize/2.0
    dist[tmp] = boxsize-dist[tmp]
    if dist.shape == (3,):
        return np.sqrt(np.sum(dist**2))
    else:
        return np.sqrt(np.sum(dist**2,axis=1))

def TagParticles(iSub,iSnap,iPids,iCat,iScale,iMass,snap_z0,hpath):
    #print snap_z0, 'snap_0'
    halopos = np.array(iSub[['posX','posY','posZ']])
    halovel = np.array(iSub[['pecVX','pecVY','pecVZ']])
    iPos = haloutils.load_partblock(hpath,iSnap,'POS ',parttype=1,ids=iPids)
    iVel = np.sqrt(iScale)*haloutils.load_partblock(hpath,iSnap,'VEL ',parttype=1,ids=iPids)  
    dr = iScale*distance(iPos,halopos,boxsize=iCat.boxsize)/iCat.h0 #in MPC physical
    peculiarVEL = iVel-halovel
    Hflow = iCat.H()*(iPos-halopos)*iScale/iCat.h0
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
    star_particles = iPids[boundsort[mask]]#sorted with most bound written first
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


# for tagging particles
class TagExtantPlugin(PluginBase):
    def __init__(self):
        super(TagExtantPlugin,self).__init__()
        self.filename='ExtantPIDs.dat'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = True
        self.xlabel='scale factor' ; self.ylabel='Mass Accreted'    # want these to be adjustable
        self.autofigname='MergerHistory'
        self.min_mass = 10**7.776

    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
        # copy tagExtant code here
        start = 0
        start_time = time.time()
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
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
        mtc = haloutils.load_mtc(hpath,haloids=[hostID])
        print 'loaded mtc'
        sys.stdout.flush()
        host = mtc.Trees[0]
        host_mb = host.getMainBranch(0)
        
        good = 0; start_pos=0; toosmall=0; sub_rank=start-1
        for subRSID in np.array(subs['id']):
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
            peakvmax = np.max(sub_mb['vmax'])
            peakvmax_snap = sub_mb[np.argmax(sub_mb['vmax'])]['snap']

            print sub_rank, 'peakmass', peakmass
            if peakmass/cat.h0 < self.min_mass: # /cat.particle_mass < self.min_particles:
                #print sub_rank, 'subhalo too small'
                sys.stdout.flush()
                toosmall+=1
                continue
            iLoc, iSnap = getInfall(sub_mb, host_mb)
            if iLoc==None:
                continue
            iMass = sub_mb['mvir'][iLoc]
            iVmax = sub_mb['vmax'][iLoc]
            iCat = haloutils.load_rscat(hpath,iSnap,unboundfrac=None,rmaxcut=False)
            iScale = sub_mb['scale'][iLoc]
            iRSID = sub_mb['origid'][iLoc]
            iSub = iCat.ix[iRSID] 
            iPids = iCat.get_all_particles_from_halo(iRSID)
            # iPids here are in bounded order
            #iPids = np.sort(iPids)
            #star_pids = TagParticles(iSub,iSnap,iPids,iCat,iScale,iMass,snap_z0,hpath)
            star_pids = iPids[0:int(np.round(len(iPids)*.1))]
            print sub_rank, len(star_pids), 'num stars tagged'
            allstars=np.r_[allstars,star_pids]
            otherdata=np.r_[otherdata,sub_rank,subRSID,iRSID,start_pos,len(star_pids),peakmass,iMass,peakvmax,peakvmax_snap, sub_mb['mvir'][0],sub_mb['vmax'][0],iSnap,peaksnap,iVmax]
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
        dtype = [('sub_rank',"float64"),('rsid',"float64"),('iRsid',"float64"),('start_pos',"float64"),('nstars',"float64"),('peakmass',"float64"),('infall_mass',"float64"),('vpeak',"float64"),('vpeak_snap',"float64"),('mvir',"float64"),('vmax',"float64"), ('isnap',"float64"),('peaksnap',"float64"),('iVmax',"float64")]
        holder = np.ndarray( (len(data)/14,), dtype=dtype )
        data2 = data.reshape(len(data)/14,14)
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
        self.min_mass = 10**7.776
        # corresponds to 10**7.776 Msun

    def _analyze(self,hpath):
        print 'in tag destroyed'
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
        #subs = cat.get_subhalos_within_halo(hostID)
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
            allstars=[]; otherdata=[]
            host_mb = host.getMainBranch(host.getMMP(cur_host_line)) # this is where it fails ultimately
            for subline in merged_subs:
                j+=1
                sub_mb = host.getMainBranch(subline)
                peakmass = np.max(sub_mb['mvir'])
                peaksnap = sub_mb[np.argmax(sub_mb['mvir'])]['snap']
                peakvmax = np.max(sub_mb['vmax'])
                peakvmax_snap = sub_mb[np.argmax(sub_mb['vmax'])]['snap']
                if peakmass/cat.h0 < self.min_mass: #/cat.particle_mass < self.min_particles:
                    sys.stdout.flush()
                    continue
                iLoc, iSnap = getInfall(sub_mb,host_mb)
                if iLoc == None:
                    print 'subhalo', j, 'is bad in MT. Skipping it'
                    sys.stdout.flush()
                    continue
                
                iMass = sub_mb['mvir'][iLoc]
                iVmax = sub_mb['vmax'][iLoc]
                iCat = haloutils.load_rscat(hpath,iSnap,unboundfrac=None,rmaxcut=False)
                iScale = iCat.scale
                iRSID = sub_mb['origid'][iLoc]
                iSub = iCat.ix[iRSID]
                iPids = iCat.get_all_particles_from_halo(iRSID)
                #iPids = np.sort(iPids)
                #star_pids = TagParticles(iSub,iSnap,iPids,iCat,iScale,iMass,snap_z0,hpath)
                star_pids = iPids[0:int(np.round(len(iPids)*.1))]
                allstars=np.r_[allstars,star_pids]
                otherdata=np.r_[otherdata,j,sub_mb['origid'][0],iRSID,start_pos,len(star_pids),peakmass,iMass,peakvmax,peakvmax_snap,iSnap,peaksnap,iVmax,i]
                start_pos+=len(star_pids)
                print j, 'halo in host level', i
                good+=1
                sys.stdout.flush()
            print i, 'host level finished. Time = ', (time.time()-start_time)/60., 'minutes'
            print good,'/',j+1,'were tagged'
            sys.stdout.flush()
            if not os.path.exists(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed'):
                subprocess.call("mkdir -p "+hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed',shell=True)
            f = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/'+'DestroyedPIDs_'+str(i)+'.dat', 'wb')
            np.array(allstars).tofile(f)
            f.close()
            g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/'+'DestroyedData_'+str(i)+'.dat','wb')
            np.array(otherdata).tofile(g)
            g.close()    
            cur_host_line = host.getMMP(cur_host_line)
            i+=1
        print 'wrote final set of data'

## convert all data into one file
    def combinefiles(self,hpath):
            i = 0; data=[]
            while os.path.exists(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/DestroyedData_'+str(i)+'.dat'):
                tmp = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/DestroyedData_'+str(i)+'.dat')              
                data = np.r_[data,tmp]
                i+=1
            dtype = [('rank',"float64"),('rsid',"float64"),('iRsid',"float64"),('start_pos',"float64"),('nstars',"float64"),('peakmass',"float64"),('infall_mass',"float64"),('vpeak',"float64"),('vpeak_snap',"float64"),('isnap',"float64"),('peaksnap',"float64"),('iVmax',"float64"),('backsnap',"float64")]
            holder = np.ndarray( (len(data)/13,), dtype=dtype )
            data2 = data.reshape(len(data)/13,13)
            for j in range(data2.shape[0]):
                holder[j]=data2[j]
            print np.sum(holder['nstars']), 'nstars from data'
            g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/DestroyedData.dat','wb')
            np.array(holder).tofile(g)
            g.close()    
           
            # Now write ids to file
            i = 0; ids=[]
            while os.path.exists(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/DestroyedPIDs_'+str(i)+'.dat'):
                tmp = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/DestroyedPIDs_'+str(i)+'.dat')
                ids = np.r_[ids,tmp]
                i+=1
            print ids[-1], 'shold not be 0'
            print len(ids), 'length of ids'
            f = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'DestroyedPIDs.dat', 'wb')
            np.array(ids, dtype=np.int64).tofile(f)
            f.close()

            #np.array(bound,dtype=np.float32)
    def _read(self,hpath):
        ids = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'DestroyedPIDs.dat',dtype=np.int64)
        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'DestroyedData.dat')
        dtype = [('rank',"float64"),('rsid',"float64"),('iRsid',"float64"),('start_pos',"float64"),('nstars',"float64"),('peakmass',"float64"),('infall_mass',"float64"),('vpeak',"float64"),('vpeak_snap',"float64"),('isnap',"float64"),('peaksnap',"float64"),('iVmax',"float64"),('backsnap',"float64")]
        holder = np.ndarray( (len(data)/13,), dtype=dtype )
        data2 = data.reshape(len(data)/13,13)
        for i in range(data2.shape[0]):
            holder[i]=data2[i]
        return np.array(ids,dtype=np.int64), holder

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return


def getScale(isnap):
    scales = np.array([0.021276596,0.025114727,0.028952858,0.032790989,0.036629120,0.040467251,0.044305382,0.048143513,0.051981644,0.055819775,0.059657906,0.063496037,0.067334168,0.071172299,0.075010430,0.078848561,0.082686692,0.086524823,0.090362954,0.094201085,0.098039216,0.101877347,0.105715478,0.109553609,0.113391740,0.117229871,0.121068002,0.124906133,0.128744264,0.132582395,0.136420526,0.140258657,0.144096788,0.147934919,0.151773050,0.155611181,0.159449312,0.163287443,0.167125574,0.170963705,0.174801836,0.178639967,0.182478098,0.186316229,0.190154360,0.193992491,0.197830622,0.201668753,0.205506884,0.209345015,0.213183146,0.217021277,0.220859408,0.224697539,0.228535670,0.232373801,0.236211932,0.240050063,0.243888194,0.247726325,0.251564456,0.255402587,0.259240718,0.263078849,0.266916980,0.270755111,0.274593242,0.278431373,0.282269504,0.286107635,0.289945766,0.293783897,0.297622028,0.301460159,0.305298290,0.309136421,0.312974552,0.316812683,0.320650814,0.324488945,0.328327076,0.332165207,0.336003338,0.339841469,0.343679599,0.347517730,0.351355861,0.355193992,0.359032123,0.362870254,0.366708385,0.370546516,0.374384647,0.378222778,0.382060909,0.385899040,0.389737171,0.393575302,0.397413433,0.401251564,0.405089695,0.408927826,0.412765957,0.416604088,0.420442219,0.424280350,0.428118481,0.431956612,0.435794743,0.439632874,0.443471005,0.447309136,0.451147267,0.454985398,0.458823529,0.462661660,0.466499791,0.470337922,0.474176053,0.478014184,0.481852315,0.485690446,0.489528577,0.493366708,0.497204839,0.501042970,0.504881101,0.508719232,0.512557363,0.516395494,0.520233625,0.524071756,0.527909887,0.531748018,0.535586149,0.539424280,0.543262411,0.547100542,0.550938673,0.554776804,0.558614935,0.562453066,0.566291197,0.570129328,0.573967459,0.577805590,0.581643721,0.585481852,0.589319983,0.593158114,0.596996245,0.600834376,0.604672507,0.608510638,0.612348769,0.616186900,0.620025031,0.623863162,0.627701293,0.631539424,0.635377555,0.639215686,0.643053817,0.646891948,0.650730079,0.654568210,0.658406341,0.662244472,0.666082603,0.669920734,0.673758865,0.677596996,0.681435127,0.685273258,0.689111389,0.692949520,0.696787651,0.700625782,0.704463913,0.708302044,0.712140175,0.715978306,0.719816437,0.723654568,0.727492699,0.731330830,0.735168961,0.739007092,0.742845223,0.746683354,0.750521485,0.754359616,0.758197747,0.762035878,0.765874009,0.769712140,0.773550271,0.777388402,0.781226533,0.785064664,0.788902795,0.792740926,0.796579057,0.800417188,0.804255319,0.808093450,0.811931581,0.815769712,0.819607843,0.823445974,0.827284105,0.831122236,0.834960367,0.838798498,0.842636629,0.846474760,0.850312891,0.854151022,0.857989153,0.861827284,0.865665415,0.869503546,0.873341677,0.877179808,0.881017939,0.884856070,0.888694201,0.892532332,0.896370463,0.900208594,0.904046725,0.907884856,0.911722987,0.915561118,0.919399249,0.923237380,0.927075511,0.930913642,0.934751773,0.938589904,0.942428035,0.946266166,0.950104297,0.953942428,0.957780559,0.961618690,0.965456821,0.969294952,0.973133083,0.976971214,0.980809345,0.984647476,0.988485607,0.992323738,0.996161869,1.000000000])
    return scales[isnap]


M10 = 11.590
M11 = 1.195
N10 = .0351
N11 = -0.0247
B10 = 1.376
B11 = -0.826
G10 = 0.608
G11 = 0.329
def getFraction(M, a):
    return 2*N(a)*( (M/M1(a) )**-beta(a) + ( M/M1(a) )**gamma(a) )**-1
    
def M1(a):
    return 10**(M10+M11*(1-a))

def N(a):
    return N10 + N11*(1-a)

def beta(a):
    return B10 + B11*(1-a)

def gamma(a):
    return G10 + G11*(1-a)

class SMFPlugin(PluginBase):
    def __init__(self,Mmin=10**1.0,Mmax=10**8.5,ymin=10**4,ymax=10**13.5):
        super(SMFPlugin,self).__init__()
        self.filename='SMF.dat'
        self.histrange = np.arange(4.0,10.5,0.2)

        self.xmin = Mmin; self.xmax = Mmax
        self.ymin = ymin;  self.ymax = ymax
        self.xlabel = r'$M_{\rm sub} (M_\odot)$'
        self.ylabel = r'$M_{\rm vir} dN/dM_{\rm sub}$'
        self.n_xmin = Mmin/10**12; self.n_xmax = Mmax/10**12
        self.n_ymin = ymin;  self.n_ymax = ymax
        self.n_xlabel = r'$M_{\rm sub}/M_{\rm vir}$'
        self.n_ylabel = self.ylabel
        self.xlog = True; self.ylog = True
        self.autofigname = 'SMF'

    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
        numsnaps = haloutils.get_numsnaps(hpath)
        rscat = haloutils.load_rscat(hpath,numsnaps-1,rmaxcut=False)
        zoomid = haloutils.load_zoomid(hpath)
        
        #subs = rscat.get_all_subhalos_within_halo(zoomid)
    
        # get stellar mass function
        TagExtant = TagExtantPlugin()
        stars, data = TagExtant.read(hpath)
        fracs = getFraction(data['infall_mass']/rscat.h0, getScale( np.array(data['isnap'],dtype=np.int32)) )
        starM = data['infall_mass']/rscat.h0*fracs
        subM = data['mvir']/rscat.h0
        
        x,y = self.MassFunc_dNdM(subM,self.histrange)
        histrange = np.arange(1.0,8.0,0.2)
        xs,ys = self.MassFunc_dNdM(starM,histrange)
        #boundM = np.array(subs['mgrav'])/rscat.h0
        #bx,by = self.MassFunc_dNdM(boundM,self.histrange)
        print np.min(starM), np.max(starM), 'stars'
        print np.min(subM), np.max(subM), 'subs'

        with open(self.get_outfname(hpath),'w') as f:
            for a,b,sa,sb in zip(x,y,xs,ys):
                f.write(str(a)+' '+str(b)+' '+str(sa)+' '+str(sb)+'\n')
    def MassFunc_dNdM(self,masses,histrange):
        """
        Adapted from Greg's MassFunctions code
        """
        numbins = len(histrange) - 1
        hist, r_array = np.histogram(np.log10(masses), bins=histrange)
        x_array = self._getMidpoints(r_array)
        dM = 10.**r_array[1:]-10.**r_array[0:numbins] #Mass size of bins in non-log space
        dNdM = hist/dM
        return 10**x_array, dNdM
    def _getMidpoints(self,bins):
        spacing = bins[1:]-bins[:-1]
        return bins[:-1]+spacing/2.0

    def _read(self,hpath):
        thisfilename = self.get_filename(hpath)
        data = asciitable.read(thisfilename,delimiter=' ')
        #don't return mvir, only mgrav
        return data['col1'],data['col2'],data['col3'],data['col4']
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        x,y,sx,sy = data
        mvir,rvir,vvir=haloutils.load_haloprops(hpath)
        y = y*mvir; sy=sy*mvir
        if normtohost:
            x = x/mvir; sx=sx/mvir
        
        if lx != None:
            #ax.plot(x,y,color=self.colordict[lx],**kwargs)
            ax.plot(sx,sy,color=self.colordict[lx],linestyle='-',**kwargs)
        else:
            #ax.plot(x,y,**kwargs)
            ax.plot(sx,sy,linestyle='--',**kwargs)


# Code taken from TagFunction
class TagMass(PluginBase):
    def __init__(self):
        super(TagMass,self).__init__()
        self.filename='ExtantPIDs_moster.dat'
        self.xmin=1;     self.xmax=400
        self.ymin=10**-5; self.ymax=10*6
        self.xlog= True; self.ylog = True
        self.xlabel='' ; self.ylabel=r''
        self.autofigname=''

    def _analyze(self,hpath):
        # RetagExtant
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
        TagExtant = TagExtantPlugin()
        stars, data = TagExtant.read(hpath)     
        fracs = getFraction(data['infall_mass']/cat.h0, getScale( np.array(data['isnap'],dtype=np.int32)) )
        nstars2 = np.round(data['infall_mass']/cat.particle_mass*.05)
        mper=(data['infall_mass']/cat.h0*fracs)/nstars2
        start_pos = np.array(data['start_pos'],dtype=np.int32)

        newstars=[]; mper_arr=[]
        for i in range(len(data)):
            newstars=np.r_[newstars, stars[start_pos[i]:start_pos[i]+nstars2[i]]]
            mper_arr=np.r_[mper_arr, [mper[i]]*nstars2[i] ]
        # rewrite data properly here
        f = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'ExtantPIDs_moster.dat','wb')
        newstars.tofile(f)
        f.close()
        g=open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'ExtantMass_moster.dat','wb')
        mper_arr.tofile(g)
        g.close()

        # Now retag destroyed data
        # THERE IS AN ERROR IN THIS! start_pos not correct after combining files
        TagDestroyed = TagDestroyedPlugin()
        stars, data = TagDestroyed.read(hpath)
        newstars=[];  mper_arr=[]
        fracs = getFraction(data['infall_mass']/cat.h0, getScale( np.array(data['isnap'],dtype=np.int32)) )
        nstars2 = np.round(data['infall_mass']/cat.particle_mass*.05)
        mper=(data['infall_mass']/cat.h0*fracs)/nstars2
        start_pos = np.array(data['start_pos'],dtype=np.int32)
        for i in range(len(data)):
            newstars=np.r_[newstars, stars[start_pos[i]:start_pos[i]+nstars2[i]]]
            mper_arr=np.r_[mper_arr, [mper[i]]*nstars2[i] ]
        f = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'DestroyedPIDs_moster.dat','wb')
        newstars.tofile(f)
        f.close()
        g=open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'DestroyedMass_moster.dat','wb')
        mper_arr.tofile(g)
        g.close()


    def _read(self,hpath):
        # extant data
        idsE = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'ExtantPIDs_moster.dat')
        massE = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'ExtantMass_moster.dat')
        # destroyed data
        idsD = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'DestroyedPIDs_moster.dat')
        massD = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'DestroyedMass_moster.dat')
        return np.array(idsE,dtype=np.int64), massE, np.array(idsD,dtype=np.int64), massD


class StellarDensProfile(PluginBase):
    def __init__(self):
        super(StellarDensProfile,self).__init__()
        self.filename='StellarDensityProfile.dat'
        self.xmin=10;     self.xmax=220
        self.ymin=10**0; self.ymax=10**6
        self.xlog= True; self.ylog = True
        self.xlabel='r [kpc]' ; self.ylabel=r'$\rho \ [M_\odot/ kpc^3]$'
        self.autofigname='Stellar_Dens_Profile'

    def _analyze(self,hpath):
        print 'SDP analyze'
        tm = TagMass()
        idsE, massE, idsD, massD = tm.read(hpath)
        # combine ids and mass
        ids = np.r_[idsE,idsD]
        mass = np.r_[massE, massD]
        argsort = np.argsort(ids)
        ids = ids[argsort]
        mass = mass[argsort]
        
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID = haloutils.load_zoomid(hpath)
        hosthalo = cat.ix[hostID]
        hostpos = np.array(hosthalo[['posX','posY','posZ']])
        part_pos = haloutils.load_partblock(hpath,snap_z0,"POS ",parttype=1,ids=ids)
        dr = distance(part_pos, hostpos,cat.boxsize)*cat.scale/cat.h0*1000. # in kpc physical
        maxr = float(hosthalo['rvir'])
        minr = 10  #maxr/50.
        binwidth = 0.09
        nbins = np.ceil((np.log10(maxr)-np.log10(minr))/binwidth)
        rarr = 10**np.linspace(np.log10(minr), np.log10(minr)+nbins*binwidth,nbins+1)
        # bin and sum the mass
        bins=np.concatenate(([0],rarr))
        masks=[(dr>bins[i])*(dr<bins[i+1]) for i in range(len(bins)-1)]
        mper = [np.sum(mass[masks[i]]) for i in range(len(masks))]
        m_lt_r = np.cumsum(mper)
        tck = interpolate.splrep(rarr,m_lt_r,k=1)
        rhoarr = interpolate.splev(rarr,tck,der=1)/(4*np.pi*rarr**2)
#        tck = interpolate.splrep(np.log10(rarr),np.log10(m_lt_r),k=3)
#        rhoarr = m_lt_r/rarr*interpolate.splev(np.log10(rarr),tck,der=1)/(4*np.pi*rarr**2)
        f = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,'wb')
        rhoarr.tofile(f)
        rarr.tofile(f)
        f.close()
    
    def _read(self, hpath):
        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename)
        rhoarr = data[0:len(data)/2]
        rarr = data[len(data)/2:]
        return [rhoarr,rarr]

    def _plot(self,hpath,datalist,ax,lx=None,labelon=False,**kwargs):
        rhoarr,rarr = datalist
        ax.plot(rarr, rhoarr)
        #ax.xticks([30,50,100,200],[30,50,100,200])


# gets all stars tagged that are not currently bound to an existing subhalo
def get_stars_not_in_subs(hpath):
    tm = TagMass()
    idsE, massE, idsD, massD = tm.read(hpath)
    ids = np.r_[idsE,idsD]
    argsort = np.argsort(ids)
    mass = np.r_[massE, massD]
    ids = ids[argsort]
    mass = mass[argsort]

    snap_z0 = haloutils.get_numsnaps(hpath)-1
    cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
    hostID = haloutils.load_zoomid(hpath)
    subs = cat.get_all_subhalos_from_halo(hostID)
    mask = subs['mvir'] > 1e8
    subids = np.array(subs[mask]['id'])

    ids_bound = []
    for rsid in subids:
        ids_bound = np.r_[ids_bound, cat.get_all_particles_from_halo(rsid)]

    insubs = np.in1d(ids, ids_bound)
    ids = ids[~insubs]
    mass = mass[~insubs]
    return ids, mass


# for just getting extant data
class ExtantDataPlugin(PluginBase):
    def __init__(self):
        super(ExtantDataPlugin,self).__init__()
        self.filename='ExtantDataOnly.dat'
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
            peak_cat = haloutils.load_rscat(hpath,peak_snap,unboundfrac=None,rmaxcut=False)
            peak_sub = peak_cat.ix[peak_rsid]

            peak_mvir = peak_sub['mvir']
            peak_mgrav = peak_sub['mgrav']
            peak_posx = peak_sub['posX']
            peak_posy = peak_sub['posY']
            peak_posz = peak_sub['posZ']
            peak_corevelx = peak_sub['corevelx']
            peak_corevely = peak_sub['corevely']
            peak_corevelz = peak_sub['corevelz']
            peak_virialratio = peak_sub['T/|U|']
            peak_rvmax = peak_sub['rvmax']
            peak_hostid_RS = peak_sub['hostID']  # rockstar ID of host, one level up
            peak_hostid_MT = sub_mb[peak_loc]['pid'] # merger tree ID of host, one level up
            peak_rvir = peak_sub['rvir']
            peak_spinbullock = peak_sub['spin_bullock']
            
            # Get infall parameters
            infall_cat = haloutils.load_rscat(hpath,iSnap,unboundfrac=None,rmaxcut=False)
            infall_snap = sub_mb[iLoc]['snap']
            infall_scale = sub_mb[iLoc]['scale']
            infall_rsid = sub_mb[iLoc]['origid']
            infall_sub = infall_cat.ix[infall_rsid] 

            infall_vmax = infall_sub['vmax']
            infall_mvir = infall_sub['mvir']
            infall_mgrav = infall_sub['mgrav']
            infall_posx = infall_sub['posX']
            infall_posy = infall_sub['posY']
            infall_posz = infall_sub['posZ']
            infall_corevelx = infall_sub['corevelx']
            infall_corevely = infall_sub['corevely']
            infall_corevelz = infall_sub['corevelz']
            infall_virialratio = infall_sub['T/|U|']
            infall_rvmax = infall_sub['rvmax']
            infall_hostid_RS = infall_sub['hostID']
            infall_hostid_MT = sub_mb[iLoc]['pid']
            infall_rvir = infall_sub['rvir']
            infall_spinbullock = infall_sub['spin_bullock']            

            otherdata=np.r_[otherdata,sub_rank,subRSID,max_mass,max_mass_snap, peak_rsid, peak_snap, peak_vmax,peak_mvir,peak_mgrav,peak_posx,peak_posy,peak_posz,peak_corevelx,peak_corevely,peak_corevelz,peak_virialratio,peak_rvmax,peak_hostid_RS,peak_hostid_MT,peak_rvir,peak_spinbullock,infall_rsid,infall_snap,infall_vmax,infall_mvir,infall_mgrav,infall_posx,infall_posy,infall_posz,infall_corevelx,infall_corevely,infall_corevelz,infall_virialratio,infall_rvmax,infall_hostid_RS,infall_hostid_MT,infall_rvir,infall_spinbullock]
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
        dtype = [('sub_rank',dt),('rsid',dt),('max_mass',dt),('max_mass_snap',dt), ('peak_rsid',dt), ('peak_snap',dt), ('peak_vmax',dt),('peak_mvir',dt),('peak_mgrav',dt),('peak_posx',dt),('peak_posy',dt),('peak_posz',dt),('peak_corevelx',dt),('peak_corevely',dt),('peak_corevelz',dt),('peak_virialratio',dt),('peak_rvmax',dt),('peak_hostid_RS',dt),('peak_hostid_MT',dt),('peak_rvir',dt),('peak_spinbullock',dt),('infall_rsid',dt),('infall_snap',dt),('infall_vmax',dt),('infall_mvir',dt),('infall_mgrav',dt),('infall_posx',dt),('infall_posy',dt),('infall_posz',dt),('infall_corevelx',dt),('infall_corevely',dt),('infall_corevelz',dt),('infall_virialratio',dt),('infall_rvmax',dt),('infall_hostid_RS',dt),('infall_hostid_MT',dt),('infall_rvir',dt),('infall_spinbullock',dt)]
        n = len(dtype)
        holder = np.ndarray( (len(data)/n,), dtype=dtype )
        data2 = data.reshape(len(data)/n,n)
        for i in range(data2.shape[0]):
            holder[i]=data2[i]
        return holder

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return




## need to change most everything to match the above line
class DestroyedDataPlugin(PluginBase):
    def __init__(self):
        super(DestroyedDataPlugin,self).__init__()
        self.filename='DestroyedDataOnly.dat'
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
                # peak_rsid could be a phantom

                peak_cat = haloutils.load_rscat(hpath,peak_snap,unboundfrac=None,rmaxcut=False)
                peak_sub = peak_cat.ix[peak_rsid]
                peak_mvir = peak_sub['mvir']
                peak_mgrav = peak_sub['mgrav']
                peak_posx = peak_sub['posX']
                peak_posy = peak_sub['posY']
                peak_posz = peak_sub['posZ']
                peak_corevelx = peak_sub['corevelx']
                peak_corevely = peak_sub['corevely']
                peak_corevelz = peak_sub['corevelz']
                peak_virialratio = peak_sub['T/|U|']
                peak_rvmax = peak_sub['rvmax']
                peak_hostid_RS = peak_sub['hostID']  # rockstar ID of host, one level up
                peak_hostid_MT = sub_mb[peak_loc]['pid'] # merger tree ID of host, one level up
                peak_rvir = peak_sub['rvir']
                peak_spinbullock = peak_sub['spin_bullock']

                # Get infall parameters
                infall_cat = haloutils.load_rscat(hpath,iSnap,unboundfrac=None,rmaxcut=False)
                infall_snap = sub_mb[iLoc]['snap']
                infall_scale = sub_mb[iLoc]['scale']
                infall_rsid = sub_mb[iLoc]['origid']
                infall_sub = infall_cat.ix[infall_rsid] 
    
                infall_vmax = infall_sub['vmax']
                infall_mvir = infall_sub['mvir']
                infall_mgrav = infall_sub['mgrav']
                infall_posx = infall_sub['posX']
                infall_posy = infall_sub['posY']
                infall_posz = infall_sub['posZ']
                infall_corevelx = infall_sub['corevelx']
                infall_corevely = infall_sub['corevely']
                infall_corevelz = infall_sub['corevelz']
                infall_virialratio = infall_sub['T/|U|']
                infall_rvmax = infall_sub['rvmax']
                infall_hostid_RS = infall_sub['hostID']
                infall_hostid_MT = sub_mb[iLoc]['pid']
                infall_rvir = infall_sub['rvir']
                infall_spinbullock = infall_sub['spin_bullock']            

                otherdata=np.r_[otherdata,j,sub_mb['origid'][0],i,max_mass,max_mass_snap, peak_rsid, peak_snap, peak_vmax,peak_mvir,peak_mgrav,peak_posx,peak_posy,peak_posz,peak_corevelx,peak_corevely,peak_corevelz,peak_virialratio,peak_rvmax,peak_hostid_RS,peak_hostid_MT,peak_rvir,peak_spinbullock,infall_rsid,infall_snap,infall_vmax,infall_mvir,infall_mgrav,infall_posx,infall_posy,infall_posz,infall_corevelx,infall_corevely,infall_corevelz,infall_virialratio,infall_rvmax,infall_hostid_RS,infall_hostid_MT,infall_rvir,infall_spinbullock]
                print j, 'halo in host level', i
                good+=1
                sys.stdout.flush()
            print i, 'host level finished. Time = ', (time.time()-start_time)/60., 'minutes'
            print good,'/',j+1,'were tagged'
            sys.stdout.flush()
            if not os.path.exists(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed'):
                subprocess.call("mkdir -p "+hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed',shell=True)
            g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/'+'DestroyedDataOnly_'+str(i)+'.dat','wb')
            np.array(otherdata).tofile(g)
            g.close()    
            cur_host_line = host.getMMP(cur_host_line)
            i+=1
        print 'wrote final set of data'

## convert all data into one file
    def combinefiles(self,hpath):
            i = 0; data=[]
            while os.path.exists(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/DestroyedDataOnly_'+str(i)+'.dat'):
                tmp = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/DestroyedDataOnly_'+str(i)+'.dat')              
                data = np.r_[data,tmp]
                i+=1
            dtype = [('rank',"float64"),('rsid',"float64"),('backsnap',"float64"),  ]
            holder = np.ndarray( (len(data)/13,), dtype=dtype )
            data2 = data.reshape(len(data)/13,13)
            for j in range(data2.shape[0]):
                holder[j]=data2[j]
            print np.sum(holder['nstars']), 'nstars from data'
            g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/DestroyedData.dat','wb')
            np.array(holder).tofile(g)
            g.close()    
           
            # Now write ids to file
            i = 0; ids=[]
            while os.path.exists(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/DestroyedPIDs_'+str(i)+'.dat'):
                tmp = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/Destroyed/DestroyedPIDs_'+str(i)+'.dat')
                ids = np.r_[ids,tmp]
                i+=1
            print ids[-1], 'shold not be 0'
            print len(ids), 'length of ids'
            f = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'DestroyedPIDs.dat', 'wb')
            np.array(ids, dtype=np.int64).tofile(f)
            f.close()

            #np.array(bound,dtype=np.float32)
    def _read(self,hpath):
        ids = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'DestroyedPIDs.dat',dtype=np.int64)
        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+'DestroyedData.dat')
        dtype = [('rank',"float64"),('rsid',"float64"),('iRsid',"float64"),('start_pos',"float64"),('nstars',"float64"),('peakmass',"float64"),('infall_mass',"float64"),('vpeak',"float64"),('vpeak_snap',"float64"),('isnap',"float64"),('peaksnap',"float64"),('iVmax',"float64"),('backsnap',"float64")]
        holder = np.ndarray( (len(data)/13,), dtype=dtype )
        data2 = data.reshape(len(data)/13,13)
        for i in range(data2.shape[0]):
            holder[i]=data2[i]
        return np.array(ids,dtype=np.int64), holder

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        return
