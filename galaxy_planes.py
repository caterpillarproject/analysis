import numpy as np
import pylab as plt
import os,sys,subprocess,time,functools
import pandas as pd
import cPickle as pickle

import haloutils
import abundmatch,stellarmass
#from SAMs_old import SimpleSAMBasePlugin
from fast_SAMs import FastSAMBasePlugin
from caterpillaranalysis import PluginBase,MassAccrPlugin
from scipy import linalg
import angmom

def plane_tabfn(hpath):
    if hpath==None: return None
    if not haloutils.check_last_rockstar_exists(hpath): return None
    names = ['hid','corr_enhance30','corr_enhance45','conc','mass','spin','scale_of_last_MM','c_to_a']
    formats = [int,float,float,float,float,float,float,float]
    samnames = ['Ni11','Np11','L0i1','L1i1','L0p1','L1p1','L0m1','L1m1']
    samvals = ['ba','ca','rperp','rpar','Nsats']
    for samname in samnames:
        names   += [samname+'_'+val for val in samvals]
        formats += [float,float,float,float,int]

    rscat = haloutils.load_rscat(hpath,haloutils.get_lastsnap(hpath))
    zoomid= haloutils.load_zoomid(hpath)
    host = rscat.ix[zoomid]
    hid = haloutils.get_parent_hid(hpath)

    plug = angmom.AngMomCorrelationPlugin()
    try:
        bins,wlist,logMpeakcutarr = plug.read(hpath)
        w = wlist[0]
        corr_enhance30 = plug.correlation_enhancement(w,30.)
        corr_enhance45 = plug.correlation_enhancement(w,45.)
    except IOError,TypeError:
        return None
    
    plug = MassAccrPlugin()
    mb = plug.read(hpath)
    mbhost = mb[-1]

    data = [hid,corr_enhance30,corr_enhance45,host['rvir']/host['rs'],host['mgrav']/rscat.h0,host['spin'],mbhost['scale_of_last_MM'],host['c_to_a2']]

    plug = SatellitePlanesPlugin()
    df,satix = plug.read(hpath)
    for samname in samnames:
        row = df.ix[samname]
        data += [row['ba'],row['ca'],row['rperp'],row['rpar'],row['n_th90']]
    return data,names,formats

def process_mcconnachie_MW():
    import astropy.coordinates as coord
    import astropy.units as u
    from astropy.table import Table
    from astropy.coordinates import SkyCoord    
    plug = SatellitePlanesPlugin()
    tab = Table.read('planes_data/J_AJ_144_4_catalog-150522.vot')    

    def get_coords(tab):
        ra = coord.Angle(tab['RAJ2000'],unit=u.hourangle)
        dec= coord.Angle(tab['DEJ2000'],unit=u.deg)
        dist = tab['D']
        c = SkyCoord(ra=ra,dec=dec,distance=dist)
        return c

    MWtab = tab[tab['SubG']=='MW'][1:]
    keepcolumns = 'Name','RAJ2000','DEJ2000','D','E_D','e_D','D_MW_','VMag','Mass','Mdyn'
    MWtab = MWtab[keepcolumns]
    MWgalnames = ['LMC','SMC','Fornax','Ursa Minor','Carina','Draco','Sagittarius dSph','Sculptor'\
,'Sextans (I)','Leo I','Leo II'] #'Canes Venatici (I)'
    iigood = np.zeros(len(MWtab)).astype(bool)
    for i,row in enumerate(MWtab):
        if row['Name'] in MWgalnames:
            #print row['Name']
            iigood[i] = True
    MWtab = MWtab[iigood]
    MWtab.sort('D_MW_')
    cMW = get_coords(MWtab).galactocentric
    pos = np.zeros((len(cMW),3))
    pos[:,0] = cMW.x
    pos[:,1] = cMW.y
    pos[:,2] = cMW.z
    vel = np.zeros((len(cMW),3))

    ba,ca,rperp,rpar,ntheta,U = plug.calculate_plane_params(pos,vel,return_vecs=True)
    with open('MW_plane.p','w') as f:
        pickle.dump({'ba':ba,'ca':ca,'rperp':rperp,'rpar':rpar,'U':U,'pos':pos,'MWtab':MWtab},f)

class SatellitePlanesPlugin(PluginBase):
    def __init__(self,usesvd=False):
        super(SatellitePlanesPlugin,self).__init__()
        self.h0 = 0.6711
        self.filename='satsamplanes2.p'

        self.mbplug = MassAccrPlugin()
        self.samplug = FastSAMBasePlugin()#SimpleSAMBasePlugin()
        self.thetaarr = [10,20,30,40,50,60,70,80,90]
        self.thetaarrnames = ['n_th'+str(theta) for theta in self.thetaarr]
        self.n1thetaarrnames = ['n1_th'+str(theta) for theta in self.thetaarr]
    
        GK14 = abundmatch.GK14AbundMatch()
        B13 = abundmatch.Behroozi13AbundMatch()
        M13 = abundmatch.Moster13AbundMatch()
        self.AMlist = [GK14,M13,B13]
        self.Mstar_cutlist = [10.**5.5,10**5] #Msun, VMag < -8.8
        self.Mz8_cut = 10.**8

        self.samnames = ['Ni11','Ni15','Ni25',
                         'Np11','Np15','Np25',
                         #'RNi11','RNi15','RNi25',
                         'L0p0','L0i0','L0p1','L0i1','L0p2','L0i2',
                         'L1p0','L1i0','L1p1','L1i1','L1p2','L1i2',
                         'L0m0','L0m1','L0m2','L1m0','L1m1','L1m2']
        self.samixfns = [self.ixfn_Ni(11),self.ixfn_Ni(15),self.ixfn_Ni(25),
                         self.ixfn_Np(11),self.ixfn_Np(15),self.ixfn_Np(25),
                         #self.ixfn_RNi(11),self.ixfn_RNi(15),self.ixfn_RNi(25),
                         self.ixfn_Lp(0,0),self.ixfn_Li(0,0),
                         self.ixfn_Lp(0,1),self.ixfn_Li(0,1),
                         self.ixfn_Lp(0,2),self.ixfn_Li(0,2),
                         self.ixfn_Lp(1,0),self.ixfn_Li(1,0),
                         self.ixfn_Lp(1,1),self.ixfn_Li(1,1),
                         self.ixfn_Lp(1,2),self.ixfn_Li(1,2),
                         self.ixfn_Lm(0,0),self.ixfn_Lm(0,1),self.ixfn_Lm(0,2),
                         self.ixfn_Lm(1,0),self.ixfn_Lm(1,1),self.ixfn_Lm(1,2)]
        self.samlist = zip(self.samnames,self.samixfns)

        if usesvd:
            self.calculate_plane_params = self.calculate_plane_params_svd
        else:
            self.calculate_plane_params = self.calculate_plane_params_eig

        if os.path.exists('MW_plane.p'):
            with open('MW_plane.p','r') as f:
                #ba,ca,rperp,rpar,U,pos,MWtab
                self.mwdat = pickle.load(f)
        else:
            self.mwdat = None

    def _ixfn_Ni(self,subs,N):
        this_subs = subs.copy()
        this_subs.sort('infall_vmax')[::-1]
        if N < len(this_subs): return this_subs.index[0:N]
        return this_subs.index
    def ixfn_Ni(self,N):
        """Returns ixfn for first N ranked by infall_vmax"""
        return functools.partial(self._ixfn_Ni,N=N)
    def _ixfn_Np(self,subs,N):
        this_subs = subs.copy()
        this_subs.sort('peak_vmax')[::-1]
        if N < len(this_subs): return this_subs.index[0:N]
        return this_subs.index
    def ixfn_Np(self,N):
        """Returns ixfn for first N ranked by infall_vmax"""
        return functools.partial(self._ixfn_Np,N=N)
    def _ixfn_RNi(self,subs,N):
        this_subs = subs.copy()
        this_subs.sort('infall_vmax')[::-1]
        this_subs = this_subs[this_subs['z8_mvir']>=self.Mz8_cut]
        if N < len(this_subs): return this_subs.index[0:N]
        return this_subs.index
    def ixfn_RNi(self,N):
        """Returns ixfn for first N ranked by infall_vmax cut by reionization (z=8) mass"""
        return functools.partial(self._ixfn_RNi,N=N)

    def _ixfn_L(self,subs,AM,whichmass,Mstar_cut):
        this_subs = subs.copy()
        if whichmass=='i':
            Mstar = AM.get_Mstar(np.array(this_subs['infall_mvir']))
        elif whichmass=='p':
            Mstar = AM.get_Mstar(np.array(this_subs['peak_mvir']))
        elif whichmass=='m':
            Mstar = AM.get_Mstar(np.array(this_subs['max_mass']))
        return this_subs.index[Mstar>Mstar_cut]
    def ixfn_Lp(self,whichMstar_cut,whichAM):
        Mstar_cut = self.Mstar_cutlist[whichMstar_cut]
        AM = self.AMlist[whichAM]
        return functools.partial(self._ixfn_L,AM=AM,whichmass='p',Mstar_cut=Mstar_cut)
    def ixfn_Li(self,whichMstar_cut,whichAM):
        Mstar_cut = self.Mstar_cutlist[whichMstar_cut]
        AM = self.AMlist[whichAM]
        return functools.partial(self._ixfn_L,AM=AM,whichmass='i',Mstar_cut=Mstar_cut)
    def ixfn_Lm(self,whichMstar_cut,whichAM):
        Mstar_cut = self.Mstar_cutlist[whichMstar_cut]
        AM = self.AMlist[whichAM]
        return functools.partial(self._ixfn_L,AM=AM,whichmass='m',Mstar_cut=Mstar_cut)

    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
        subs = self.samplug.read(hpath)
        try:
            if subs==None:
                raise IOError("No SAMs")
        except TypeError:
            pass
        mb = self.mbplug.read(hpath)
        host = mb[-1]
        hpos = np.array([host['x'],host['y'],host['z']])
        hvel = np.array([host['vx'],host['vy'],host['vz']])
        hLmom= np.array([host['Jx'],host['Jy'],host['Jz']])
        hA   = np.array([host['A[x]'],host['A[y]'],host['A[z]']])
        
        names = ['ba','ca','rperp','rpar']+['n_th'+str(theta) for theta in self.thetaarr]+['n1_th'+str(theta) for theta in self.thetaarr]
        data = dict(zip(names,[[] for name in names]))
        satindexdict = {}
        for samname,ixfn in self.samlist:
            satix = ixfn(subs)
            satindexdict[samname] = satix
            this_subs = subs.ix[satix].copy()
            spos = (np.array(this_subs[['posX','posY','posZ']])-hpos).reshape((-1,3))
            svel = (np.array(this_subs[['pecVX','pecVY','pecVZ']])-hvel).reshape((-1,3))
            ba,ca,rperp,rpar,ntheta,n1theta = self.calculate_plane_params(spos,svel,getn1theta=True)
            rperp = rperp*1000/self.h0; rpar = rpar*1000/self.h0
            data['ba'].append(ba); data['ca'].append(ca)
            data['rperp'].append(rperp); data['rpar'].append(rpar)
            for i,theta in enumerate(self.thetaarr):
                data['n_th'+str(theta)].append(ntheta[i])
            for i,theta in enumerate(self.thetaarr):
                data['n1_th'+str(theta)].append(n1theta[i])
        mydict = dict(zip(names,data))
        df = pd.DataFrame(data=data,columns=names,index=self.samnames)
        #df.to_csv(self.get_outfname(hpath))
        with open(self.get_outfname(hpath),'w') as f:
            pickle.dump([df,satindexdict],f)

    def _read(self,hpath):
        #df = pd.read_csv(self.get_outfname(hpath),index_col=0)
        with open(self.get_outfname(hpath),'r') as f:
            df,satindexdict = pickle.load(f)
        return df,satindexdict

    def calc_eig(self,satpos):
        #eigenvalues \propto 1/a^2, 1/b^2, 1/c^2 in that order
        I = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                I[i,j] = np.sum(satpos[:,i]*satpos[:,j])
        evals, evecs = linalg.eigh(I)
        return evals,evecs

    def calc_svd(self,satpos,center=False):
        N = float(len(pos))
        mean = np.sum(pos,axis=0)/N
        if center:
            pos = pos - mean
        U,s,Vh = linalg.svd(pos.T)
        if len(s) != 3: 
            raise ValueError("Degenerate SVD (probably too few satellites)")
        return s,U
        
    def calculate_plane_params_eig(self,pos,vel,return_vecs=False,getn1theta=False):
        if len(pos)==1: return np.nan,np.nan,np.nan,np.nan,[1 if theta==90 else -1 for theta in self.thetaarr]
        evals,evecs = self.calc_eig(pos)
        c,b,a = np.sqrt(evals) #the "wrong" way
        vc = evecs[:,0]; vb = evecs[:,1]; va = evecs[:,2]
        ba,ca,rperp,rpar = self._calculate_plane_thickness(pos,vel,a,b,c,va,vb,vc)
        if getn1theta:
            ntheta,n1theta = self._calculate_Lmom_deviation(pos,vel,vc,getn1theta=True)
        else:
            ntheta = self._calculate_Lmom_deviation(pos,vel,vc)
        if return_vecs:
            U = np.vstack([va,vb,vc]).T #U[:,0] is va

        if return_vecs and getn1theta:
            return ba,ca,rperp,rpar,ntheta,n1theta,U
        elif getn1theta:
            return ba,ca,rperp,rpar,ntheta,n1theta
        elif return_vecs:
            return ba,ca,rperp,rpar,ntheta,U
        else:
            return ba,ca,rperp,rpar,ntheta
        
    def calculate_plane_params_svd(self,pos,vel,return_vecs=False):
        if len(pos)==1: return np.nan,np.nan,np.nan,np.nan,[-1 for theta in self.thetaarr]
        s,U = self.calc_svd(pos)
        a,b,c = np.sqrt(s)
        va = U[:,0]; vb = U[:,1]; vc = U[:,2]
        ba,ca,rperp,rpar = self._calculate_plane_thickness(pos,vel,a,b,c,va,vb,vc)
        ntheta = self._calculate_Lmom_deviation(pos,vel,vc)
        if return_vecs:
            return ba,ca,rperp,rpar,ntheta,U
        return ba,ca,rperp,rpar,ntheta
        
    def _calculate_plane_thickness(self,pos,vel,a,b,c,va,vb,vc):
        N = float(len(pos))
        ba = b/a; ca = c/a
        U = np.vstack([va,vb,vc]).T #U[:,i] is vi
        rotcoord = (U.T.dot(pos.T)).T #U^T * x changes coordinates
        u1 = rotcoord[:,0]; u2 = rotcoord[:,1]; u3 = rotcoord[:,2]
        rperp_rms =np.sqrt(np.sum(u3**2)/N)
        rpar_rms  =np.sqrt(np.sum(u1**2)/N)
        return ba,ca,rperp_rms,rpar_rms

    def _calculate_Lmom_deviation(self,pos,vel,pole,getn1theta=False):
        pole = pole/np.sqrt(np.sum(pole**2))
        Lmom = np.cross(pos,vel)
        Lmomnorm = np.sqrt(np.sum(Lmom**2,1))
        theta_deg = np.arccos(np.abs(np.dot(Lmom,pole)/Lmomnorm))*180./np.pi
        n_theta = [np.sum(theta_deg < theta0) for theta0 in self.thetaarr]
        if getn1theta:
            theta_deg = np.arccos(np.dot(Lmom,pole)/Lmomnorm)*180./np.pi
            n1_theta = [max(np.sum(theta_deg < theta0),np.sum(theta_deg > 180-theta0)) for theta0 in self.thetaarr]
            return n_theta,n1_theta
        return n_theta

    def _calculate_Lmomtot_deviation(self,pos,vel,pole):
        pole = pole/np.sqrt(np.sum(pole**2))
        Lmom = np.cross(pos,vel)
        Lmomtot = np.sum(Lmom,axis=0)
        Lmomtot = Lmomtot/np.sqrt(np.sum(Lmomtot**2))
        theta_deg = np.arccos(np.dot(pole,Lmomtot))*180./np.pi
        return theta_deg

    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        raise NotImplementedError

class StochasticSatellitePlanesPlugin(SatellitePlanesPlugin):
    def __init__(self,usesvd=False):
        super(StochasticSatellitePlanesPlugin,self).__init__()
        self.filename='stochsatsams.p'
        self.Nmc = 10**3
    def _analyze(self,hpath):
        # TODO
        # for each in Nmc
        #    get random Mstar
        #    for Mstar_cut in Mstar_cutlist
        #        get satellites
        #        compute/save planes
        raise NotImplementedError

class PlaneEvolutionPlugin(SatellitePlanesPlugin):
    def __init__(self,whichsam='L0p1',verbose=False,usesvd=False):
        super(PlaneEvolutionPlugin,self).__init__()
        self.whichsam = whichsam
        self.filename='planeevol_'+whichsam+'.p'
        self.satplug = SatellitePlanesPlugin()
        self.verbose=verbose

        self.xmin = 0; self.xmax = 14
        self.ymin = 0; self.ymax = 1
        self.xlabel = r'$t\ (\rm{Gyr})$'
        self.ylabel = r'$c/a,\ b/a$'
        self.xlog=False; self.ylog=False
        self.autofigname='planeevolution'

    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
        df,satix = self.satplug.read(hpath)
        row = df.ix[self.whichsam]
        satids = satix[self.whichsam]
        if self.verbose: print "loading mtc"; start = time.time()
        mtc = haloutils.load_zoom_mtc(hpath,indexbyrsid=True)
        zoomid = haloutils.load_zoomid(hpath)
        satmbs = []
        max_minsnap = 0
        if self.verbose: print "done! {0:.1f}".format(time.time()-start)
        if self.verbose: print "finding satmb's"; start = time.time()
        for satid in satids:
            mb = mtc[satid].getMainBranch()
            satmbs.append(mb)
            max_minsnap = max(max_minsnap,np.min(mb['snap']))
        for j,(satid,satmb) in enumerate(zip(satids,satmbs)):
            ii = satmb['snap'] >= max_minsnap
            satmb = satmb[ii]
            satmbs[j] = satmb[::-1] #order from small to large snap
        hostmb = mtc[zoomid].getMainBranch()
        ii = hostmb['snap'] >= max_minsnap
        hostmb = hostmb[ii][::-1] #order from small to large snap

        if self.verbose: print "done! {0:.1f}".format(time.time()-start)
        if self.verbose: print "calculating planes"
        numdata = haloutils.get_numsnaps(hpath)-max_minsnap
        all_snaps = np.zeros(numdata)
        all_satpos = np.zeros((numdata,len(satids),3))
        all_ba = np.zeros(numdata)
        all_ca = np.zeros(numdata)
        all_rperp = np.zeros(numdata)
        all_rpar = np.zeros(numdata)
        vel = np.zeros((len(satids),3))
        for i in range(numdata):
            all_snaps[i] = i+max_minsnap
            host = hostmb[i]
            hpos = np.array([host['posX'],host['posY'],host['posZ']])
            for j,satmb in enumerate(satmbs):
                sub = satmb[i]
                spos = np.array([sub['posX'],sub['posY'],sub['posZ']])-hpos
                all_satpos[i,j,:] = spos
            satpos = all_satpos[i,:,:]
            ba,ca,rperp,rpar,ntheta = self.calculate_plane_params(satpos,vel)
            all_ba[i] = ba; all_ca[i] = ca
            all_rperp[i] = rperp; all_rpar[i] = rpar

        all_infall_snaps = np.zeros(len(satids),dtype=int)
        subs = self.samplug.read(hpath)
        for j,satid in enumerate(satids):
            all_infall_snaps[j] = int(subs.ix[satid]['infall_snap'])

        with open(self.get_outfname(hpath),'w') as f:
            pickle.dump([all_snaps,all_ba,all_ca,all_rperp,all_rpar,all_satpos,all_infall_snaps],f)

    def _read(self,hpath):
        with open(self.get_outfname(hpath),'r') as f:
            data = pickle.load(f)
        return data

    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        assert lx==14 or lx==None
        all_snaps,all_ba,all_ca,all_rperp,all_rpar,all_satpos,all_infall_snaps = data
        t = haloutils.get_t_snap(hpath,all_snaps)
        max_infall_time = haloutils.get_t_snap(hpath,np.max(all_infall_snaps))
        ax.plot(t,all_ca,**kwargs)
        ax.plot(t,all_ba,**kwargs)
        ax.plot([self.xmin,self.xmax],[0.2,0.2],'k:')
        ax.plot([max_infall_time,max_infall_time],[self.ymin,self.ymax],'r:')

if __name__=="__main__":
    RECALC=True
    plug = SatellitePlanesPlugin()
    alldf = {}
    allix = {}
    for hid in haloutils.cid2hid.values():
        hpath = haloutils.get_hpath_lx(hid,14)
        if hpath==None: continue
        try:
            df,satindexdict = plug.read(hpath,recalc=RECALC)
        except TypeError:
            continue
        catnum = haloutils.hid2catnum[hid]
        alldf[catnum]=df
        allix[catnum]=satindexdict
