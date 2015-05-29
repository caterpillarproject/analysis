import numpy as np
import pylab as plt
import os,sys,subprocess,time,functools
import pandas as pd
import cPickle as pickle

import haloutils
import abundmatch,stellarmass
from SAMs_old import SimpleSAMBasePlugin
from caterpillaranalysis import PluginBase,MassAccrPlugin
from scipy import linalg

def plane_tabfn(hpath):
    if hpath==None: return None
    if not haloutils.check_last_rockstar_exists(hpath): return None
    names = ['hid',]

    plug = SatellitePlanesPlugin()
    df = plug.read(hpath)

    for samname in ['Ni11','L0i1','L1i1']:
        row = df.ix[samname]
        Nsats = row['n_th90']
        ba = row['ba']; ca = row['ca']
        rperp = row['rperp']; rpar = row['rpar']

def process_mcconnachie_MW():
    import astropy.coordinates as coord
    import astropy.units as u
    from astropy.table import Table
    from astropy.coordinates import SkyCoord    
    plug = SatellitePlanesPlugin()
    tab = Table.read('J_AJ_144_4_catalog-150522.vot')    

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
        pickle.dump([ba,ca,rperp,rpar,U,pos,MWtab],f)

class SatellitePlanesPlugin(PluginBase):
    def __init__(self,usesvd=False):
        super(SatellitePlanesPlugin,self).__init__()
        self.h0 = 0.6711
        self.filename='satsamplanes.p'

        self.mbplug = MassAccrPlugin()
        self.samplug = SimpleSAMBasePlugin()
        self.thetaarr = [10,20,30,40,50,60,70,80,90]
        self.thetaarrnames = ['n_th'+str(theta) for theta in self.thetaarr]
    
        GK14 = abundmatch.GK14AbundMatch()
        B13 = abundmatch.Behroozi13AbundMatch()
        M13 = abundmatch.Moster13AbundMatch()
        self.AMlist = [GK14,M13,B13]
        self.Mstar_cutlist = [10.**5.5,10**5] #Msun, VMag < -8.8
        self.Mz8_cut = 10.**8

        self.samnames = ['Ni11','Ni15','Ni25',
                         'RNi11','RNi15','RNi25',
                         'L0p0','L0i0','L0p1','L0i1','L0p2','L0i2',
                         'L1p0','L1i0','L1p1','L1i1','L1p2','L1i2']
        self.samixfns = [self.ixfn_Ni(11),self.ixfn_Ni(15),self.ixfn_Ni(25),
                         self.ixfn_RNi(11),self.ixfn_RNi(15),self.ixfn_RNi(25),
                         self.ixfn_Lp(0,0),self.ixfn_Li(0,0),
                         self.ixfn_Lp(0,1),self.ixfn_Li(0,1),
                         self.ixfn_Lp(0,2),self.ixfn_Li(0,2),
                         self.ixfn_Lp(1,0),self.ixfn_Li(1,0),
                         self.ixfn_Lp(1,1),self.ixfn_Li(1,1),
                         self.ixfn_Lp(1,2),self.ixfn_Li(1,2)]
        self.samlist = zip(self.samnames,self.samixfns)

        if usesvd:
            self.calculate_plane_params = self.calculate_plane_params_svd
        else:
            self.calculate_plane_params = self.calculate_plane_params_eig

    def _ixfn_Ni(self,subs,N):
        this_subs = subs.copy()
        this_subs.sort('infall_vmax')[::-1]
        if N < len(this_subs): return this_subs.index[0:N]
        return this_subs.index
    def ixfn_Ni(self,N):
        """Returns ixfn for first N ranked by infall_vmax"""
        return functools.partial(self._ixfn_Ni,N=N)
    def _ixfn_RNi(self,subs,N):
        this_subs = subs.copy()
        this_subs.sort('infall_vmax')[::-1]
        this_subs = this_subs[this_subs['z8_mvir']>=self.Mz8_cut]
        if N < len(this_subs): return this_subs.index[0:N]
        return this_subs.index
    def ixfn_RNi(self,N):
        """Returns ixfn for first N ranked by infall_vmax cut by reionization (z=8) mass"""
        return functools.partial(self._ixfn_RNi,N=N)

    def _ixfn_L(self,subs,AM,useinfall,Mstar_cut):
        this_subs = subs.copy()
        if useinfall:
            Mstar = AM.get_Mstar(np.array(this_subs['infall_mvir']))
        else:
            Mstar = AM.get_Mstar(np.array(this_subs['peak_mvir']))
        return this_subs.index[Mstar>Mstar_cut]
    def ixfn_Lp(self,whichMstar_cut,whichAM):
        Mstar_cut = self.Mstar_cutlist[whichMstar_cut]
        AM = self.AMlist[whichAM]
        return functools.partial(self._ixfn_L,AM=AM,useinfall=False,Mstar_cut=Mstar_cut)
    def ixfn_Li(self,whichMstar_cut,whichAM):
        Mstar_cut = self.Mstar_cutlist[whichMstar_cut]
        AM = self.AMlist[whichAM]
        return functools.partial(self._ixfn_L,AM=AM,useinfall=True,Mstar_cut=Mstar_cut)

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
        
        names = ['ba','ca','rperp','rpar']+['n_th'+str(theta) for theta in self.thetaarr]
        data = dict(zip(names,[[] for name in names]))
        satindexdict = {}
        for samname,ixfn in self.samlist:
            satix = ixfn(subs)
            satindexdict[samname] = satix
            this_subs = subs.ix[satix].copy()
            spos = (np.array(this_subs[['posX','posY','posZ']])-hpos).reshape((-1,3))
            svel = (np.array(this_subs[['pecVX','pecVY','pecVZ']])-hvel).reshape((-1,3))
            ba,ca,rperp,rpar,ntheta = self.calculate_plane_params(spos,svel)
            rperp = rperp*1000/self.h0; rpar = rpar*1000/self.h0
            data['ba'].append(ba); data['ca'].append(ca)
            data['rperp'].append(rperp); data['rpar'].append(rpar)
            for i,theta in enumerate(self.thetaarr):
                data['n_th'+str(theta)].append(ntheta[i])
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
        
    def calculate_plane_params_eig(self,pos,vel,return_vecs=False):
        if len(pos)==1: return np.nan,np.nan,np.nan,np.nan,[1 if theta==90 else -1 for theta in self.thetaarr]
        evals,evecs = self.calc_eig(pos)
        c,b,a = np.sqrt(evals) #the "wrong" way
        vc = evecs[:,0]; vb = evecs[:,1]; va = evecs[:,2]
        ba,ca,rperp,rpar = self._calculate_plane_thickness(pos,vel,a,b,c,va,vb,vc)
        ntheta = self._calculate_Lmom_deviation(pos,vel,vc)
        if return_vecs:
            U = np.vstack([va,vb,vc]).T #U[:,0] is va
            return ba,ca,rperp,rpar,ntheta,U
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

    def _calculate_Lmom_deviation(self,pos,vel,pole):
        pole = pole/np.sqrt(np.sum(pole**2))
        Lmom = np.cross(pos,vel)
        Lmomnorm = np.sqrt(np.sum(Lmom**2,1))
        theta_deg = np.arccos(np.abs(np.dot(Lmom,pole)/Lmomnorm))*180./np.pi
        n_theta = [np.sum(theta_deg < theta0) for theta0 in self.thetaarr]
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

if __name__=="__main__":
    RECALC=False
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
        alldf[hid]=df
        allix[hid]=satindexdict
