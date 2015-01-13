import numpy as np
import pylab as plt
import asciitable
import readsnapshots.readsnapHDF5_greg as rsg
import haloutils
import sys

import profilefit
from scipy.optimize import curve_fit

from caterpillaranalysis import *

def test_fit(rbin,dr,mpart):
    rmid = 10**((np.log10(rbin[1:])+np.log10(rbin[:-1]))/2.)
    #Varr = 4*np.pi/3 * (rbin[1:]-rbin[:-1])**3 #kpc^3
    #h,x = np.histogram(dr,bins=rbin)
    #Marr = h*mpart #Msun
    #rhoarr = Marr/Varr # Msun/kpc^3
    #logrhoarr = np.log10(rhoarr)
    rhoarr = profilefit.calc_rhoarr(rbin,dr,mpart)
    
    p0NFW = [.5,10]
    p0EIN = [.5,10,.18]

    pNFW0,pNFW1,Q2NFW = profilefit.fitNFW(rbin,rhoarr,p0NFW,retQ2=True)
    pEIN0,pEIN1,pEIN2,Q2EIN = profilefit.fitEIN(rbin,rhoarr,p0EIN,retQ2=True)
    pNFW = [pNFW0,pNFW1]
    pEIN = [pEIN0,pEIN1,pEIN2]

    print pNFW,Q2NFW
    print pEIN,Q2EIN
    #pNFW = curve_fit(profilefit.logNFWprofile,rmid,logrhoarr,p0=p0NFW)[0]
    #pEIN = curve_fit(profilefit.logEINprofile,rmid,logrhoarr,p0=p0EIN)[0]

    #print pNFW,pEIN
    #def Q2(y1,y2):
    #    assert len(y1)==len(y2)
    #    Nbin = len(y1)
    #    return np.sum((y1-y2)**2)/Nbin
    #print Q2(logrhoarr,profilefit.logNFWprofile(rmid,pNFW[0],pNFW[1]))
    #print Q2(logrhoarr,profilefit.logEINprofile(rmid,pEIN[0],pEIN[1],pEIN[2]))

    plt.figure()
    plt.plot(rmid,rmid**2*rhoarr,'sk')
    plt.plot(rmid,rmid**2*profilefit.NFWprofile(rmid,pNFW[0],pNFW[1]),'b:')
    plt.plot(rmid,rmid**2*profilefit.EINprofile(rmid,pEIN[0],pEIN[1],pEIN[2]),'r--')
    plt.loglog()
    plt.show()

class SubProfileSoftPlugin(ProfilePlugin):
    def __init__(self,rmin=10**-2,rmax=10**3,ymin=10**-1.5,ymax=10**2.5):
        self.filename='subprofilesoft.npz'
        self.nr = 50
        self.nrfit = 20
        self.rminfit = .291  #kpc, Draco rvmax
        self.mmin = 10**8

        self.xmin = rmin; self.xmax = rmax
        self.ymin = ymin;  self.ymax = ymax
        self.xlabel = r'$r\ (kpc)$' #$h^{-1}$ 
        self.ylabel = r'$r^2 \rho(r)\ [10^{10}\ M_\odot\ Mpc^{-1}]$'
        self.xlog = True; self.ylog = True
        self.autofigname = 'subprofsoft'

    def _analyze(self,hpath):
        snap = 255
        zoomid = haloutils.load_zoomid(hpath)
        rscat = haloutils.load_bound_rscat(hpath,snap)
        subs = rscat.get_all_subhalos_within_halo(zoomid)
        subs = subs[subs['mgrav']/rscat.h0 > self.mmin]
        subids = np.array(subs['id'])
        nsubs = len(subs)
        nr = self.nr
        
        idarr = np.zeros(nsubs)
        rvirarr = np.zeros(nsubs)
        rvmaxarr = np.zeros(nsubs)
        mgravarr = np.zeros(nsubs)
        allmltrarr  = np.zeros((nsubs,nr))
        allmltrsoftarr  = np.zeros((nsubs,nr))
        Q2arr = np.zeros(nsubs)

        snapstr = str(snap).zfill(3)
        snapfile = hpath+'/outputs/snapdir_'+snapstr+'/snap_'+snapstr
        header = rsg.snapshot_header(snapfile+'.0')
        mpart = header.massarr[1]*1e10/header.hubble

        for i,subid in enumerate(subids):
            i_rvir = np.array(subs.ix[subid]['rvir']/rscat.h0) #kpc
            i_rvmax= np.array(subs.ix[subid]['rvmax']/rscat.h0) #kpc
            i_mgrav= np.array(subs.ix[subid]['mgrav']/rscat.h0) #Msun
            idarr[i] = subid
            rvirarr[i] = i_rvir
            rvmaxarr[i]= i_rvmax
            mgravarr[i]= i_mgrav

            rarr = self.get_scaled_rarr(i_rvir) #Mpc
            rarr,mltr,p03rmin,halorvir,r200c,halomass,dr = self.compute_one_profile(rarr,hpath,rscat,subid,snap,header,calcp03r=True,calcr200=True,retdr=True,usebound=True)
            dr   *= 1000. #kpc
            rarr *= 1000. #kpc
            allmltrarr[i,:] = mltr
            if i_rvmax >= .5:
                EINmltr,Q2 = self.compute_mltr_soft(dr,i_rvir,i_rvmax,mpart) 
                if Q2==None: Q2=-1
                elif Q2 < .1:
                    ii = (rarr < self.rminfit)
                    mltr[ii] = EINmltr(rarr[ii])
            else: Q2 = -1
            Q2arr[i] = Q2
            allmltrsoftarr[i,:] = mltr
        np.savez(self.get_outfname(hpath),rsid=idarr,rvir=rvirarr,rvmax=rvmaxarr,
                 mgrav=mgravarr,mltr=allmltrarr,mltrsoft=allmltrsoftarr,Q2=Q2arr)
            
    def compute_mltr_soft(self,dr,rvir,rvmax,mpart):
        rbin = self.get_fit_rarr(rvmax) #kpc
        rhoarr = profilefit.calc_rhoarr(rbin,dr,mpart) #Msun/kpc^3
        try:
            p0,p1,p2,Q2 = profilefit.fitEIN(rbin,rhoarr,[.5,10,.2],retQ2=True)
            EINmltr = lambda r: profilefit.EINmltr(r,p0,p1,p2)
        except RuntimeError as e:
            # if 'maxfev' in error msg
            EINmltr=None; Q1=-1
            # else raise e
        return EINmltr,Q2

    def _read(self,hpath):
        thisfilename = self.get_filename(hpath)
        d = np.load(thisfilename) #d['rsid']
        rarr = self.get_scaled_rarr(d['rvir'])
        return rarr,d['rsid'],d['rvir'],d['rvmax'],d['mgrav'],d['mltr'],d['mltrsoft'],d['Q2']
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        rarr, rsid, rvir, rvmax, mgrav, mltr, mltrsoft, Q2 = data
        mltr /= 1e10; mltrsoft /= 1e10
        rhoarr = np.zeros(mltr.shape)
        rhosoftarr = np.zeros(mltrsoft.shape)
        for i in range(len(rsid)):
            rhoarr[i,:]     = self.mltr_to_rho(rarr[i],mltr[i,:]) #Msun/Mpc^3
            rhosoftarr[i,:] = self.mltr_to_rho(rarr[i],mltrsoft[i,:])
        plotqty = (rarr)**2 * rhoarr #Msun/Mpc
        plotqtysoft = (rarr)**2 * rhosoftarr #Msun/Mpc
        rarr = rarr*1000 #kpc

        eps = 1000*haloutils.load_soft(hpath)
        if lx != None:
            color = self.colordict[lx]
        for i in xrange(len(rsid)):
            if i==5: break
            ii = rarr[i,:] >= eps
            if np.sum(ii) == 0: continue
            ax.plot(rarr[i,ii], plotqty[i,ii], '-', color=color, **kwargs)
            ax.plot(rarr[i,ii], plotqtysoft[i,ii], ':', color=color, **kwargs)

    def get_fit_rarr(self,rvmax):
        rlo = self.rminfit
        rhi = min(3,1.5*rvmax) #kpc
        return np.logspace(np.log10(rlo),np.log10(rhi),self.nrfit)
        
    def get_scaled_rarr(self,rvir):
        """ 50 logspaced bins (3e-5 to 3)*rvir. Input rvir in kpc, return Mpc """
        out = 3*rvir.reshape(-1,1)/1000.*np.logspace(-5,0,self.nr).reshape(1,-1)
        if out.shape[0]==1: return out[0]
        return out

