import numpy as np
import pylab as plt
import asciitable
import readsnapshots.readsnapHDF5_greg as rsg
import haloutils
import sys

import profilefit

from caterpillaranalysis import *

class SubProfileSoftPlugin(ProfilePlugin):
    def __init__(self,rmin=10**-2,rmax=10**3,ymin=10**-1.5,ymax=10**1.5):
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
        rscat = haloutils.load_rscat(hpath,snap)
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
            rarr,mltr,p03rmin,halorvir,r200c,halomass,dr = self.compute_one_profile(rarr,hpath,rscat,subid,snap,header,calcp03r=True,calcr200=True,retdr=True)
            dr   *= 1000. #kpc
            rarr *= 1000. #kpc
            allmltrarr[i,:] = mltr
            Marr = self.mltr_to_Marr(mltr)
            if i_rvmax >= .5:
                EINmltr,Q2 = self.compute_mltr_soft(dr,i_rvir,i_rvmax,mpart) 
                if Q2==None: Q2=-1
                elif Q2 < .1:
                    ii = (rarr < self.rminfit)
                    mltr[ii] = EINmltr(rarr[ii])
                    iilast = np.max(np.where(ii)[0])
                    mltrlast = mltr[iilast]
                    mltr[~ii] = np.cumsum(Marr[~ii])+mltrlast
            else: Q2 = -1
            Q2arr[i] = Q2
            allmltrsoftarr[i,:] = mltr
        np.savez(self.get_outfname(hpath),rsid=idarr,rvir=rvirarr,rvmax=rvmaxarr,
                 mgrav=mgravarr,mltr=allmltrarr,mltrsoft=allmltrsoftarr,Q2=Q2arr)
            
    def compute_rho_soft(self,dr,rvir,rvmax,mpart):
        rbin = self.get_fit_rarr(rvmax) #kpc
        rhoarr = profilefit.calc_rhoarr(rbin,dr,mpart) #Msun/kpc^3
        try:
            p0,p1,p2,Q2 = profilefit.fitEIN(rbin,rhoarr,[.5,10,.2],retQ2=True)
            EINprof = lambda r: profilefit.EINprofile(r,p0,p1,p2)
        except RuntimeError as e:
            # if 'maxfev' in error msg
            EINprof=None; Q1=-1
            # else raise e
        return EINprof,Q2
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
            if i==10: break
            ii = rarr[i,:] >= eps
            if np.sum(ii) == 0: continue
            ax.plot(rarr[i,ii], plotqty[i,ii], ':', color=color, **kwargs)
            ax.plot(rarr[i,ii], plotqtysoft[i,ii], '-', color=color, **kwargs)

    def get_fit_rarr(self,rvmax):
        rlo = self.rminfit
        rhi = min(3,1.5*rvmax) #kpc
        return np.logspace(np.log10(rlo),np.log10(rhi),self.nrfit)
        
    def get_scaled_rarr(self,rvir):
        """ 50 logspaced bins (3e-5 to 3)*rvir. Input rvir in kpc, return Mpc """
        out = 3*rvir.reshape(-1,1)/1000.*np.logspace(-5,0,self.nr).reshape(1,-1)
        if out.shape[0]==1: return out[0]
        return out


class SubVelocityProfileSoftPlugin(SubProfileSoftPlugin):
    def __init__(self,rmin=10**-1.5,rmax=10**3,vmin=10**0,vmax=10**2.3):
        super(SubVelocityProfileSoftPlugin,self).__init__()
        self.xmin = rmin; self.xmax = rmax
        self.ymin = vmin;  self.ymax = vmax
        self.xlabel = r'$r\ (kpc)$'
        self.ylabel = r'$v_{circ}\ (km/s)$'
        self.n_xmin = 10**-3.9; self.n_xmax = 10**0.5
        self.n_ymin = 10**-2.0; self.n_ymax = 10**0.2
        self.n_xlabel = r'$r/r_{\rm vir,host}$'
        self.n_ylabel = r'$v_{\rm circ}/v_{\rm vir,host}$'
        self.xlog = True; self.ylog = True
        self.autofigname = 'subvcircsoft'
    def _read(self,hpath):
        thisfilename = self.get_filename(hpath)
        d = np.load(thisfilename) #d['rsid']
        rsid = d['rsid']
        rarr = self.get_scaled_rarr(d['rvir'])
        rvir = d['rvir']
        mltrarr = d['mltr']; mltrsoftarr = d['mltrsoft']
        vcircarr = np.zeros(mltrarr.shape)
        vcircsoftarr = np.zeros(mltrsoftarr.shape)
        for i in range(len(rsid)):
            vcircarr[i,:] = self.mltr_to_vcirc(rarr[i],mltrarr[i,:])
            vcircsoftarr[i,:] = self.mltr_to_vcirc(rarr[i],mltrsoftarr[i,:])
        return rsid,rarr,rvir,vcircarr,vcircsoftarr
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,alpha=.2,color='k',**kwargs):
        rsid,rarr,rvir,vcircarr,vcircsoftarr = data
        rarr = rarr*1000 #kpc
        eps = 1000*haloutils.load_soft(hpath)
        if normtohost:
            mvir,rvir,vvir=haloutils.load_haloprops(hpath)
            rarr = rarr/rvir
            vcircarr = vcircarr/vvir
            vcircsoftarr = vcircsoftarr/vvir
            eps = eps/rvir
        if lx != None:
            color = self.colordict[lx]
        for i in xrange(len(rsid)):
            ii = rarr[i,:] >= eps
            if np.sum(ii) == 0: continue
            ax.plot(rarr[i,ii], vcircarr[i,ii], color=color, alpha=alpha, **kwargs)
