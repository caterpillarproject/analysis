import matplotlib; matplotlib.use('Agg')
import numpy as np
import pylab as plt
import os,subprocess,sys,time
import cPickle as pickle
from scipy import special
from scipy.spatial.distance import pdist
import pandas as pd

import haloutils
import rotations
from caterpillaranalysis import PluginBase
from SAMs import SimpleSAMBasePlugin

from seaborn.apionly import cubehelix_palette
chcmap = cubehelix_palette(as_cmap=True,start=.5,rot=-1.5,hue=1.0,gamma=1.0)

def plot_mollweide_L(logMpeakcut=None,lx=14,tag='A'):
    hids = haloutils.cid2hid.values()
    for hid in hids:
        hpath = haloutils.get_hpath_lx(hid,lx)
        if hpath==None: continue
        if not haloutils.check_last_rockstar_exists(hpath): continue
        print haloutils.hidstr(hid)
        fig,ax = plt.subplots(subplot_kw={'projection':'mollweide'})
        sc = _plot_mollweide_SAM('angmom',hpath,ax,tag,logMpeakcut=logMpeakcut)
        fig.colorbar(sc,orientation='horizontal')
        if logMpeakcut != None:
            figfilename = '5-19/mollweideL'+tag+'_Mpeak'+str(logMpeakcut)+'_'+haloutils.hidstr(hid)+'.png'
        else:
            figfilename = '5-19/mollweideL'+tag+'_'+haloutils.hidstr(hid)+'.png'

        fig.savefig(figfilename,bbox_inches='tight')
        plt.close('all')

def plot_mollweide_infall(logMpeakcut=None,lx=14,tag='A'):
    hids = haloutils.cid2hid.values()
    for hid in hids:
        hpath = haloutils.get_hpath_lx(hid,lx)
        if hpath==None: continue
        if not haloutils.check_last_rockstar_exists(hpath): continue
        print haloutils.hidstr(hid)
        fig,ax = plt.subplots(subplot_kw={'projection':'mollweide'})
        sc = _plot_mollweide_SAM('infallpos',hpath,ax,tag,logMpeakcut=logMpeakcut)
        fig.colorbar(sc,orientation='horizontal')
        if logMpeakcut != None:
            figfilename = '5-19/mollweideIn'+tag+'_Mpeak'+str(logMpeakcut)+'_'+haloutils.hidstr(hid)+'.png'
        else:
            figfilename = '5-19/mollweideIn'+tag+'_'+haloutils.hidstr(hid)+'.png'

        fig.savefig(figfilename,bbox_inches='tight')
        plt.close('all')

def _plot_mollweide_SAM(whatdata,hpath,ax,tag,logMpeakcut=None):
    # TODO assert ax is mollweide projection
    assert whatdata in ['infallpos','angmom']

    plug = SimpleSAMBasePlugin()
    subs = plug.read(hpath)
    
    rscat = haloutils.load_rscat(hpath,haloutils.get_numsnaps(hpath)-1)
    zoomid= haloutils.load_zoomid(hpath)
    hpos = np.array(rscat.ix[zoomid][['posX','posY','posZ']])
    hvel = np.array(rscat.ix[zoomid][['pecVX','pecVY','pecVZ']])
    hLmom = np.array(rscat.ix[zoomid][['Jx','Jy','Jz']])
    hA = np.array(rscat.ix[zoomid][['A2[x]','A2[y]','A2[z]']])
    
    tagdict = {'Z':np.array([0,0,1]),
               'J':hLmom,
               'A':hA}
    Ldisk = tagdict[tag]
    rotmat = rotations.rotate_to_z(Ldisk)
    
    subs.sort(columns='infall_vmax')
    if whatdata=='angmom':
        spos = np.array(subs[['posX','posY','posZ']])-hpos
        svel = np.array(subs[['pecVX','pecVY','pecVZ']])-hvel
        sLmom = np.cross(spos,svel)
        plot_pos = rotmat.dot(sLmom.T).T
    elif whatdata=='infallpos':
        mtc = haloutils.load_zoom_mtc(hpath,indexbyrsid=True)
        hostmb = mtc[zoomid].getMainBranch()
        mbhostpos = hostmb[['posX','posY','posZ']][::-1].view((np.float,3))
        mbhostsnap= hostmb['snap'][::-1]
        maxsnap = mbhostsnap[-1]; minsnap = mbhostsnap[0]
        assert len(mbhostsnap)==maxsnap-minsnap+1
        assert np.all(mbhostsnap==np.sort(mbhostsnap))
        iigood = ~np.array(np.isnan(subs['infall_snap']))
        infall_ix  = np.zeros(len(subs)).astype(int)
        infall_hpos= np.zeros((len(subs),3))+np.nan
        infall_ix[iigood] = subs['infall_snap'][iigood]-minsnap
        infall_hpos[iigood] = mbhostpos[infall_ix[iigood],:]
        infall_pos = subs[['infall_posx','infall_posy','infall_posz']]-infall_hpos
        plot_pos = rotmat.dot(infall_pos.T).T

    min_vmax_size=0.
    max_vmax_size=100.
    min_vmax=0.
    max_vmax=100.
    normed_vmax = np.array((subs['infall_vmax']-min_vmax)/(max_vmax-min_vmax))
    sizes_vmax = normed_vmax * (max_vmax_size - min_vmax_size) + min_vmax_size

    plot_size= sizes_vmax

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    theta,phi = rotations.xyz2thetaphi(plot_pos,rotate_for_mollweide=True)
    ii = np.array(~np.isnan(subs['infall_scale']))
    if logMpeakcut != None:
        logMpeak = np.log10(np.array(subs['peak_mvir']))
        ii = ii & (logMpeak > logMpeakcut)

    theta = theta[ii]; phi = phi[ii]; infall_scale = np.array(subs['infall_scale'])[ii]
    plot_size = plot_size[ii]
    
    theta = theta[::-1]; phi = phi[::-1]; infall_scale = infall_scale[::-1]; plot_size = plot_size[::-1]
    sc = ax.scatter(phi,theta,c=infall_scale,vmin=0,vmax=1,
                    s=plot_size,linewidth=0,cmap=chcmap)
    return sc

class AngMomCorrelationPlugin(PluginBase):
    def __init__(self):
        super(AngMomCorrelationPlugin,self).__init__()
        self.filename='angmom_corr.p'

        self.xmin = -1; self.xmax = 1
        self.ymin = -0.2; self.ymax = 0.3
        self.xlabel = r'$\cos \theta$'
        self.ylabel = r'$w(\theta)$'
        self.xlog=False; self.ylog=False
        self.autofigname='angmom_corr'

        self.nbins=40; self.lmax=20
        self.samplug = SimpleSAMBasePlugin()
        self.logMpeakcutarr = [6,7,8,9]

    def get_bins(self):
        return np.linspace(-1,1,self.nbins+1)
        
    def angular_correlation(self,pos):
        numpoints = len(pos)
        cosd = 1-pdist(pos,'cosine')
        bins = self.get_bins()
        dd,x = np.histogram(cosd,bins=bins)
        x = (x[:-1]+x[1:])/2.
        rr = float(len(cosd))/len(dd) #uniform in cos(theta)
        w = dd/rr - 1
        return w

    def a_lm(self,l,m,theta,phi):
        #Note special.sph_harm documentation has opposite definitions of theta/phi
        return 1./(np.pi * 4.0) * np.sum(np.conj(special.sph_harm(m,l,phi,theta)))
    def C_l(self,l,theta,phi,rnorm=1):
        alm_arr = np.zeros(2*l+1,dtype=complex)
        for i in range(2*l+1):
            m = i-l
            alm_arr[i] = self.a_lm(l,m,theta,phi)
        sumsquares = np.sum(alm_arr * np.conj(alm_arr))
        return sumsquares/(4*np.pi*(2*l+1)*rnorm**2)
    def angular_powspec(self,pos,lmax,lmin=0):
        Cl_arr = np.zeros(lmax-lmin+1,dtype=complex)
        theta,phi = rotations.xyz2thetaphi(pos)
        for i,l in enumerate(np.arange(lmin,lmax+1)):
            Cl_arr[i] = self.C_l(l,theta,phi)
        return Cl_arr

    def angular_powspec_from_corr(self,w,lmax):
        bins = self.get_bins()
        x = (bins[1:]+bins[:-1])/2.
        dx = x[1]-x[0]
        Cl_arr = np.zeros(lmax+1)
        for l in range(lmax+1):
            # 2pi * integral(w(x) P_l(x) dx, -1, 1)
            Cl_arr[l] = 2*np.pi * np.sum(2*np.pi*w*special.legendre(l)(x)) * dx
        return Cl_arr

    def _analyze(self,hpath):
        subs = self.samplug.read(hpath)
        try:
            badsubs = subs==None
        except:
            badsubs = False
        if badsubs:
            raise IOError('No SimpleSAMs.p')
        rscat = haloutils.load_rscat(hpath,haloutils.get_numsnaps(hpath)-1)
        zoomid= haloutils.load_zoomid(hpath)
        hpos = np.array(rscat.ix[zoomid][['posX','posY','posZ']])
        hvel = np.array(rscat.ix[zoomid][['pecVX','pecVY','pecVZ']])
        hLmom = np.array(rscat.ix[zoomid][['Jx','Jy','Jz']])
        hA = np.array(rscat.ix[zoomid][['A2[x]','A2[y]','A2[z]']])
        
        spos = np.array(subs[['posX','posY','posZ']])-hpos
        svel = np.array(subs[['pecVX','pecVY','pecVZ']])-hvel
        sLmom = np.cross(spos,svel)
        logMpeak = np.log10(np.array(subs['peak_mvir']))

        bins = self.get_bins()

        wlist = []
        Cllist = []
        for logMpeakcut in self.logMpeakcutarr:
            ii = logMpeak > logMpeakcut
            w = self.angular_correlation(sLmom[ii,:])
            wlist.append(w)

            #Cl = self.angular_powspec_from_corr(w,self.nbins)
            Cl = self.angular_powspec(sLmom[ii,:],self.lmax)
            Cllist.append(Cl)

        with open(self.get_outfname(hpath),'w') as f:
            pickle.dump([bins,wlist,Cllist,self.logMpeakcutarr],f)

    def _read(self,hpath):
        try:
            with open(self.get_outfname(hpath),'r') as f:
                bins,wlist,Cllist,logMpeakcutarr = pickle.load(f)
            return bins,wlist,logMpeakcutarr
        except IOError:
            return None

    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        bins,wlist,logMpeakcutarr = data
        x = (bins[1:]+bins[:-1])/2.
        if lx != None:
            for w in wlist:
                ax.plot(x,w,color=self.colordict[lx],**kwargs)
        else:
            for w in wlist:
                ax.plot(x,w,**kwargs)

class AngMomPowerSpectrumPlugin(AngMomCorrelationPlugin):
    def __init__(self):
        super(AngMomPowerSpectrumPlugin,self).__init__()

        self.xmin = 0; self.xmax = self.lmax
        self.ymin = 1e-4; self.ymax = 10**0.2
        self.xlabel = r'$\ell$'
        self.ylabel = r'$C_\ell/C_0$'
        self.xlog=False; self.ylog=True
        self.autofigname='angmom_powspec'

    def _read(self,hpath):
        try:
            with open(self.get_outfname(hpath),'r') as f:
                bins,wlist,Cllist,logMpeakcutarr = pickle.load(f)
            return Cllist,logMpeakcutarr
        except IOError:
            return None

    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,maxlines=None,**kwargs):
        Cllist,logMpeakcutarr = data
        l = np.arange(len(Cllist[0]))
        numlines=0
        if lx != None:
            for Cl in Cllist:
                ax.plot(l,Cl/Cl[0],'s-',color=self.colordict[lx],**kwargs)
                numlines += 1
                if (maxlines != None) and (numlines == maxlines): break
        else:
            for Cl in Cllist:
                ax.plot(l,Cl/Cl[0],'s-',**kwargs)
                numlines += 1
                if (maxlines != None) and (numlines == maxlines): break

