import matplotlib; matplotlib.use('Agg')
import numpy as np
import pylab as plt
import os,subprocess,sys,time,warnings
import cPickle as pickle
from scipy import special
from scipy.spatial.distance import pdist
import pandas as pd

import haloutils
import rotations
from caterpillaranalysis import PluginBase,MassAccrPlugin
from SAMs import SimpleSAMBasePlugin
#from fast_SAMs import FastSAMBasePlugin as SimpleSAMBasePlugin

from seaborn.apionly import cubehelix_palette,color_palette
chcmap = cubehelix_palette(as_cmap=True,start=.5,rot=-1.5,hue=1.0,gamma=1.0)
default_colors = color_palette('muted')

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
            figfilename = '5-29/mollweideL/mollweideL'+tag+'_Mpeak'+str(logMpeakcut)+'_'+haloutils.hidstr(hid)+'.png'
        else:
            figfilename = '5-29/mollweideL/mollweideL'+tag+'_'+haloutils.hidstr(hid)+'.png'

        fig.savefig(figfilename,bbox_inches='tight')
        plt.close('all')

def plot_mollweide_pos(logMpeakcut=None,lx=14,tag='A'):
    hids = haloutils.cid2hid.values()
    for hid in hids:
        hpath = haloutils.get_hpath_lx(hid,lx)
        if hpath==None: continue
        if not haloutils.check_last_rockstar_exists(hpath): continue
        print haloutils.hidstr(hid)
        fig,ax = plt.subplots(subplot_kw={'projection':'mollweide'})
        sc = _plot_mollweide_SAM('satpos',hpath,ax,tag,logMpeakcut=logMpeakcut)
        fig.colorbar(sc,orientation='horizontal')
        if logMpeakcut != None:
            figfilename = '5-29/mollweideX/mollweideX'+tag+'_Mpeak'+str(logMpeakcut)+'_'+haloutils.hidstr(hid)+'.png'
        else:
            figfilename = '5-29/mollweideX/mollweideX'+tag+'_'+haloutils.hidstr(hid)+'.png'

        fig.savefig(figfilename,bbox_inches='tight')
        plt.close('all')

def plot_mollweide_infall(logMpeakcut=None,lx=14,tag='A'):
    raise NotImplementedError
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

def plot_mollweide_host(lx=14):
    hids = haloutils.cid2hid.values()
    plug = MassAccrPlugin()
    for hid in hids:
        hpath = haloutils.get_hpath_lx(hid,lx)
        if hpath==None: continue
        if not haloutils.check_last_rockstar_exists(hpath): continue
        print haloutils.hidstr(hid)

        hostmb = plug.read(hpath)
        hostA   = hostmb[['A[x]','A[y]','A[z]']].view((np.float,3))
        hostJ   = hostmb[['Jx','Jy','Jz']].view((np.float,3))
        hostscale = hostmb['scale']

        thetaA,phiA = rotations.xyz2thetaphi(hostA,rotate_for_mollweide=True)
        thetaJ,phiJ = rotations.xyz2thetaphi(hostJ,rotate_for_mollweide=True)

        fig,ax = plt.subplots(subplot_kw={'projection':'mollweide'})
        ax.grid(True)
        ax.set_xticklabels(['' for i in ax.get_xticklabels()])
        ax.set_yticklabels(['' for i in ax.get_yticklabels()])
        sc = ax.scatter(phiA,thetaA,c=hostscale,vmin=0,vmax=1,
                        marker='s',linewidth=0,cmap=chcmap)
        fig.colorbar(sc,orientation='horizontal')
        figfilename = '5-20/mollweideHostA'+'_'+haloutils.hidstr(hid)+'.png'
        fig.savefig(figfilename,bbox_inches='tight')

        fig,ax = plt.subplots(subplot_kw={'projection':'mollweide'})
        ax.grid(True)
        ax.set_xticklabels(['' for i in ax.get_xticklabels()])
        ax.set_yticklabels(['' for i in ax.get_yticklabels()])
        sc = ax.scatter(phiJ,thetaJ,c=hostscale,vmin=0,vmax=1,
                        marker='o',linewidth=0,cmap=chcmap)
        fig.colorbar(sc,orientation='horizontal')
        figfilename = '5-20/mollweideHostJ'+'_'+haloutils.hidstr(hid)+'.png'
        fig.savefig(figfilename,bbox_inches='tight')

        fig,ax = plt.subplots(subplot_kw={'projection':'mollweide'})
        ax.grid(True)
        ax.set_xticklabels(['' for i in ax.get_xticklabels()])
        ax.set_yticklabels(['' for i in ax.get_yticklabels()])
        sc = ax.scatter(phiA,thetaA,c=hostscale,vmin=0,vmax=1,
                        marker='s',linewidth=0,cmap=chcmap)
        sc = ax.scatter(phiJ,thetaJ,c=hostscale,vmin=0,vmax=1,
                        marker='o',linewidth=0,cmap=chcmap)
        fig.colorbar(sc,orientation='horizontal')
        figfilename = '5-20/mollweideHostAJ'+'_'+haloutils.hidstr(hid)+'.png'
        fig.savefig(figfilename,bbox_inches='tight')

        plt.close('all')

def plot_mollweide_time(hids=None,lx=14):
    if hids==None:
        hids = haloutils.cid2hid.values()
    plug = MassAccrPlugin()
    for hid in hids:
        hpath = haloutils.get_hpath_lx(hid,lx)
        if hpath==None: continue
        if not haloutils.check_last_rockstar_exists(hpath): continue
        print haloutils.hidstr(hid)

        hostmb = plug.read(hpath)
        hostA   = hostmb[['A[x]','A[y]','A[z]']].view((np.float,3))
        hostJ   = hostmb[['Jx','Jy','Jz']].view((np.float,3))
        hostscale = hostmb['scale']
        hostsnaps = hostmb['snap']
        hostpos   = hostmb[['x','y','z']].view((np.float,3))
        hostrvir  = hostmb['rvir']
        hostzoomid= hostmb['origid']


        figfilenamebase = '5-20/'+haloutils.hidstr(hid)
        subprocess.call(['mkdir -p '+figfilenamebase],shell=True)
        for snap,zoomid,hA,hJ,hpos,hrvir in zip(hostsnaps,hostzoomid,hostA,hostJ,hostpos,hostrvir):
            snapstr = str(snap).zfill(3)

            rscat = haloutils.load_rscat(hpath,snap,rmaxcut=False)
            subs = rscat.get_all_subhalos_within_halo(zoomid)
            if len(subs)==0: continue
            spos = np.array(subs[['posX','posY','posZ']])-hpos
            svmax= np.array(subs['vmax'])
            sdr  = 1000.*np.sqrt(np.sum(spos**2,1))/hrvir

            min_vmax_size=0.
            max_vmax_size=100.
            min_vmax=0.
            max_vmax=100.
            normed_vmax = np.array((svmax-min_vmax)/(max_vmax-min_vmax))
            sizes_vmax = normed_vmax * (max_vmax_size - min_vmax_size) + min_vmax_size

            for tag,direction in zip(['Z','J','A'],[np.array([0,0,1]),hJ,hA]):
                figfilename = figfilenamebase+'/'+tag+snapstr+'.png'
                rotmat = rotations.rotate_to_z(direction)
                thispos = rotmat.dot(spos.T).T
                theta,phi = rotations.xyz2thetaphi(thispos,rotate_for_mollweide=True)

                fig,ax = plt.subplots(subplot_kw={'projection':'mollweide'})
                ax.grid(True)
                ax.set_xticklabels(['' for i in ax.get_xticklabels()])
                ax.set_yticklabels(['' for i in ax.get_yticklabels()])
                sc = ax.scatter(phi,theta,c=sdr,vmin=0,vmax=1,
                                s=sizes_vmax,linewidth=0,cmap=chcmap)
                fig.colorbar(sc,orientation='horizontal')
                fig.savefig(figfilename)
            plt.close('all')
        break

def _plot_mollweide_SAM(whatdata,hpath,ax,tag,logMpeakcut=None):
    # TODO assert ax is mollweide projection
    assert whatdata in ['infallpos','angmom','satpos']

    plug = SimpleSAMBasePlugin()
    subs = plug.read(hpath)
    mbplug = MassAccrPlugin()
    
    mb = mbplug.read(hpath)
    host = mb[-1]
    hpos = np.array([host['x'],host['y'],host['z']])
    hvel = np.array([host['vx'],host['vy'],host['vz']])
    hLmom = np.array([host['Jx'],host['Jy'],host['Jz']])
    hA = np.array([host['A[x]'],host['A[y]'],host['A[z]']])
    #rscat = haloutils.load_rscat(hpath,haloutils.get_numsnaps(hpath)-1)
    #zoomid= haloutils.load_zoomid(hpath)
    #hpos = np.array(rscat.ix[zoomid][['posX','posY','posZ']])
    #hvel = np.array(rscat.ix[zoomid][['pecVX','pecVY','pecVZ']])
    #hLmom = np.array(rscat.ix[zoomid][['Jx','Jy','Jz']])
    #hA = np.array(rscat.ix[zoomid][['A2[x]','A2[y]','A2[z]']])
    
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
    elif whatdata=='satpos':
        spos = np.array(subs[['posX','posY','posZ']])-hpos
        plot_pos = rotmat.dot(spos.T).T
    elif whatdata=='infallpos':
        raise NotImplementedError("rotating to A is not a very stable coordinate frame")
        hostmb = mbplug.read(hpath)
        mbhostsnap= hostmb['snap'][::-1]
        maxsnap = mbhostsnap[-1]; minsnap = mbhostsnap[0]
        assert len(mbhostsnap)==maxsnap-minsnap+1
        assert np.all(mbhostsnap==np.sort(mbhostsnap))
        mbhostpos = hostmb[['x','y','z']][::-1].view((np.float,3))
        mbhostA   = hostmb[['A[x]','A[y]','A[z]']][::-1].view((np.float,3))
        rotmatlist = [rotations.rotate_to_z(np.array(Ldisk)) for Ldisk in mbhostA]
        iigood = ~np.array(np.isnan(subs['infall_snap']))
        infall_ix  = np.zeros(len(subs)).astype(int)
        infall_hpos= np.zeros((len(subs),3))+np.nan
        infall_ix[iigood] = subs['infall_snap'][iigood]-minsnap
        infall_hpos[iigood] = mbhostpos[infall_ix[iigood],:]
        infall_pos = subs[['infall_posx','infall_posy','infall_posz']]-infall_hpos
        for snap,rotmat in zip(range(minsnap,maxsnap+1),rotmatlist):
            iisnap = subs['infall_snap']==snap
            thispos = np.array(infall_pos[iisnap])
            infall_pos[iisnap] = rotmat.dot(thispos.T).T
        plot_pos = np.array(infall_pos)

    min_vmax_size=0.
    max_vmax_size=100.
    min_vmax=0.
    max_vmax=100.
    normed_vmax = np.array((subs['infall_vmax']-min_vmax)/(max_vmax-min_vmax))
    sizes_vmax = normed_vmax * (max_vmax_size - min_vmax_size) + min_vmax_size

    plot_size= sizes_vmax

    ax.grid(True)
    ax.set_xticklabels(['' for i in ax.get_xticklabels()])
    ax.set_yticklabels(['' for i in ax.get_yticklabels()])
    theta,phi = rotations.xyz2thetaphi(plot_pos,rotate_for_mollweide=True)
    ii = np.array(~np.isnan(subs['infall_scale']))
    if logMpeakcut != None:
        #logMpeak = np.log10(np.array(subs['infall_mvir']))
        logMpeak = np.log10(np.array(subs['max_mass']))
        ii = ii & (logMpeak > logMpeakcut)

    theta = theta[ii]; phi = phi[ii]; infall_scale = np.array(subs['infall_scale'])[ii]
    plot_size = plot_size[ii]
    
    theta = theta[::-1]; phi = phi[::-1]; infall_scale = infall_scale[::-1]; plot_size = plot_size[::-1]
    sc = ax.scatter(phi,theta,c=infall_scale,vmin=0,vmax=1,
                    s=plot_size,linewidth=0,cmap=chcmap)
    return sc

class AngMomCorrelationPlugin(PluginBase):
    def __init__(self,enhancethresh=1.06):
        super(AngMomCorrelationPlugin,self).__init__()
        self.filename='angmom_corr.p'

        self.xmin = -1; self.xmax = 1
        self.ymin = -0.2; self.ymax = 0.3
        self.xlabel = r'$\cos\, \theta_{\mathbf{L}}$'
        self.ylabel = r'$w(\theta_{\mathbf{L}})$'
        self.xlog=False; self.ylog=False
        self.autofigname='angmom_corr'

        self.nbins=40; self.lmax=20
        self.samplug = SimpleSAMBasePlugin()
        self.logMpeakcutarr = [6,7,8,9]

        x = self.get_bins()
        self.x = (x[1:]+x[:-1])/2.

        self.enhancethresh = enhancethresh

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

    def correlation_endfraction(self,w,theta_deg):
        ii = np.abs(self.x) > np.cos(theta_deg*np.pi/180.)
        return np.sum(1.+w[ii])/np.sum(1.+w)
    def correlation_enhancement(self,w,theta_deg):
        ii = np.abs(self.x) > np.cos(theta_deg*np.pi/180.)
        return np.sum(1.+w[ii])/np.sum(ii)
        
    def powspec_ratio(self,Cl):
        return Cl[2]/np.sum(Cl)

    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
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
        #logMpeak = np.log10(np.array(subs['infall_mvir']))
        logMpeak = np.log10(np.array(subs['max_mass']))

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

    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,maxlines=1,
              color_by_enhanced=False,
              color1=default_colors[0],color2=default_colors[1],**kwargs):
        bins,wlist,logMpeakcutarr = data
        x = (bins[1:]+bins[:-1])/2.
        numlines = 0
        hid = haloutils.get_parent_hid(hpath)
        for w in wlist:
            enhancement = self.correlation_enhancement(w,45.)
            if color_by_enhanced:
                if enhancement > self.enhancethresh: color = color1
                else: color = color2
            else:
                if 'color' in kwargs: 
                    color = kwargs.pop('color')
                else: color = None
            ax.plot(x,w,color=color,**kwargs)
            numlines += 1
            if numlines >= maxlines: break

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

class CorrelationPowSpecSnapsPlugin(AngMomCorrelationPlugin):
    def __init__(self,verbose=False):
        super(CorrelationPowSpecSnapsPlugin,self).__init__()
        self.filename = 'corr_powspec_snaps.p'

        self.xmin = None; self.xmax = None
        self.ymin = None; self.ymax = None
        self.xlabel = None
        self.ylabel = None
        self.xlog = False
        self.ylog = False
        self.autofigname='corr_snaps'

        self.vmaxcutarr = [0,5,10,15,20]
        self.mbplug = MassAccrPlugin()
        self.lmax = 60

        self.verbose = verbose

    def _analyze(self,hpath):
        numsnaps = haloutils.get_numsnaps(hpath)
        hostmb = self.mbplug.read(hpath)
        hostmb = self.mbplug.get_phantomless_mb(hostmb)

        hostsnaps = hostmb['snap']
        hostpos = hostmb[['x','y','z']].view((np.float,3))
        hostvel = hostmb[['vx','vy','vz']].view((np.float,3))
        hostrsid = hostmb['origid']

        alldata = []
        if self.verbose: print haloutils.hidstr(haloutils.get_parent_hid(hpath))
        for snap,hpos,hvel,zoomid in zip(hostsnaps,hostpos,hostvel,hostrsid):
            if self.verbose: start = time.time(); print snap
            rscat = haloutils.load_rscat(hpath,snap)
            subs = rscat.get_all_subhalos_within_halo(zoomid)
            svmax= np.array(subs['vmax'])
            spos = np.array(subs[['posX','posY','posZ']])-hpos
            svel = np.array(subs[['pecVX','pecVY','pecVZ']])-hvel
            sLmom = np.cross(spos,svel)
            if len(subs)==0: continue
            snapdata_Lw = []
            snapdata_LCl = []
            snapdata_Xw = []
            snapdata_XCl = []
            for vmaxcut in self.vmaxcutarr:
                iicut = svmax > vmaxcut
                Lw = self.angular_correlation(sLmom[iicut,:])
                LCl= self.angular_powspec(sLmom[iicut,:],self.lmax)
                Xw = self.angular_correlation(spos[iicut,:])
                XCl= self.angular_powspec(spos[iicut,:],self.lmax)
                snapdata_Lw.append(Lw)
                snapdata_LCl.append(LCl)
                snapdata_Xw.append(Xw)
                snapdata_XCl.append(XCl)
            alldata.append([snapdata_Lw,snapdata_LCl,snapdata_Xw,snapdata_XCl])
            if self.verbose:
                print " snap {0}: {1:.1f} sec".format(snap,time.time()-start)
                sys.stdout.flush()

        bins = self.get_bins(); lmax = self.lmax
        with open(self.get_outfname(hpath),'w') as f:
            pickle.dump([bins,lmax,self.vmaxcutarr,hostsnaps,alldata],f)

    def _read(self,hpath):
        try:
            with open(self.get_outfname(hpath),'r') as f:
                bins,lmax,vmaxcutarr,hostsnaps,alldata = pickle.load(f)
            assert lmax==self.lmax
            assert len(self.vmaxcutarr) == len(alldata[0][0])
        except IOError:
            return None
        return alldata

    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,maxlines=None,**kwargs):
        raise NotImplementedError


class AngMomCorrelationSnapsPlugin(CorrelationPowSpecSnapsPlugin):
    def __init__(self,**kwargs):
        super(AngMomCorrelationSnapsPlugin,self).__init__(**kwargs)
        self.xmin = 0; self.xmax = 1
        self.ymin = 0.9; self.ymax = 1.25
        self.xlog = False; self.ylog = False
        self.xlabel = 'scale'
        self.ylabel = '45 deg enhancement'
    def _read(self,hpath):
        try:
            with open(self.get_outfname(hpath),'r') as f:
                bins,lmax,vmaxcutarr,hostsnaps,alldata = pickle.load(f)
            assert lmax==self.lmax
            assert len(self.vmaxcutarr) == len(alldata[0][0])
            angmomcorr = [snapdata[0] for snapdata in alldata]
        except IOError:
            return None
        return angmomcorr,bins
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,whichlines=None,**kwargs):
        angmomcorr,bins = data
        numvmaxcut = len(self.vmaxcutarr)
        hostmb = self.mbplug.read(hpath)
        hostmb = self.mbplug.get_phantomless_mb(hostmb)
        hostscale = hostmb['scale']
        output = np.zeros((numvmaxcut,len(hostscale))) + np.nan
        for i,vmaxcut in enumerate(self.vmaxcutarr):
            for j,wlist in enumerate(angmomcorr):
                w = wlist[i]
                output[i,j] = self.correlation_enhancement(w,45.)
        if lx != None:
            for i in range(numvmaxcut):
                if whichlines != None and i not in whichlines: continue
                yplot = output[i,:]
                ax.plot(hostscale,yplot,color=self.colordict[lx],**kwargs)
        else:
            for i in range(numvmaxcut):
                if whichlines != None and i not in whichlines: continue
                yplot = output[i,:]
                ax.plot(hostscale,yplot,**kwargs)

class AngMomPowSpecSnapsPlugin(CorrelationPowSpecSnapsPlugin):
    def __init__(self,**kwargs):
        super(AngMomPowSpecSnapsPlugin,self).__init__(**kwargs)
        self.xmin = 0; self.xmax = 1
        self.ymin = 1e-3; self.ymax = 1
        self.xlog = False; self.ylog = True
        self.xlabel = 'scale'
        self.ylabel = r'$C_2/\sum_\ell C_\ell$'
    def _read(self,hpath):
        try:
            with open(self.get_outfname(hpath),'r') as f:
                bins,lmax,vmaxcutarr,hostsnaps,alldata = pickle.load(f)
            assert lmax==self.lmax
            assert len(self.vmaxcutarr) == len(alldata[0][0])
            angmompowspec = [snapdata[1] for snapdata in alldata]
        except IOError:
            return None
        return angmompowspec,lmax
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,whichlines=None,**kwargs):
        angmompowspec,lmax = data
        numvmaxcut = len(self.vmaxcutarr)
        hostmb = self.mbplug.read(hpath)
        hostmb = self.mbplug.get_phantomless_mb(hostmb)
        hostscale = hostmb['scale']
        output = np.zeros((numvmaxcut,len(hostscale))) + np.nan
        for i,vmaxcut in enumerate(self.vmaxcutarr):
            for j,Cllist in enumerate(angmompowspec):
                Cl = Cllist[i]
                output[i,j] = self.powspec_ratio(Cl)
        if lx != None:
            for i in range(numvmaxcut):
                if whichlines != None and i not in whichlines: continue
                yplot = output[i,:]
                ax.plot(hostscale,yplot,color=self.colordict[lx],**kwargs)
        else:
            for i in range(numvmaxcut):
                if whichlines != None and i not in whichlines: continue
                yplot = output[i,:]
                ax.plot(hostscale,yplot,**kwargs)
        ax.plot([0,1],[1,1],'k:')

class InPosCorrelationSnapsPlugin(CorrelationPowSpecSnapsPlugin):
    def __init__(self,**kwargs):
        super(InPosPowSpecSnapsPlugin,self).__init__(**kwargs)
    def _read(self,hpath):
        try:
            with open(self.get_outfname(hpath),'r') as f:
                bins,lmax,vmaxcutarr,hostsnaps,alldata = pickle.load(f)
            assert lmax==self.lmax
            assert len(self.vmaxcutarr) == len(alldata[0][0])
            inposcorr = [snapdata[2] for snapdata in alldata]
        except IOError:
            return None
        return inposcorr,bins
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,maxlines=None,**kwargs):
        inposcorr,bins = data
        numvmaxcut = len(self.vmaxcutarr)
        hostmb = self.mbplug.read(hpath)
        hostmb = self.mbplug.get_phantomless_mb(hostmb)
        hostscale = hostmb['scale']
        output = np.zeros((numvmaxcut,len(hostscale))) + np.nan
        for i,vmaxcut in enumerate(self.vmaxcutarr):
            for j,wlist in enumerate(inposcorr):
                w = wlist[i]
                output[i,j] = self.correlation_endfraction(w)
        if lx != None:
            for i in range(numvmaxcut):
                yplot = output[i,:]
                ax.plot(hostscale,yplot,color=self.colordict[lx],**kwargs)
        else:
            for i in range(numvmaxcut):
                yplot = output[i,:]
                ax.plot(hostscale,yplot,**kwargs)

class InPosPowSpecSnapsPlugin(CorrelationPowSpecSnapsPlugin):
    def __init__(self,**kwargs):
        super(InPosPowSpecSnapsPlugin,self).__init__(**kwargs)
    def _read(self,hpath):
        try:
            with open(self.get_outfname(hpath),'r') as f:
                bins,lmax,vmaxcutarr,hostsnaps,alldata = pickle.load(f)
            assert lmax==self.lmax
            assert len(self.vmaxcutarr) == len(alldata[0][0])
            inpospowspec = [snapdata[3] for snapdata in alldata]
        except IOError:
            return None
        return inpospowspec
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,maxlines=None,**kwargs):
        inpospowspec,lmax = data
        numvmaxcut = len(self.vmaxcutarr)
        hostmb = self.mbplug.read(hpath)
        hostmb = self.mbplug.get_phantomless_mb(hostmb)
        hostscale = hostmb['scale']
        output = np.zeros((numvmaxcut,len(hostscale))) + np.nan
        for i,vmaxcut in enumerate(self.vmaxcutarr):
            for j,Cllist in enumerate(inpospowspec):
                Cl = Cllist[i]
                output[i,j] = self.powspec_ratio(Cl)
        if lx != None:
            for i in range(numvmaxcut):
                yplot = output[i,:]
                ax.plot(hostscale,yplot,color=self.colordict[lx],**kwargs)
        else:
            for i in range(numvmaxcut):
                yplot = output[i,:]
                ax.plot(hostscale,yplot,**kwargs)
        ax.plot([0,1],[1,1],'k:')
