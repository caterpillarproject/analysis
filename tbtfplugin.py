import numpy as np
import pylab as plt
import os,subprocess,sys,time
import asciitable
import pynbody as pnb
import readsnapshots.readsnapHDF5_greg as rsg
import haloutils
import pandas as pd

from caterpillaranalysis import *
from subprofileplugin import SubVelocityProfileSoftPlugin

import MTanalysis

class TBTFPlugin(SubVelocityProfilePlugin):
    def __init__(self,minlum=2.e5):
        super(TBTFPlugin,self).__init__()
        self.xmin = .1; self.xmax = 2
        self.ymin = 8;  self.ymax = 50
        self.autofigname = 'tbtf'
        self.minlum = minlum
        self.loadwolf10(minlum=self.minlum)
        self.vcut = 30.
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,alpha=.2,**kwargs):
        assert lx != None
        color = self.colordict[lx]

        #TODO factor out grabbing vpeak?
        extantplug = MTanalysis.TagExtantPlugin()
        extantdata = extantplug.read(hpath,autocalc=False)
        if extantdata == None:
            print "ERROR: %s does not have Greg's Extant data (not plotting)" % (haloutils.get_foldername(hpath))
            return
        extantids,extantdata = extantdata
        vpeakarr = extantdata['vpeak']
        keeprsids = (extantdata['rsid'][vpeakarr >= self.vcut]).astype(int)

        rsid,rarr,rvir,vcircarr = data
        iicut = np.in1d(rsid,keeprsids,assume_unique=True)
        rsid = rsid[iicut]
        rvir = rvir[iicut]
        vcircarr = vcircarr[iicut,:]

        nbound = len(rsid)
        nall = len(keeprsids)
        #if nall != nbound:
        #    print "ERROR %s: bound halos have %i >= %3.1f, but extant has %i" % (haloutils.get_foldername(hpath),
        #                                                                         len(rsid),self.vcut,len(keeprsids))

        rarr = rarr*1000 #kpc
        eps = 1000*haloutils.load_soft(hpath)
        if normtohost:
            mvir,rvir,vvir=haloutils.load_haloprops(hpath)
            rarr = rarr/rvir
            vcircarr = vcircarr/vvir
            eps = eps/rvir
            normtohost=(rvir,vvir)
        self.plotwolf10(ax,normtohost=normtohost)
        for i in xrange(len(rsid)):
            ii = rarr[i,:] >= eps
            if np.sum(ii) == 0: continue
            ax.plot(rarr[i,ii], vcircarr[i,ii], color=color, lw=2, alpha=alpha, **kwargs)
            xmin,xmax,ymin,ymax,xlog,ylog,xlabel,ylabel = self.get_plot_params(normtohost)
            logxoff = np.log10(xmax/xmin)*.05
            xlabel  = xmax * 10**(-logxoff)
            logyoff = np.log10(ymax/ymin)*.1
            ylabel  = ymax * 10**(-logyoff)
            ax.text(xlabel,ylabel,r"$%i/%i$" % (nbound,nall),ha='right')

    def plotwolf10(self,ax,bigerr=False,normtohost=False,**kwargs):
        r = self.rhalf; v = self.vhalf
        xerr = self.rerr
        if bigerr:
            yerr = self.verr2
        else:
            yerr = self.verr
        if normtohost != False:
            rvir,vvir = normtohost
            r = r/rvir
            v = v/vvir
            xerr = xerr/rvir
            yerr = yerr/vvir
        ax.scatter(r,v,s=self.loglum,marker='s',color='k',**kwargs)
        ax.errorbar(r,v,xerr=xerr,yerr=yerr,ecolor='k',fmt='none')
    def loadwolf10(self,minlum=None):
        tab = asciitable.read('/bigbang/data/AnnaGroup/literaturedata/Wolf2010.csv')
        if minlum != None:
            tab = tab[tab['lum']>=minlum]

        G = 4.34e-6 #(km/s)^2 kpc/Msun
        rhalf = tab['rhalf']/1000.
        rhalf_errlo = -tab['rhalf_err2']/1000.
        rhalf_errhi =  tab['rhalf_err1']/1000.
        rhalf_lo = rhalf-rhalf_errlo
        rhalf_hi = rhalf+rhalf_errhi
        
        Mhalf = tab['Mhalf']
        vhalf = np.sqrt(G*Mhalf/rhalf)
        
        Mhalf_hi = tab['Mhalf']+tab['Mhalf_err1']
        Mhalf_lo = tab['Mhalf']+tab['Mhalf_err2']
        Mhalf_errlo = Mhalf-Mhalf_lo
        Mhalf_errhi = Mhalf_hi-Mhalf
        
        vhalf_lo = np.sqrt(G*Mhalf_lo/rhalf)
        vhalf_hi = np.sqrt(G*Mhalf_hi/rhalf)
        vhalf_lolo = np.sqrt(G*Mhalf_lo/rhalf_hi)
        vhalf_hihi = np.sqrt(G*Mhalf_hi/rhalf_lo)
        vhalf_errlo = vhalf-vhalf_lo
        vhalf_errhi = vhalf_hi-vhalf
        vhalf_errlolo = vhalf-vhalf_lolo
        vhalf_errhihi = vhalf_hihi-vhalf
        
        loglum = np.log10(tab['lum'])
        
        xerr = [rhalf_errlo,rhalf_errhi]
        yerr = [vhalf_errlo,vhalf_errhi]
        yerr2 = [vhalf_errlolo,vhalf_errhihi]

        self.wolf10tab = tab
        self.rhalf = rhalf
        self.rerr  = xerr
        self.vhalf = vhalf
        self.verr  = yerr
        self.verr2 = yerr2
        self.loglum= loglum
        
class TBTFSoftPlugin(SubVelocityProfileSoftPlugin,TBTFPlugin):
    def __init__(self,minlum=2.e5):
        TBTFPlugin.__init__(self,minlum=minlum)
        super(TBTFSoftPlugin,self).__init__()
        self.autofigname = 'tbtfsoft'
        self.xmin = .1; self.xmax = 2
        self.ymin = 8;  self.ymax = 50
        self.vcut = 30.
        self.vLMC = 60.
        self.vSMC = 40.

    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,alpha=.3,**kwargs):
        assert lx != None
        color = self.colordict[lx]

        #TODO factor out grabbing vpeak?
        extantplug = MTanalysis.TagExtantPlugin()
        extantdata = extantplug.read(hpath,autocalc=False)
        if extantdata == None:
            print "ERROR: %s does not have Greg's Extant data (not plotting)" % (haloutils.get_foldername(hpath))
            return
        extantids,extantdata = extantdata
        vpeakarr = extantdata['vpeak']
        keeprsids = (extantdata['rsid'][vpeakarr >= self.vcut]).astype(int)
        vpeak = pd.Series(vpeakarr,index=extantdata['rsid'])

        rsid,rarr,rvir,vcircarr,vcircsoftarr = data
        vcircarr = vcircsoftarr
        #rsid,rarr,rvir,vcircarr = data
        iicut = np.in1d(rsid,keeprsids,assume_unique=True)
        rsid = rsid[iicut]
        rvir = rvir[iicut]
        vpeak = np.array(vpeak.ix[rsid])
        nbound = len(rsid)
        nall = len(keeprsids)
        vcircarr = vcircarr[iicut,:]
        
        iicut2 = vpeak <= min(self.vLMC,self.vSMC)
        nNotMC = np.sum(iicut2)

        mc_rarr = rarr[~iicut2,:]
        mc_vcircarr = vcircarr[~iicut2,:]
        mc_rsid = rsid[~iicut2]
        mc_rvir = rvir[~iicut2]
        mc_vpeak= vpeak[~iicut2]

        rarr = rarr[iicut2,:]
        vcircarr = vcircarr[iicut2,:]
        rsid = rsid[iicut2]
        rvir = rvir[iicut2]
        vpeak = vpeak[iicut2]

        rarr = rarr*1000 #kpc
        mc_rarr = mc_rarr*1000
        eps = 1000*haloutils.load_soft(hpath)
        if normtohost:
            mvir,rvir,vvir=haloutils.load_haloprops(hpath)
            rarr = rarr/rvir
            vcircarr = vcircarr/vvir
            eps = eps/rvir
            normtohost=(rvir,vvir)
        self.plotwolf10(ax,normtohost=normtohost)
        for i in xrange(len(rsid)):
            ii = rarr[i,:] >= eps
            if np.sum(ii) == 0: continue
            ax.plot(rarr[i,ii], vcircarr[i,ii], color=color, alpha=alpha, **kwargs)
        for i in xrange(len(mc_rsid)):
            ii = mc_rarr[i,:] >= eps
            if np.sum(ii) == 0: continue
            ax.plot(mc_rarr[i,ii], mc_vcircarr[i,ii], color=color, linestyle=':', alpha=alpha, **kwargs)
        xmin,xmax,ymin,ymax,xlog,ylog,xlabel,ylabel = self.get_plot_params(normtohost)
        logxoff = np.log10(xmax/xmin)*.05
        xlabel  = xmax * 10**(-logxoff)
        logyoff = np.log10(ymax/ymin)*.1
        ylabel  = ymax * 10**(-logyoff)
        ax.text(xlabel,ylabel,r"$%i/%i/%i$" % (nNotMC,nbound,nall),ha='right')

class NvinfallPlugin(MTanalysis.TagExtantPlugin,NvmaxPlugin):
    def __init__(self):
        # First init Nvmax, then init TagExtantPlugin
        NvmaxPlugin.__init__(self)
        super(NvinfallPlugin,self).__init__()

        self.autofigname = 'Nvinfall'
        self.vcut = 30.
        self.vLMC = 60.
        self.vSMC = 40.
    
        #self.logvmin = -1.
        #self.logvmax = 3.
        #self.dlogv = 0.05
        #self.vmaxbins = 10.**np.arange(self.logvmin,self.logvmax+self.dlogv,self.dlogv)

        self.xmin=20; self.xmax=80
        self.ymin=0; self.ymax=20
        self.xlog= False; self.ylog = False
        self.xlabel=r'$v_{\rm peak}\ (km/s)$' ; self.ylabel=r'$N(>v_{\rm peak})$'
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        ids,tab = data
        vpeak = tab['vpeak']
        minv = np.min(vpeak)
        Nvinfall = self.calcNvmax(vpeak) #uses self.vmaxbins
        vplot = self.vmaxbins[1:]
        ii = vplot >= minv

        if lx != None:
            ax.plot(vplot[ii],Nvinfall[ii],color=self.colordict[lx],**kwargs)
        else:
            ax.plot(vplot[ii],Nvinfall[ii],color=self.colordict[lx],**kwargs)
