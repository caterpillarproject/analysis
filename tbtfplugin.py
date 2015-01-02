import numpy as np
import pylab as plt
import os,subprocess,sys,time
import asciitable
import pynbody as pnb
import readsnapshots.readsnapHDF5_greg as rsg
import haloutils

from caterpillaranalysis import *

class TBTFPlugin(SubVelocityProfilePlugin):
    def __init__(self):
        super(TBTFPlugin,self).__init__()
        self.xmin = .1; self.xmax = 2
        self.ymin = 8;  self.ymax = 50
        self.autofigname = 'tbtf'

        self.loadwolf10(minlum=2*10**5)
        self.vcut = 20
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,alpha=.2,**kwargs):
        assert lx != None
        color = self.colordict[lx]

        #TODO load Vpeak > 30km/s
        rsid,rarr,rvir,vcircarr = data
        vmaxarr = np.max(vcircarr,axis=1)
        iicut = vmaxarr >= self.vcut
        rsid = rsid[iicut]
        rvir = rvir[iicut]
        vcircarr = vcircarr[iicut,:]

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
        
