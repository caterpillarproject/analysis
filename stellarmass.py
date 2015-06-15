import numpy as np
import pylab as plt
import os,sys,subprocess,time

import haloutils
from SAMs import SimpleSAMBasePlugin
#from fast_SAMs import FastSAMBasePlugin as SimpleSAMBasePlugin
import abundmatch

class AMStellarMassPlugin(SimpleSAMBasePlugin):
    def __init__(self,verbose=False):
        super(AMStellarMassPlugin,self).__init__(verbose=verbose)

        self.xmin = 1e3; self.xmax = 1e11
        self.ymin = 1; self.ymax = 300
        self.xlabel = r'$M_{\rm star}\ (M_\odot)$'
        self.ylabel = r'$N(>M_{\rm star})$'
        self.xlog=True; self.ylog=True
        self.autofigname='AMStellarMass'

        M13AM = abundmatch.Moster13AbundMatch()
        B13AM = abundmatch.Behroozi13AbundMatch()
        G14AM = abundmatch.GK14AbundMatch()
        self.AMlist = [M13AM,B13AM,G14AM]
        self.namelist = ['M13','B13','GK14']
        self.linestyles = ['-','--',':']
        #self.yticks = np.concatenate((np.arange(9)+1,10*(np.arange(9)+1),[100,200,300]))

    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        subs = data
        if normtohost: raise NotImplementedError
        Minfall = subs['max_mass']
        Mstarlist = [AM.get_Mstar(Minfall) for AM in self.AMlist]
        thiscol = None
        for Mstar,linestyle in zip(Mstarlist,self.linestyles):
            Mstar = np.sort(np.array(Mstar[~np.isnan(Mstar)]))[::-1]
            count = np.arange(len(Mstar))+1
            if lx != None:
                ax.plot(Mstar,count,linestyle=linestyle,color=self.colordict[lx],**kwargs)
            else:
                if thiscol==None:
                    l, = ax.plot(Mstar,count,linestyle=linestyle,**kwargs)
                    thiscol = plt.getp(l,'color')
                    kwargs.pop('color',None)
                else:
                    ax.plot(Mstar,count,linestyle=linestyle,color=thiscol,**kwargs)

class AMStellarMassz8CutPlugin(AMStellarMassPlugin):
    def __init__(self,**kwargs):
        super(AMStellarMassz8CutPlugin,self).__init__(**kwargs)
        self.Mcut = 1.e8 #Msun/h

    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        subs = data
        if normtohost: raise NotImplementedError
        Minfall = subs['infall_mvir']
        ii = subs['z8_mvir']>self.Mcut
        Mstarlist = [AM.get_Mstar(Minfall[ii]) for AM in self.AMlist]
        thiscol = None
        for Mstar,linestyle in zip(Mstarlist,self.linestyles):
            Mstar = np.sort(np.array(Mstar[~np.isnan(Mstar)]))[::-1]
            count = np.arange(len(Mstar))+1
            if lx != None:
                ax.plot(Mstar,count,linestyle=linestyle,color=self.colordict[lx],**kwargs)
            else:
                if thiscol==None:
                    l, = ax.plot(Mstar,count,linestyle=linestyle,**kwargs)
                    thiscol = plt.getp(l,'color')
                    kwargs.pop('color',None)
                else:
                    ax.plot(Mstar,count,linestyle=linestyle,color=thiscol,**kwargs)

if __name__=="__main__":
    #bins = np.arange(-0.5,12.5+.5,.5)

    hpath = haloutils.get_hpath_lx(5320,14)
    plug = SimpleSAMBasePlugin()
    M13AM = abundmatch.Moster13AbundMatch()
    B13AM = abundmatch.Behroozi13AbundMatch()
    G14AM = abundmatch.GK14AbundMatch()
    AMlist = [M13AM,B13AM,G14AM]
    namelist = ['M13','B13','GK14']
    
    subs = plug.read(hpath)
    Minfall = subs['infall_mvir']
    Mstarlist = [AM.get_Mstar(Minfall) for AM in AMlist]
    fig,ax = plt.subplots()
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'$M_{\rm star}\ (M_\odot)$'); ax.set_ylabel(r'$N(>M_{\rm star})$')
    ax.set_ylim((1,300)); ax.set_xlim((1e3,1e11))
    yticks = np.concatenate((np.arange(9)+1,10*(np.arange(9)+1),[100,200,300]))
    ax.set_yticks(yticks)
    for Mstar,name in zip(Mstarlist,namelist):
        Mstar = np.sort(np.array(Mstar[~np.isnan(Mstar)]))[::-1]
        count = np.arange(len(Mstar))+1
        ax.plot(Mstar,count,label=name)
        #h,x = np.histogram(np.log10(Mstar),bins=bins)
        #x = (x[1:]+x[:-1])/2.
        #ax.plot(x,h,drawstyle='steps-mid',label=name)
    ax.legend()
    plt.show()
