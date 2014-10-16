import pylab as plt
import numpy as np
import os,asciitable
import haloutils

class PlotterBase(object):
    colordict = {11:'blue',12:'red',13:'green',14:'purple'}

class NvmaxPlotter(sheetplot.PlotterBase):
    def __init__(self):
        self.fprefix = 'Nvmax'
        self.fpostfix = 'rs'
    def __call__(self,ax,data,fignum=None,**kwargs):
        hid,lxlist,vlist,Nlist,minvlist,sNlist,sminvlist,Nplist,sNplist = data
        vmin = .3; vmax = 100
        Nmin = 1; Nmax = 10**4.5
        ax.set_xscale('log'); ax.set_yscale('log')
        for i,lx in enumerate(lxlist):
            color = self.colordict[lx]
            v = vlist[i]; N = Nlist[i]
            minv = minvlist[i]
            ii = np.where(v >= minv)
            ax.plot(v[ii],N[ii],color=color,**kwargs)
        ax.set_xlabel(r'$V_{\rm max}$ [km/s]')
        ax.set_ylabel(r'N($>V_{\rm max}$)')
        ax.set_xlim([vmin,vmax])
        ax.set_ylim([Nmin,Nmax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(vmin*10**0.1,Nmax*10**-.3,plotlabel,color='black',fontsize='small')

class NvmaxpPlotter(sheetplot.PlotterBase):
    def __init__(self):
        self.fprefix = 'Nvmaxp'
        self.fpostfix = 'rs'
    def __call__(self,ax,data,fignum=None,**kwargs):
        hid,lxlist,vlist,Nlist,minvlist,sNlist,sminvlist,Nplist,sNplist = data
        vmin = .3; vmax = 100
        Nmin = 1; Nmax = 10**4.5
        ax.set_xscale('log'); ax.set_yscale('log')
        for i,lx in enumerate(lxlist):
            color = self.colordict[lx]
            v = vlist[i]; N = Nplist[i]
            minv = minvlist[i]
            ii = np.where(v >= minv)
            ax.plot(v[ii],N[ii],color=color,**kwargs)
        ax.set_xlabel(r'$V_{\rm max}$ [km/s]')
        ax.set_ylabel(r'N($>V_{\rm max}$)')
        ax.set_xlim([vmin,vmax])
        ax.set_ylim([Nmin,Nmax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(vmin*10**0.1,Nmax*10**-.3,plotlabel,color='black',fontsize='small')

class sNvmaxPlotter(sheetplot.PlotterBase):
    def __init__(self):
        self.fprefix = 'Nvmax'
        self.fpostfix = 'subf'
    def __call__(self,ax,data,fignum=None,**kwargs):
        hid,lxlist,vlist,Nlist,minvlist,sNlist,sminvlist,Nplist,sNplist = data
        vmin = .3; vmax = 100
        Nmin = 1; Nmax = 10**4.5
        ax.set_xscale('log'); ax.set_yscale('log')
        for i,lx in enumerate(lxlist):
            color = self.colordict[lx]
            v = vlist[i]; N = sNlist[i]
            minv = sminvlist[i]
            ii = np.where(v >= minv)
            ax.plot(v[ii],N[ii],color=color,**kwargs)
        ax.set_xlabel(r'$V_{\rm max}$ [km/s]')
        ax.set_ylabel(r'N($>V_{\rm max}$)')
        ax.set_xlim([vmin,vmax])
        ax.set_ylim([Nmin,Nmax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(vmin*10**0.1,Nmax*10**-.3,plotlabel,color='black',fontsize='small')

class sNvmaxpPlotter(sheetplot.PlotterBase):
    def __init__(self):
        self.fprefix = 'Nvmaxp'
        self.fpostfix = 'subf'
    def __call__(self,ax,data,fignum=None,**kwargs):
        hid,lxlist,vlist,Nlist,minvlist,sNlist,sminvlist,Nplist,sNplist = data
        vmin = .3; vmax = 100
        Nmin = 1; Nmax = 10**4.5
        ax.set_xscale('log'); ax.set_yscale('log')
        for i,lx in enumerate(lxlist):
            color = self.colordict[lx]
            v = vlist[i]; N = sNplist[i]
            minv = sminvlist[i]
            ii = np.where(v >= minv)
            ax.plot(v[ii],N[ii],color=color,**kwargs)
        ax.set_xlabel(r'$V_{\rm max}$ [km/s]')
        ax.set_ylabel(r'N($>V_{\rm max}$)')
        ax.set_xlim([vmin,vmax])
        ax.set_ylim([Nmin,Nmax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(vmin*10**0.1,Nmax*10**-.3,plotlabel,color='black',fontsize='small')

class SHMFPlotter(sheetplot.PlotterBase):
    def __init__(self):
        self.fprefix = 'SHMF'
        self.fpostfix = 'rs'
    def __call__(self,ax,data,fignum=None,**kwargs):
        hid,lxlist,xlist,ylist,sxlist,sylist = data
        xmin = 10**5; xmax = 10**10.6
        ymin = 10**-10; ymax = 10**-2
        ax.set_xscale('log'); ax.set_yscale('log')
        for i,lx in enumerate(lxlist):
            color = self.colordict[lx]
            x = xlist[i]; y = ylist[i]
            ax.plot(x,y,color=color,**kwargs)
        ax.set_xlabel(r'$M_{\rm sub} [h^{-1} M_\odot]$')
        ax.set_ylabel(r'$dn/dM_{\rm sub} [h/M_\odot]$')
        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(xmin*10**0.1,ymax*10**-.3,plotlabel,color='black',fontsize='small')

class sSHMFPlotter(sheetplot.PlotterBase):
    def __init__(self):
        self.fprefix = 'SHMF'
        self.fpostfix = 'subf'
    def __call__(self,ax,data,fignum=None,**kwargs):
        hid,lxlist,xlist,ylist,sxlist,sylist = data
        xmin = 10**5; xmax = 10**10.6
        ymin = 10**-10; ymax = 10**-2
        ax.set_xscale('log'); ax.set_yscale('log')
        for i,lx in enumerate(lxlist):
            color = self.colordict[lx]
            x = sxlist[i]; y = sylist[i]
            ax.plot(x,y,color=color,**kwargs)
        ax.set_xlabel(r'$M_{\rm sub} [h^{-1} M_\odot]$')
        ax.set_ylabel(r'$dn/dM_{\rm sub} [h/M_\odot]$')
        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(xmin*10**0.1,ymax*10**-.3,plotlabel,color='black',fontsize='small')
