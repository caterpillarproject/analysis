import pylab as plt
import matplotlib
import numpy as np
import os,asciitable
import haloutils
from brendanlib.grifflib import makecolormap
import seaborn as sns
from profiles.profilefit import NFWprofile,fitNFW
sns.set_context("poster")
sns.set_style('ticks',{'xtick.direction':'in','ytick.direction':'in'})
snsdict = dict(zip(plt.rcParamsDefault['axes.color_cycle'],plt.rcParams['axes.color_cycle']))

class PlotterBase(object):
    colordict = {11:snsdict['b'],12:snsdict['r'],13:snsdict['g'],14:snsdict['m']}

class NvmaxPlotter(PlotterBase):
    def __init__(self):
        self.fprefix = 'Nvmax'
        self.fpostfix = 'rs'
    def __call__(self,ax,data,fignum=None,**kwargs):
        try:
            hid,lxlist,vlist,Nlist,minvlist,sNlist,sminvlist,Nplist,sNplist = data
        except TypeError:
            hid = data; lxlist = []
        vmin = .3; vmax = 100
        Nmin = 1; Nmax = 10**4.5
        ax.set_xscale('log'); ax.set_yscale('log')
        for i,lx in enumerate(lxlist):
            color = self.colordict[lx]
            v = vlist[i]; N = Nlist[i]
            if v==None: continue
            minv = minvlist[i]
            ii = np.where(v >= minv)
            ax.plot(v[ii],N[ii],color=color,**kwargs)
        ax.set_xlabel(r'$V_{\rm max}$ [km/s]')
        ax.set_ylabel(r'N($>V_{\rm max}$)')
        ax.set_xlim([vmin,vmax])
        ax.set_ylim([Nmin,Nmax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(vmin*10**0.1,Nmax*10**-.5,plotlabel,color='black',fontsize='medium')

class NvmaxpPlotter(PlotterBase):
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
            if v==None: continue
            minv = minvlist[i]
            ii = np.where(v >= minv)
            ax.plot(v[ii],N[ii],color=color,**kwargs)
        ax.set_xlabel(r'$V_{\rm max}$ [km/s]')
        ax.set_ylabel(r'N($>V_{\rm max}$)')
        ax.set_xlim([vmin,vmax])
        ax.set_ylim([Nmin,Nmax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(vmin*10**0.1,Nmax*10**-.5,plotlabel,color='black',fontsize='medium')

class sNvmaxPlotter(PlotterBase):
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
            if v==None: continue
            minv = sminvlist[i]
            ii = np.where(v >= minv)
            ax.plot(v[ii],N[ii],color=color,**kwargs)
        ax.set_xlabel(r'$V_{\rm max}$ [km/s]')
        ax.set_ylabel(r'N($>V_{\rm max}$)')
        ax.set_xlim([vmin,vmax])
        ax.set_ylim([Nmin,Nmax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(vmin*10**0.1,Nmax*10**-.5,plotlabel,color='black',fontsize='medium')

class sNvmaxpPlotter(PlotterBase):
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
            if v==None: continue
            minv = sminvlist[i]
            ii = np.where(v >= minv)
            ax.plot(v[ii],N[ii],color=color,**kwargs)
        ax.set_xlabel(r'$V_{\rm max}$ [km/s]')
        ax.set_ylabel(r'N($>V_{\rm max}$)')
        ax.set_xlim([vmin,vmax])
        ax.set_ylim([Nmin,Nmax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(vmin*10**0.1,Nmax*10**-.5,plotlabel,color='black',fontsize='medium')

class SHMFPlotter(PlotterBase):
    def __init__(self):
        self.fprefix = 'SHMF'
        self.fpostfix = 'rs'
    def __call__(self,ax,data,fignum=None,**kwargs):
        try:
            hid,lxlist,xlist,ylist,sxlist,sylist = data
        except TypeError:
            hid = data; lxlist = []
        xmin = 10**4.5; xmax = 10**10.6
        ymin = 10**-10; ymax = 10**-1.0
        ax.set_xscale('log'); ax.set_yscale('log')
        for i,lx in enumerate(lxlist):
            color = self.colordict[lx]
            x = xlist[i]; y = ylist[i]
            if x==None: continue
            #ii = y>ymin
            #ax.plot(x[ii],y[ii],color=color,**kwargs)
            ax.plot(x,y,color=color,**kwargs)
        ax.set_xlabel(r'$M_{\rm sub} [h^{-1} M_\odot]$')
        ax.set_ylabel(r'$dn/dM_{\rm sub} [h/M_\odot]$')
        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(xmin*10**0.2,ymax*10**-.5,plotlabel,color='black',fontsize='medium')

class sSHMFPlotter(PlotterBase):
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
            if x==None: continue
            #ii = y>ymin
            #ax.plot(x[ii],y[ii],color=color,**kwargs)
            ax.plot(x,y,color=color,**kwargs)
        ax.set_xlabel(r'$M_{\rm sub} [h^{-1} M_\odot]$')
        ax.set_ylabel(r'$dn/dM_{\rm sub} [h/M_\odot]$')
        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(xmin*10**0.2,ymax*10**-.5,plotlabel,color='black',fontsize='medium')

class ProfilePlotter(PlotterBase):
    def __init__(self):
        self.fprefix = 'rhor2'
        self.fpostfix= 'rs'
    def __call__(self,ax,data,fignum=None,**kwargs):
        try:
            hid,lxlist,rlist,rholist,p03rlist,rvirlist,r200clist = data
        except TypeError:
            hid = data; lxlist = []
        ymin = 10**-1.5; ymax = 10**2.5
        ax.set_xscale('log'); ax.set_yscale('log')
        for i,lx in enumerate(lxlist):
            color = self.colordict[lx]
            r = rlist[i]; rho = rholist[i]
            p03r = p03rlist[i]; rvir = rvirlist[i]; r200c = r200clist[i]
            if r == None: continue
            ax.plot(r,(r/1000.)**2 * rho, color=color, **kwargs)
            ax.plot([p03r,p03r],[ymin,ymax],color=color,ls='--',**kwargs)
            ax.plot([rvir,rvir],[ymin,ymax],'k-.')
            ax.plot([r200c,r200c],[ymin,ymax],'k:')

            good_r = (r > p03r) & (r < rvir)
            rs,rhos = fitNFW(r[good_r],rho[good_r],x0=[20,6],bounds=[(1,300),(4,8)])
            ax.plot(r,(r/1000.)**2 * NFWprofile(r,rs,rhos),ls=':',lw='0.5',**kwargs)

            yexponent = 1.3 + 0.2*(14-lx)
            ax.text(10**-1.8,10**yexponent,r"LX%i $r_s=$%3.2f kpc" % (lx,rs),
                    color=color,fontsize='x-small')

        ax.set_xlabel(r'r [$h^{-1}$ kpc]')
        ax.set_ylabel(r'$r^2 \rho(r)$ [$10^{10} M_\odot$ Mpc$^{-1}$]')
        ax.set_xlim([10**-2,10**3])
        ax.set_ylim([ymin,ymax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(10**-1.8,10**2.2,plotlabel,color='black',fontsize='small')
        
class ProjPlotter(PlotterBase):
    def __init__(self,lx,snap=255,vmin=10**6,vmax=10**9):
        self.lx=lx
        self.snap=snap
        self.fprefix='proj'+str(lx)
        if snap==255:
            self.fpostfix='xy'
        else:
            self.fpostfix='xy'+str(snap).zfill(3)
        self.cmap = makecolormap()
        self.linthresh=None
        self.vmin=vmin; self.vmax=vmax
    def __call__(self,ax,data,fignum=None,**kwargs):
        """Replicate pynbody.plot.sph.image"""
        hid,lx,im,width = data
        assert lx==self.lx
        try:
            im[np.where(im==0)] = abs(im[np.where(abs(im != 0))]).min()
        except ValueError:
            raise ValueError, "Failed to make a sensible logarithmic image. This probably means there are no particles in the view."
        if (im < 0).any():
            if self.linthresh is None:
                linthresh = np.nanmax(abs(im)) / 1000.
            else: linthresh = self.linthresh
            norm = matplotlib.colors.SymLogNorm(
                linthresh, vmin=self.vmin, vmax=self.vmax)
        else:
            norm = matplotlib.colors.LogNorm(vmin=self.vmin, vmax=self.vmax)
        ax.imshow(im[::-1,:].view(np.ndarray), extent=(-width/2.,width/2.,-width/2.,width/2.),
                  vmin=self.vmin,vmax=self.vmax,cmap=self.cmap, norm=norm)
        #aspect='auto' #fills up box
        ax.text(-width/2.+width*.05,width/2.-width*.1,haloutils.hidstr(hid),color='white')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim([-width/2.,width/2.]); ax.set_ylim([-width/2.,width/2.])

class MassAccrPlotter(PlotterBase):
    def __init__(self):
        self.fprefix='MT'
        self.fpostfix=''
    def __call__(self,ax,data,fignum=None,**kwargs):
        try:
            hid,lxlist,tablist = data
        except TypeError as e:
            hid = data; lxlist = []
        xmin = 0; xmax = 1
        ymin = 10**6; ymax = 10.**13
        ax.set_yscale('log')
        for i,lx in enumerate(lxlist):
            color = self.colordict[lx]
            tab = tablist[i]
            if tab==None: continue
            x = tab['scale']
            y = tab['mvir']
            ax.plot(x,y,color=color,**kwargs)
        ax.set_xlabel('scale')
        ax.set_ylabel(r'$M$ [$M_\odot$]')
        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(0.1,ymax*10**-.5,plotlabel,color='black',fontsize='medium')
