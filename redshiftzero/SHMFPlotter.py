import numpy as np
import pylab as plt
import os,asciitable
import haloutils,sheetplot

class SHMFReader(sheetplot.ReaderBase):
    def __init__(self,ictype="BB",nv=4):
        self.ictype = ictype;
        self.nv = nv
        self.filename = 'SHMF.dat'
    def __call__(self,hid):
        ictype = self.ictype; nv = self.nv
        lxlist = self.get_lxlist(hid)
        
        xlist = []; ylist = []
        sxlist = []; sylist = []
        for lx in lxlist:
            thisfilename = self.get_filename(hid,lx)
            data = asciitable.read(thisfilename,delimiter=' ')
            xlist.append(data['col1'])
            ylist.append(data['col2'])
            sxlist.append(data['col3'])
            sylist.append(data['col4'])
        return hid,lxlist,xlist,ylist,sxlist,sylist

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
