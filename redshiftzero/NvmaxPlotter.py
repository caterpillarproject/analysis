import numpy as np
import pylab as plt
import os,asciitable
import haloutils,sheetplot

class NvmaxReader(sheetplot.ReaderBase):
    def __init__(self,ictype="BB",nv=4):
        self.ictype = ictype;
        self.nv = nv
        self.filename = 'Nvmax.dat'
    def __call__(self,hid):
        ictype = self.ictype; nv = self.nv
        lxlist = self.get_lxlist(hid)
        
        vlist = []; Nlist = []; sNlist = []
        Nplist = []; sNplist = []
        minvlist = []; sminvlist = []
        for lx in lxlist:
            thisfilename = self.get_filename(hid,lx)
            data = asciitable.read(thisfilename,delimiter=' ',data_start=1)
            
            vlist.append(data['col1'])
            Nlist.append(data['col2'])
            sNlist.append(data['col3'])
            Nplist.append(data['col4'])
            sNplist.append(data['col5'])
            f = open(thisfilename,'r')
            split = f.readline().split(" ")
            minvlist.append(float(split[0])); sminvlist.append(float(split[1]))
            f.close()
        return hid,lxlist,vlist,Nlist,minvlist,sNlist,sminvlist,Nplist,sNplist

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

