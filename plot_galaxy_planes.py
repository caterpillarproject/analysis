import matplotlib; matplotlib.use('Agg')
import numpy as np
import pylab as plt
import pandas as pd
import os,sys,subprocess,time,functools

import haloutils
import caterpillarplot
import abundmatch,stellarmass
from SAMs_old import SimpleSAMBasePlugin
from galaxy_planes import SatellitePlanesPlugin

from seaborn.apionly import cubehelix_palette

global_plot_sams = ['Ni11','RNi11','L0i1','L1i1']

class SatellitePlanesAnglePlotter(SatellitePlanesPlugin):
    def __init__(self,plot_sams=global_plot_sams):
        super(SatellitePlanesAnglePlotter,self).__init__()
        self.xmin = 0; self.xmax = 90
        self.ymin = 0; self.ymax = 20
        self.xlabel = r'$\theta\ \rm{(deg)}$'
        self.ylabel = r'$N(<\theta)$'
        self.xlog=False; self.ylog=False
        self.autofigname='satplanesangle'
        self.plot_sams = plot_sams
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        assert lx==14 or lx==None
        df,satix = data
        xaxis = np.concatenate([[0],self.thetaarr])
        isotropic = np.sin(np.array(xaxis)*np.pi/180.)
        for sam in self.plot_sams:
            row = df.ix[sam]
            N_lt_theta = np.concatenate([[0],row[self.thetaarrnames]])
            l, = ax.plot(xaxis,N_lt_theta,label=sam,**kwargs)
            color = plt.getp(l,'color')
            ax.plot(xaxis,N_lt_theta[-1]*isotropic,':',color=color,**kwargs)

class SatellitePlanesRadialPlotter(SatellitePlanesPlugin):
    def __init__(self,plot_sams=global_plot_sams):
        super(SatellitePlanesRadialPlotter,self).__init__()
        self.xmin = 0; self.xmax = 300
        self.ymin = 0; self.ymax = 1
        self.xlabel = r'$r\ \rm{(kpc)}$'
        self.ylabel = r'$\rm{fraction}(<r)$'
        self.xlog=False; self.ylog=False
        self.autofigname='satplanesradial'
        self.plot_sams = plot_sams
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        assert lx==14 or lx==None
        df,satix = data
        subs = self.samplug.read(hpath)
        for sam in self.plot_sams:
            ix = satix[sam]
            sats = subs.ix[ix]
            dr = np.sort(np.array(sats['dr']/self.h0))
            Nsats = float(len(dr))
            ax.plot(dr,(np.arange(Nsats)+1)/Nsats,label=sam,**kwargs)

class SatellitePlanesBACAPlotter(SatellitePlanesPlugin):
    def __init__(self,plot_sams=global_plot_sams):
        super(SatellitePlanesBACAPlotter,self).__init__()
        self.xmin = 0; self.xmax = 1
        self.ymin = 0; self.ymax = 1
        self.xlabel = r'b/a'
        self.ylabel = r'c/a'
        self.xlog=False; self.ylog=False
        self.autofigname='satplanesbaca'
        self.plot_sams = plot_sams
        self.cm_df = pd.read_csv('concmass14.csv',index_col=0)
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,labelconc=True,**kwargs):
        assert lx==14 or lx==None
        df,satix = data
        for sam in self.plot_sams:
            row = df.ix[sam]
            ba = row['ba']; ca = row['ca']
            ax.plot(ba,ca,label=sam,**kwargs)
        ax.plot([0,1],[0,1],'k:')
        if labelconc:
            try:
                hid = haloutils.get_parent_hid(hpath)
                conc = np.log10(self.cm_df['conc'].ix[hid])
                ax.text(.1,.8,r'$\log\ c={0:4.2f}$'.format(conc),fontsize='small')
            except KeyError as e:
                print 'Missing HID:',e,'skipping conc...'

class SatellitePlanesPlanesizePlotter(SatellitePlanesPlugin):
    def __init__(self,plot_sams=global_plot_sams):
        super(SatellitePlanesPlanesizePlotter,self).__init__()
        self.xmin = 0; self.xmax = 200
        self.ymin = 0; self.ymax = 100
        self.xlabel = r'$r_{\rm par}\ \rm{(kpc)}$'
        self.ylabel = r'$r_{\rm perp}\ \rm{(kpc)}$'
        self.xlog=False; self.ylog=False
        self.autofigname='satplanesplanesize'
        self.plot_sams = plot_sams
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        assert lx==14 or lx==None
        df,satix = data
        for sam in self.plot_sams:
            row = df.ix[sam]
            rperp = row['rperp']; rpar = row['rpar']
            ax.plot(rpar,rperp,label=sam,marker='o',**kwargs)
        ax.plot([0,250],[0,250],'k:')

class SatellitePlanesFaceOnPlotter(SatellitePlanesPlugin):
    def __init__(self):
        super(SatellitePlanesFaceOnPlotter,self).__init__()
        self.xmin = -250; self.xmax = 250
        self.ymin = -250; self.ymax = 250
        self.xlabel = r'$\rm{Major\ Axis}$'
        self.ylabel = r'$\rm{Middle\ Axis}$'
        self.xlog=False; self.ylog=False
        self.autofigname='satplanesfaceon'
        self.cmap = cubehelix_palette(as_cmap=True,start=.5,rot=-1.5,hue=1.0,gamma=1.0)
    def process_subs(self,subs,hpos,hvel,U=None):
        """ Rotate and normalize positions"""
        #Center
        spos = np.array(subs[['posX','posY','posZ']])-hpos
        svel = np.array(subs[['pecVX','pecVY','pecVZ']])-hvel
        #Rotate using U
        if U==None:
            ba,ca,rperp,rpar,ntheta,U = self.calculate_plane_params(spos,svel,return_vecs=True)
        spos = U.T.dot(spos.T).T
        svel = U.T.dot(svel.T).T
        #spos in kpc from host center
        spos *= 1000./self.h0
        #svel normalized to point how far it goes in 1Myr
        #1 km/s -> 0.1kpc/Myr
        svel *= 0.1
        return spos,svel,U
    def vmax_to_size(self,svmax):
        min_vmax_size=0.
        max_vmax_size=100.
        min_vmax=0.
        max_vmax=100.
        normed_vmax = np.array((svmax-min_vmax)/(max_vmax-min_vmax))
        return normed_vmax * (max_vmax_size - min_vmax_size) + min_vmax_size
    def logMstar_to_color(self,logMstar):
        min_color=0.
        max_color=255.
        min_logMstar=5.
        max_logMstar=10.
        normed_logMstar = np.array((logMstar-min_logMstar)/(max_logMstar-min_logMstar))
        return normed_logMstar * (max_color - min_color) + min_color

    def _plot_axes(self,xax,yax,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        assert lx==14 or lx==None
        df,satix = data
        subs = self.samplug.read(hpath)
        mb = self.mbplug.read(hpath)
        host = mb[-1]
        hpos = np.array([host['x'],host['y'],host['z']])
        hvel = np.array([host['vx'],host['vy'],host['vz']])

        #plot first 11 by Vinfall
        Nsubs = subs.ix[satix['Ni11']]
        Nspos,Nsvel,U = self.process_subs(Nsubs,hpos,hvel,U=None)
        ax.scatter(Nspos[:,xax],Nspos[:,yax],c='k',s=4)
        for pos,vel in zip(Nspos,Nsvel):
            ax.arrow(pos[xax],pos[yax],vel[xax],vel[yax],head_width=5, head_length=5, fc='k', ec='k')

        #plot Mstar > 10^5 Msun with Moster et al 2013 AM
        Lsubs = subs.ix[satix['L1i1']]
        Lspos,Lsvel,U = self.process_subs(Lsubs,hpos,hvel,U=U)

        logMstar = np.log10(self.AMlist[1].get_Mstar(np.array(Lsubs['infall_mvir'])))
        colors= self.logMstar_to_color(logMstar)
        #Lsvmax = np.array(Lsubs['vmax'])
        #sizes = self.vmax_to_size(Lsvmax)
        ax.scatter(Lspos[:,xax],Lspos[:,yax],c=colors,cmap=self.cmap)
        for pos,vel in zip(Lspos,Lsvel):
            ax.arrow(pos[xax],pos[yax],vel[xax],vel[yax],head_width=5, head_length=5, fc='k', ec='k')
        
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        self._plot_axes(0,1,hpath,data,ax,lx,labelon,normtohost,**kwargs)

class SatellitePlanesEdgeOnPlotter(SatellitePlanesFaceOnPlotter):
    def __init__(self):
        super(SatellitePlanesEdgeOnPlotter,self).__init__()
        self.xmin = -250; self.xmax = 250
        self.ymin = -250; self.ymax = 250
        self.xlabel = r'$\rm{Major\ Axis}$'
        self.ylabel = r'$\rm{Minor\ Axis}$'
        self.xlog=False; self.ylog=False
        self.autofigname='satplanesedgeon'
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        self._plot_axes(0,2,hpath,data,ax,lx,labelon,normtohost,**kwargs)

def plot_ca_2x2():
    pass

if __name__=="__main__":
    plug = SatellitePlanesPlugin()
    angle = SatellitePlanesAnglePlotter()
    radial = SatellitePlanesRadialPlotter()
    baca = SatellitePlanesBACAPlotter()
    planesize = SatellitePlanesPlanesizePlotter()
    faceon = SatellitePlanesFaceOnPlotter()
    edgeon = SatellitePlanesEdgeOnPlotter()

#    fig,axarr = plt.subplots(2,2,figsize=(12,12))
#    axlist = np.ravel(axarr)
#    hids = haloutils.cid2hid.values(); hids.remove(94687)
#    for ax,sam in zip(axlist,global_plot_sams):
#        ax.set_title(sam)
#        thisplotter = SatellitePlanesBACAPlotter(plot_sams=[sam])
#        caterpillarplot.stackplot(hids,14,thisplotter,ax=ax,labelconc=False,color='k')
#    fig.savefig('5-28/baca_stack.png')

#    fig = caterpillarplot.plot_5x5(faceon)
#    fig.savefig('5-28/faceon.png')
#
#    fig = caterpillarplot.plot_5x5(edgeon)
#    fig.savefig('5-28/edgeon.png')
#
    fig = caterpillarplot.plot_5x5(planesize)
    fig.axes[4].legend(fontsize='x-small',loc='upper right')
    fig.savefig('5-28/planesize.png')
#    
#    fig = caterpillarplot.plot_5x5(baca)
#    fig.axes[4].legend(fontsize='x-small',loc='upper right')
#    fig.savefig('5-28/baca.png')
#
#    fig = caterpillarplot.plot_5x5(angle)
#    fig.axes[4].legend(fontsize='x-small',loc='upper right')
#    fig.savefig('5-28/angle.png')
#
#    fig = caterpillarplot.plot_5x5(radial)
#    fig.axes[4].legend(fontsize='x-small',loc='lower right')
#    fig.savefig('5-28/radial.png')

