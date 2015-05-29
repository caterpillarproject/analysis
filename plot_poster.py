import matplotlib
import seaborn.apionly as sns
import numpy as np
import pylab as plt

sns.set_context('poster')
sns.set_style('ticks')
sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
colors = [(0.2980392156862745, 0.4470588235294118, 0.6901960784313725),
          (0.3333333333333333, 0.6588235294117647, 0.40784313725490196),
          (0.7686274509803922, 0.3058823529411765, 0.3215686274509804),
          (0.5058823529411764, 0.4470588235294118, 0.6980392156862745),
          (0.8, 0.7254901960784313, 0.4549019607843137),
          (0.39215686274509803, 0.7098039215686275, 0.803921568627451)]

import haloutils
import caterpillarplot

from angmom import AngMomCorrelationPlugin,_plot_mollweide_SAM
from plot_galaxy_planes import SatellitePlanesBACAPlotter,SatellitePlanesEdgeOnPlotter
hids = haloutils.cid2hid.values()
g_hids = list(hids); g_hids.remove(94687)

u_hids = [447649,1130025,1387186,581180,1725372,1354437,94638,95289,1422331,94687]
f_hids = [1631506,264569,1725139,5320,581141,1725272,1195448,1599988,196589,1268839]

def plot_subs_w_angmom():
    plug = AngMomCorrelationPlugin()
    fig,ax = plt.subplots(figsize=(8,8))
    caterpillarplot.stackplot(hids+[1422331,1631506],14,plug,ax=ax,color1=colors[0],color2=colors[1],color3=colors[2])
    fig.savefig('POSTER_PLOTS/angmom_corr.png',bbox_inches='tight')
    return fig
def plot_two_mollweide_angmom():
    hid1 = 1631506
    hid2 = 1422331
    hpath1 = haloutils.get_hpath_lx(hid1,14)
    hpath2 = haloutils.get_hpath_lx(hid2,14)
    fig,(ax1,ax2) = plt.subplots(2,1,figsize=(6,8),subplot_kw={'projection':'mollweide'})
    fig.subplots_adjust(bottom=0.2)
    whatdata = 'angmom'; tag = 'J'; logMpeakcut = 6
    sc = _plot_mollweide_SAM(whatdata,hpath1,ax1,tag,logMpeakcut=logMpeakcut)
    _plot_mollweide_SAM(whatdata,hpath2,ax2,tag,logMpeakcut=logMpeakcut)
    ax1.set_title(r'$\rm{Cat}-1$')
    ax2.set_title(r'$\rm{Cat}-22$')

    cbar_ax = fig.add_axes([.125,.1,.75,.06])
    fig.colorbar(sc,cax=cbar_ax,orientation='horizontal' )
    cbar_ax.set_xlabel('infall scale')
    fig.savefig('POSTER_PLOTS/angmom_mollweide.png',bbox_inches='tight')
    return fig

def plot_sat_baca():
    fig,ax = plt.subplots(figsize=(8,8))
    plot_sams = ['Ni11','L0i1','L1i1']
    markers = ['o','s','^']
    for i,sam in enumerate(plot_sams):
        plug = SatellitePlanesBACAPlotter(plot_sams=[sam])
        caterpillarplot.stackplot(g_hids,14,plug,color = colors[i],marker=markers[i],labelconc=False,ax=ax)
    import matplotlib.lines as mlines
    l0 = mlines.Line2D([],[],color=colors[0],marker=markers[0],lw=0,label='A')
    l1 = mlines.Line2D([],[],color=colors[1],marker=markers[1],lw=0,label='B1')
    l2 = mlines.Line2D([],[],color=colors[2],marker=markers[2],lw=0,label='B2')
    ax.legend(handles=[l0,l1,l2],loc='upper left',fontsize='xx-large')
    fig.savefig('POSTER_PLOTS/sat_baca.png',bbox_inches='tight')
    return fig

def plot_H5320_thin():
    #TODO label with c/a, r_perp, r_par, N corotating
    #make symbols larger
    fig,ax = plt.subplots(figsize=(8,8))
    hid = 5320; hpath = haloutils.get_hpath_lx(hid,14)
    plug = SatellitePlanesEdgeOnPlotter()
    plug.plot(hpath,ax)
    return fig
    
if __name__=="__main__":


    #fig = plot_subs_w_angmom()
    #fig = plot_two_mollweide_angmom()
    #fig = plot_sat_baca()
    fig = plot_H5320_thin()
    plt.show()
