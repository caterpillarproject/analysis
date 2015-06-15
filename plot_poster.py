import matplotlib; matplotlib.use('Agg')
import seaborn.apionly as sns
import numpy as np
import pylab as plt
import pandas as pd

from matplotlib.patches import Rectangle

sns.set_context('poster')
sns.set_style('ticks')
sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
colors = [(0.2980392156862745, 0.4470588235294118, 0.6901960784313725),
          (0.3333333333333333, 0.6588235294117647, 0.40784313725490196),
          (0.7686274509803922, 0.3058823529411765, 0.3215686274509804),
          (0.5058823529411764, 0.4470588235294118, 0.6980392156862745),
          (0.8, 0.7254901960784313, 0.4549019607843137),
          (0.39215686274509803, 0.7098039215686275, 0.803921568627451)]
matplotlib.rcParams.update({'axes.color_cycle':colors})
samlabels = [r'$\rm{First\ 11}$',r'$M_* > 10^{5.5}M_\odot$',r'$M_* > 10^{5}M_\odot$']
shortsamlabels = [r'$\rm{First\ 11}$',r'$> 10^{5.5}M_\odot$',r'$> 10^{5}M_\odot$']

import haloutils
import caterpillarplot

from angmom import AngMomCorrelationPlugin,_plot_mollweide_SAM
from plot_galaxy_planes import SatellitePlanesBACAPlotter,SatellitePlanesEdgeOnPlotter,SatellitePlanesRadialPlotter
from galaxy_planes import SatellitePlanesPlugin,PlaneEvolutionPlugin,plane_tabfn
from caterpillaranalysis import MassAccrPlugin
hids = haloutils.cid2hid.values()
g_hids = list(hids); g_hids.remove(94687)

u_hids = [447649,1130025,1387186,581180,1725372,1354437,94638,95289,1422331,94687]
f_hids = [1631506,264569,1725139,5320,581141,1725272,1195448,1599988,196589,1268839]

def plot_5x4(plug,lx=14,figfilename=None,usecatnum=True,share=False,**kwargs):
    if share:
        fig,axarr = plt.subplots(5,4,figsize=(12,14.4),sharex=True,sharey=True)
    else:
        fig,axarr = plt.subplots(5,4,figsize=(12,14.4))
    fig.subplots_adjust(hspace=0,wspace=0)
    if usecatnum: label='catnum'
    else: label=None
    catnums = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,19,20,22,23,24]
    hids = [haloutils.cid2hid[catnum] for catnum in catnums]
    for i,ax in enumerate(np.ravel(axarr)):
        hid = hids[i]
        hpath = haloutils.get_hpath_lx(hid,lx)
        plug.plot(hpath,ax,**kwargs)
        if hpath==None: continue
        plug.label_plot(hpath,ax,label=label,fontsize='x-large')
    for i in range(4):
        for j in range(1,4):
            axarr[i,j].set_xlabel('')
            axarr[i,j].set_ylabel('')
    for i in range(4):
        axarr[i,0].set_xlabel('')
    for j in range(1,4):
        axarr[4,j].set_ylabel('')
    if figfilename != None:
        fig.savefig(figfilename,bbox_inches='tight')
    else:
        plt.show()
    return fig,axarr

def plot_subs_w_angmom():
    plug = AngMomCorrelationPlugin()
    fig,ax = plt.subplots(figsize=(8,8))
    caterpillarplot.stackplot(hids+[1422331,1631506],14,plug,ax=ax,color1=colors[0],color2=colors[1],color3=colors[2])
    fig.savefig('6-5/angmom_corr.png',bbox_inches='tight')
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
    cbar_ax.set_xlabel(r'$\rm{infall\ scale}$')
    fig.savefig('6-5/angmom_mollweide.png',bbox_inches='tight')
    return fig

def plot_sat_baca():
    fig,ax = plt.subplots(figsize=(8,8))
    #plot_sams = ['Np11','L0p1','L1p1']
    plot_sams = ['Np11','L0m1','L1m1']
    markers = ['o','s','^']
    
    ax.add_patch(Rectangle((0,0),1,0.18,facecolor='lightskyblue',alpha=0.5,linewidth=0))
    for i,sam in enumerate(plot_sams):
        plug = SatellitePlanesBACAPlotter(plot_sams=[sam])
        caterpillarplot.stackplot(g_hids,14,plug,color = colors[i],marker=markers[i],labelconc=False,ax=ax,plotMWM31=False)
    import matplotlib.lines as mlines
    l0 = mlines.Line2D([],[],color=colors[0],marker=markers[0],lw=0,label=samlabels[0])
    l1 = mlines.Line2D([],[],color=colors[1],marker=markers[1],lw=0,label=samlabels[1])
    l2 = mlines.Line2D([],[],color=colors[2],marker=markers[2],lw=0,label=samlabels[2])
    ax.legend(handles=[l0,l1,l2],loc='upper left',fontsize='xx-large')


    fig.savefig('6-5/sat_baca.png',bbox_inches='tight')
    return fig

def plot_H5320_thin():
    hid = 5320; hpath = haloutils.get_hpath_lx(hid,14)
    plug = SatellitePlanesPlugin()
    df,satix = plug.read(hpath)
    #row = df.ix['L0p1']
    row = df.ix['L0m1']
    rperp = row['rperp']
    
    fig,ax = plt.subplots(figsize=(8,8))
    fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
    ax.add_patch(Rectangle((-250,-rperp),500,2*rperp,facecolor='gray',alpha=0.3,linewidth=0))
    ax.text(-225,225,r'$\rm{Cat}-5$',fontsize='x-large')
    plug = SatellitePlanesEdgeOnPlotter()
    plug.plot(hpath,ax,plottracks=True)
    ax.axis('off')
    fig.savefig('thincatplanetracks.png',bbox_inches='tight')

    fig,ax = plt.subplots(figsize=(8,8))
    ax.add_patch(Rectangle((-250,-rperp),500,2*rperp,facecolor='gray',alpha=0.3,linewidth=0))
    ax.text(-225,225,r'$\rm{Cat}-5$',fontsize='x-large')
    plug = SatellitePlanesEdgeOnPlotter()
    plug.plot(hpath,ax,plottracks=True)
    ax.text(90,200,r'$c/a={0:.2f}$'.format(row['ca']),fontsize='x-large')
    ax.text(90,175,r'$r_{{\rm \bot,rms}}={0:.1f}\,\rm{{kpc}}$'.format(row['rperp']),fontsize='x-large')
    ax.text(90,150,r'${0}\ \rm{{satellites}}$'.format(int(row['n_th90'])),fontsize='x-large')
    ax.text(90,125,r'${0}\ \rm{{aligned\ at\ }}30^\circ$'.format(int(row['n_th30'])),fontsize='x-large')
    fig.savefig('6-5/H5320_planetracks.png',bbox_inches='tight')

    fig,ax = plt.subplots(figsize=(8,8))
    ax.add_patch(Rectangle((-250,-rperp),500,2*rperp,facecolor='gray',alpha=0.3,linewidth=0))
    ax.text(-225,225,r'$\rm{Cat}-5$',fontsize='x-large')
    plug = SatellitePlanesEdgeOnPlotter()
    plug.plot(hpath,ax,plottracks=False)
    ax.text(90,200,r'$c/a={0:.2f}$'.format(row['ca']),fontsize='x-large')
    ax.text(90,175,r'$r_{{\rm \bot,rms}}={0:.1f}\,\rm{{kpc}}$'.format(row['rperp']),fontsize='x-large')
    ax.text(90,150,r'${0}\ \rm{{satellites}}$'.format(int(row['n_th90'])),fontsize='x-large')
    ax.text(90,125,r'${0}\ \rm{{aligned\ at\ }}30^\circ$'.format(int(row['n_th30'])),fontsize='x-large')
    fig.savefig('6-5/H5320_notracks.png',bbox_inches='tight')
    return fig
    
def plot_baca_time():
    class PlotPlaneEvol(PlaneEvolutionPlugin):
        def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,
                  color='k',**kwargs):
            assert lx==14 or lx==None
            all_snaps,all_ba,all_ca,all_rperp,all_rpar,all_satpos,all_infall_snaps = data
            t = haloutils.get_t_snap(hpath,all_snaps)
            max_infall_time = haloutils.get_t_snap(hpath,np.max(all_infall_snaps))
            ax.plot(t,all_ca,'-',color=color,**kwargs)
            ax.plot(t,all_ba,'--',color=color,**kwargs)
            ax.plot([self.xmin,self.xmax],[0.2,0.2],'k:')
            ax.plot([max_infall_time,max_infall_time],[self.ymin,self.ymax],':',color=color)

            plug = MassAccrPlugin()
            a_lastMM = plug.read(hpath)[-1]['scale_of_last_MM']
            import brendanlib.conversions as bconversions
            t_lastMM = bconversions.GetTime(a_lastMM,OmegaM=.3175,OmegaL=.6825,h=.6711)
            ax.plot([t_lastMM,t_lastMM],[self.ymin,self.ymax],'k-.')
    #plotter0 = PlotPlaneEvol('L0p1')
    #plotter1 = PlotPlaneEvol('L1p1')
    plotter0 = PlotPlaneEvol('L0m1')
    plotter1 = PlotPlaneEvol('L1m1')

    fig,axarr = plt.subplots(5,4,figsize=(12,14.4))
    fig.subplots_adjust(hspace=0,wspace=0)
    catnums = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,19,20,22,23,24]
    hids = [haloutils.cid2hid[catnum] for catnum in catnums]
    for i,ax in enumerate(np.ravel(axarr)):
        hid = hids[i]
        hpath = haloutils.get_hpath_lx(hid,14)
        plotter0.plot(hpath,ax,color=colors[1])
        plotter1.plot(hpath,ax,color=colors[2])
        plotter0.label_plot(hpath,ax,label='catnum',fontsize='x-large')
    for i in range(4):
        for j in range(1,4):
            axarr[i,j].set_xlabel('')
            axarr[i,j].set_ylabel('')
    for i in range(4):
        axarr[i,0].set_xlabel('')
    for j in range(1,4):
        axarr[4,j].set_ylabel('')

    for i in range(5):
        for j in range(4):
            ax = axarr[i,j]
            ax.set_yticks([.2,.4,.6,.8,1])
            #ax.set_xticks([7,8,9,10,11,12,13])
            ax.set_xticks([0,2,4,6,8,10,12])
            if j>0: ax.set_yticklabels(['','','','',''])
            if i<4: ax.set_xticklabels(['','','','','','',''])
    for j in range(4):
        axarr[4,j].set_yticks([0,.2,.4,.6,.8,1])
        if j==0: axarr[4,j].set_yticklabels(['0.0','0.2','0.4','0.6','0.8','1.0'])
    fig.savefig('6-5/baca_time3.png',bbox_inches='tight')
    return fig

def plot_radial_grid():
    plug = SatellitePlanesRadialPlotter()
    fig,axarr = plot_5x4(plug,colorarr=colors[0:3],markerarr=['o','s','^'],samlabels=shortsamlabels)
    for i in range(5):
        for j in range(4):
            ax = axarr[i,j]
            ax.set_yticks([.2,.4,.6,.8,1])
            ax.set_xticks([0,50,100,150,200,250])
            if j>0: ax.set_yticklabels(['','','','',''])
            if i<4: ax.set_xticklabels(['','','','','',''])
    for j in range(4):
        axarr[4,j].set_yticks([0,.2,.4,.6,.8,1])
        if j==0: axarr[4,j].set_yticklabels(['0.0','0.2','0.4','0.6','0.8','1.0'])
    #axarr[0,0].legend(loc='lower right',fontsize='small')
    fig.savefig('6-5/radial_grid.png',bbox_inches='tight')
    return fig

def plot_host_ca():
    fig,axarr = plt.subplots(2,2,figsize=(12,12),sharey=True)
    df = haloutils.tabulate(plane_tabfn,numprocs=8,exclude_hids=[94687],savefile='plane_table.csv')
    df = pd.read_csv('plane_table.csv',index_col=0)

    #MW sorted by D_MW
    MWdat = np.array([
            [3 , 0.058 ], [4 , 0.051 ], [5 , 0.041 ], [6 , 0.175 ], [7 , 0.241 ], 
            [8 , 0.224 ], [9 , 0.274 ], [10 , 0.221 ], [11 , 0.183 ]])
    #MW sorted by VMag
    MWdat2 = np.array([
        [  3, 0.058 ], [  4, 0.134 ], [  5, 0.161 ], [  6, 0.158 ], [  7, 0.121 ],
        [  8, 0.123 ], [  9, 0.144 ], [ 10, 0.172 ], [ 11, 0.183 ], [ 12, 0.173 ]])
    
    MWdat3 = np.array([
            [3 , 0.058 ], [4 , 0.134 ], [5 , 0.161 ], [6 , 0.158 ], [7 , 0.121 ],
            [8 , 0.123 ], [9 , 0.144 ], [10 , 0.162 ], [11 , 0.183 ], [12 , 0.173 ],
            [13 , 0.214 ], [14 , 0.217 ], [15 , 0.233 ], [16 , 0.232 ], [17 , 0.256 ],
            [18 , 0.258 ], [19 , 0.268 ], [20 , 0.258 ], [21 , 0.264 ], [22 , 0.263 ],
            [23 , 0.267 ], [24 , 0.267 ], [25 , 0.275 ], [26 , 0.276 ]])
    axarr[1,1].plot(MWdat3[:,0],MWdat3[:,1],'k--',label='MW')

    #plot_sams = ['Np11','L0p1','L1p1']
    plot_sams = ['Np11','L0m1','L1m1']
    labels = samlabels
    markers = ['o','s','^']
    mycolors = colors[0:3]
    s = 50
    for sam,color,marker,label in zip(plot_sams,mycolors,markers,labels):
        ax = axarr[0,0]
        ax.set_xlabel(r'$\rm{log\ host\ concentration}$')
        ax.set_ylabel(r'$c/a$')
        ax.scatter(np.log10(df['conc']),df[sam+'_ca'],color=color,marker=marker,s=s)
        ax.set_xlim((.8,1.2)); ax.set_ylim((0,1))
        
        ax = axarr[0,1]
        #ax.set_xlabel(r'$\rm{host\ spin}$')
        #ax.scatter(df['spin'],df[sam+'_ca'],color=color,marker=marker,s=s)
        ax.set_xlabel(r'$\rm{log\ host\ mass}$')
        ax.scatter(np.log10(df['mass']),df[sam+'_ca'],color=color,marker=marker,s=s,label=label)
        #ax.set_xlabel(r'$\rm{host}\ c/a$')
        #ax.scatter(df['c_to_a'],df[sam+'_ca'],color=color,marker=marker,s=s)
        
        ax = axarr[1,0]
        ax.set_ylabel(r'$c/a$')
        ax.set_xlabel(r'$\rm{scale\ of\ last\ major\ merger}$')
        ax.scatter(df['scale_of_last_MM'],df[sam+'_ca'],color=color,marker=marker,s=s)
        ax.set_xlim((0,1))
        
        ax = axarr[1,1]
        ax.set_xlabel(r'$\rm{number\ satellites}$')
        ax.scatter(df[sam+'_Nsats'],df[sam+'_ca'],color=color,marker=marker,s=s)

    axarr[0,1].legend(loc='upper left')
    axarr[1,1].legend(loc='upper left')
    fig.savefig('6-5/host_ca.png',bbox_inches='tight')

def plot_host_cordev(theta_deg=30):
    assert theta_deg in [30,45]

    fig,axarr = plt.subplots(2,2,figsize=(12,12),sharey=True)
    #df = haloutils.tabulate(plane_tabfn,numprocs=8,savefile='plane_table.csv')
    df = pd.read_csv('plane_table.csv',index_col=0)


    color = 'k'
    marker = 'o'
    s = 50

    ax = axarr[0,0]
    ax.set_xlabel(r'$\rm{log\ host\ concentration}$')
    ax.set_ylabel(r'$\rm{correlation\ deviation}$')
    ax.scatter(np.log10(df['conc']),df['corr_enhance'+str(theta_deg)],color=color,marker=marker,s=s)
    ax.set_xlim((.8,1.2))#; ax.set_ylim((0,1))
    
    ax = axarr[0,1]
    ax.set_xlabel(r'$\rm{log\ host\ mass}$')
    ax.scatter(np.log10(df['mass']),df['corr_enhance'+str(theta_deg)],color=color,marker=marker,s=s)
        #ax.set_xlabel(r'$\rm{host}\ c/a$')
        #ax.scatter(df['c_to_a'],df['corr_enhance'+str(theta_deg)],color=color,marker=marker,s=s)
    
    ax = axarr[1,0]
    ax.set_ylabel(r'$\rm{correlation\ deviation}$')
    ax.set_xlabel(r'$\rm{scale\ of\ last\ major\ merger}$')
    ax.scatter(df['scale_of_last_MM'],df['corr_enhance'+str(theta_deg)],color=color,marker=marker,s=s)
    ax.set_xlim((0,1))
    
    ax = axarr[1,1]
    ax.set_xlabel(r'$\rm{host\ spin}$')
    ax.scatter(df['spin'],df['corr_enhance'+str(theta_deg)],color=color,marker=marker,s=s)
    fig.savefig('6-5/host_cordev'+str(theta_deg)+'.png',bbox_inches='tight')

def plot_5x4_all():
    from plot_galaxy_planes import SatellitePlanesAnglePlotter,SatellitePlanesFaceOnPlotter,SatellitePlanesEdgeOnPlotter

    angle = SatellitePlanesAnglePlotter()
    plot_5x4(angle,share=True,figfilename='angle_grid.png')

    edgeon = SatellitePlanesEdgeOnPlotter()
    plot_5x4(edgeon,share=True,figfilename='edgeon_grid.png')
    
    faceon = SatellitePlanesFaceOnPlotter()
    plot_5x4(faceon,share=True,figfilename='faceon_grid.png')
    
    

if __name__=="__main__":
    plot_5x4_all()

    fig = plot_subs_w_angmom()
    fig = plot_two_mollweide_angmom()
    fig = plot_sat_baca()
    fig = plot_H5320_thin()
    fig = plot_baca_time()
    fig = plot_radial_grid()
    fig = plot_host_ca()
    fig = plot_host_cordev(30)
    fig = plot_host_cordev(45)
    #plt.show()
