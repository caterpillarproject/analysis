import numpy as np
import pylab as plt
import os,subprocess,sys,time
import asciitable

import haloutils
from caterpillaranalysis import PluginBase

class SubPhaseRadialPlugin(PluginBase):
    def __init__(self):
        super(SubPhaseRadialPlugin,self).__init__()
        self.filename='subphaseradial.npz'
        self.nbinsx = 20
        self.nbinsy = 20

        self.logMmin = 6.; self.logMmax = 11.
        self.cmap = plt.get_cmap('cubehelix')

        self.xmin = 0; self.xmax = 300.
        self.ymin = -600.; self.ymax = 600.

        self.n_xmin = 0; self.n_xmax = 1.0
        self.n_ymin = -3;  self.n_ymax = 3
        self.n_xlabel = r'$r/r_{vir}$'
        self.n_ylabel = r'$v_r/v_{vir}$'
        self.xlog = False; self.ylog = False
        self.n_xlog = False; self.n_ylog = False
        self.autofigname = 'subphaserad'
        
    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No Rockstar")
        zoomid = haloutils.load_zoomid(hpath)
        rscat = haloutils.load_rscat(hpath,haloutils.get_numsnaps(hpath)-1)
        hpos = np.array(rscat.ix[zoomid][['posX','posY','posZ']])
        hvel = np.array(rscat.ix[zoomid][['pecVX','pecVY','pecVZ']])
        G = 1.326*10**11 # in km^3/s^2/Msun
        rvir = rscat.ix[zoomid]['rvir'] #kpc/h
        rvirkm = rvir * 3.086*10**16 #km/h
        vvir = np.sqrt(2.*G*rscat.ix[zoomid]['mgrav']/rvirkm) #km/s

        subs = self.get_rssubs(rscat,zoomid)
        spos = np.array(subs[['posX','posY','posZ']])
        svel = np.array(subs[['pecVX','pecVY','pecVZ']]) #pecVX, corevelx, bulkvelx
        smass= np.array(subs['mgrav'])/rscat.h0 #Msun
        dist = 1000.*self.distance(spos,hpos)/rscat.h0 #kpc
        unit_r = self.row_norm(spos-hpos)
        vrad = self.row_dot(svel-hvel, unit_r).T[0] #km/s, not accounting for hubble flow
        vtan = np.sqrt(np.sum((svel-hvel)**2,1)-vrad**2)

        np.savez(self.get_outfname(hpath),
                 rvir=rvir/rscat.h0, vvir=vvir,
                 spos=spos,svel=svel,smass=smass,
                 dist=dist,vrad=vrad,vtan=vtan)
    def _read(self,hpath):
        thisfilename = self.get_filename(hpath)
        data = np.load(thisfilename)
        rvir=data['rvir']; vvir=data['vvir']
        spos=data['spos']; svel=data['svel']; smass=data['smass']
        dist=data['dist']; vrad=data['vrad']; vtan=data['vtan']
        return spos,svel,smass,dist,vrad,vtan,rvir,vvir
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,colorbar=False,**kwargs):
        if colorbar:
            raise NotImplementedError
        spos,svel,smass,dist,vrad,vtan,rvir,vvir = data
        if normtohost:
            dist = dist/rvir #dist and rvir are both units of kpc
            vrad = vrad/vvir #vrad and vvir are both units of km/s
            vtan = vtan/vvir
        
        logsmass = np.log10(smass)
        ii = (logsmass > self.logMmin) & (logsmass < self.logMmax)
        dist = dist[ii]; vrad=vrad[ii]; logsmass=logsmass[ii]
        sc = ax.scatter(dist[::-1],vrad[::-1],s=(logsmass[::-1]/6)**10,
                        cmap=self.cmap,c=logsmass[::-1],vmin=self.logMmin,vmax=self.logMmax,
                        **kwargs)

if __name__=='__main__':
    from caterpillarplot import get_haloidlist
    hids = np.concatenate((get_haloidlist(1),get_haloidlist(2)))
    plug = SubPhaseRadialPlugin()

    all_dist = np.array([])
    all_vrad = np.array([])
    all_vtan = np.array([])
    all_mass = np.array([])
    for hid in hids:
        hpath = haloutils.get_hpath_lx(hid,14)
        if hpath==None: continue
        if not haloutils.check_last_rockstar_exists(hpath): continue
        mvir,rvir,vvir = haloutils.load_haloprops(hpath)
        spos,svel,smass,dist,vrad,vtan,rvir,vvir = plug.read(hpath)
        dist = dist/rvir
        vrad = vrad/vvir
        vtan = vtan/vvir
        smass = smass/mvir
        
        all_dist = np.concatenate((all_dist,dist))
        all_vrad = np.concatenate((all_vrad,vrad))
        all_vtan = np.concatenate((all_vtan,vtan))
        all_mass = np.concatenate((all_mass,smass))
    #import cPickle as pickle
    #with open('data_phase.p','w') as f:
    #    pickle.dump([all_dist,all_vrad,all_mass],f)
    #with open('data_phase.p','r') as f:
    #    all_dist,all_vrad,all_mass = pickle.load(f)

    all_mass = np.log10(all_mass)
    ii = np.isfinite(all_mass)
    all_dist = all_dist[ii]
    all_vrad = all_vrad[ii]
    all_vtan = all_vtan[ii]
    all_mass = all_mass[ii]

    import seaborn as sns
    sns.set_context('poster')
    sns.set_style('ticks')
    fig1,axarr1 = plt.subplots(3,2,figsize=(8,11),sharex=True,sharey=True)
    fig1.subplots_adjust(hspace=.25,wspace=.25,top=.95,right=.95)
    fig2,axarr2 = plt.subplots(3,2,figsize=(8,11),sharex=True,sharey=True)
    fig2.subplots_adjust(hspace=.25,wspace=.25,top=.95,right=.95)
    fig3,axarr3 = plt.subplots(3,2,figsize=(8,11),sharex=True,sharey=True)
    fig3.subplots_adjust(hspace=.25,wspace=.25,top=.95,right=.95)
    axlist1 = np.ravel(axarr1)
    axlist2 = np.ravel(axarr2)
    axlist3 = np.ravel(axarr3)

    Mlabel = [r'$\log \frac{{M_{{\rm sub}}}}{{M_{{\rm vir}}}} = {0}$'.format(x) for x in [-6,-5,-4,-3,-2,-1]]
    logMbins = np.array([-8,-5.5,-4.5,-3.5,-2.5,-1.5,0])

    nbinsx = 20
    nbinsy = 20
    for i in range(6):
        logMmin = logMbins[i]; logMmax = logMbins[i+1]
        ii = (all_mass > logMmin) & (all_mass <= logMmax)
        this_dist = all_dist[ii]
        this_vrad = all_vrad[ii]
        this_vtan = all_vtan[ii]
        this_mass = all_mass[ii]
        
        ax1 = axlist1[i]
        ax1.set_xlim((0,1))
        ax1.set_ylim((-3,3))
        ax1.text(.1,2,Mlabel[i]+', N={0}'.format(np.sum(ii)))
        if i % 2 == 0: ax1.set_ylabel(r'$v_r/v_{\rm vir}$')
        if i > 3: ax1.set_xlabel(r'$r/r_{\rm vir}$')

        ax2 = axlist2[i]
        ax2.set_xlim((0,1))
        ax2.set_ylim((0,2.2))
        ax2.text(.1,2,Mlabel[i]+', N={0}'.format(np.sum(ii)))
        if i % 2 == 0: ax2.set_ylabel(r'$|v_t|/v_{\rm vir}$')
        if i > 3: ax2.set_xlabel(r'$r/r_{\rm vir}$')

        ax3 = axlist3[i]
        ax3.set_xlim((-3,3))
        ax3.set_ylim((0,2.2))
        ax3.text(-2.5,2,Mlabel[i]+', N={0}'.format(np.sum(ii)))
        if i % 2 == 0: ax3.set_ylabel(r'$|v_t|/v_{\rm vir}$')
        if i > 3: ax3.set_xlabel(r'$v_r/v_{\rm vir}$')

        X,Y,Z,levels = plug.density_contour(this_dist,this_vrad,nbinsx,nbinsy)
        cs = ax1.contour(X,Y,Z,levels=levels,colors='k')
        if i>2:
            ax1.scatter(this_dist,this_vrad)

        X,Y,Z,levels = plug.density_contour(this_dist,this_vtan,nbinsx,nbinsy)
        cs = ax2.contour(X,Y,Z,levels=levels,colors='k')
        if i>2:
            ax2.scatter(this_dist,this_vtan)
        
        X,Y,Z,levels = plug.density_contour(this_vrad,this_vtan,nbinsx,nbinsy)
        cs = ax3.contour(X,Y,Z,levels=levels,colors='k')
        if i>2:
            ax3.scatter(this_vrad,this_vtan)
        

    #fig1.savefig('subphase_r_vr.png',bbox_inches='tight')
    #fig2.savefig('subphase_r_vt.png',bbox_inches='tight')
    #fig3.savefig('subphase_vr_vt.png',bbox_inches='tight')
    #plt.show()
