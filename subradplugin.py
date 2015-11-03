import numpy as np
import pylab as plt
import os,subprocess,sys,time
import asciitable
import readsnapshots.readsnapHDF5_greg as rsg
import haloutils

from caterpillaranalysis import *

class SubhaloRadialPlugin(PluginBase):
    def __init__(self,rmin=1,rmax=1000,ymin=10**-2,ymax=10**2):
        super(SubhaloRadialPlugin,self).__init__()
        #self.filename='subradial.dat'
        self.filename='subradial2.dat'

        self.xmin = rmin; self.xmax = rmax
        self.ymin = ymin; self.ymax = ymax
        self.xlabel = r'$r\ (kpc)$'
        self.ylabel = r'$n(r)/<n>$'
        self.xlog = True; self.ylog = True
        self.autofigname = 'subradial'
    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No Rockstar")
        zoomid = haloutils.load_zoomid(hpath)
        rscat = haloutils.load_rscat(hpath,haloutils.get_numsnaps(hpath)-1)
        hpos = np.array(rscat.ix[zoomid][['posX','posY','posZ']])
        subs = self.get_rssubs(rscat,zoomid)
        subs['directsub'] = (subs['hostID'] == zoomid) #new in 2
        spos = np.array(subs[['posX','posY','posZ']])
        dist = 1000*self.distance(spos,hpos)/rscat.h0 #kpc
        iisort = np.argsort(dist)
        sortedids = subs.index[iisort]
        sorteddist= dist[iisort]
        submass = np.array(subs.ix[sortedids]['mvir'])
        submgrav= np.array(subs.ix[sortedids]['mgrav'])
        subrvir = np.array(subs.ix[sortedids]['rvir'])
        subofhost = np.array(subs.ix[sortedids]['directsub']).astype(int) #new in 2
        
        names = ['id','dist','mass','mgrav','rvir','directsub']
        outdict=dict(zip(names,[sortedids,sorteddist,submass,submgrav,subrvir,subofhost]))
        asciitable.write(outdict,self.get_outfname(hpath),names=names)
    def _read(self,hpath):
        thisfilename = self.get_filename(hpath)
        tab = asciitable.read(thisfilename,header_start=0)
        return tab['id'],tab['dist'],tab['mass'],tab['mgrav'],tab['rvir'],tab['directsub'] #new in 2
    def get_Varr(self,rarr):
        rarr = np.concatenate(([0],rarr))
        return 4*np.pi/3 * (rarr[1:]**3 - rarr[:-1]**3)
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        # TODO refactor histogramming into read?
        subid,dist,submass,submgrav,subrvir,directsub = data
        ii = submgrav/submass > 0.0
        dist = dist[ii]; submass = submass[ii]; submgrav = submgrav[ii]; subrvir = subrvir[ii]

        rarr = np.logspace(-2,3,25)
        Varr = self.get_Varr(rarr)
        h_r, x_r = np.histogram(dist, bins=np.concatenate(([0],rarr)))
        # units are number/kpc^3
        n_of_r = h_r/Varr
        mvir,rvir,vvir=haloutils.load_haloprops(hpath)
        n_mean = float(len(dist))/(4*np.pi/3. * rvir**3)
        n_ratio= n_of_r/n_mean

        minii = np.min(np.where(n_ratio > 0)[0])
        if normtohost:
            rarr = rarr/rvir
        if lx != None:
            ax.plot(rarr[minii:],n_ratio[minii:],color=self.colordict[lx],**kwargs)
        else:
            ax.plot(rarr[minii:],n_ratio[minii:],**kwargs)
class SubhaloRadialByMassPlugin(SubhaloRadialPlugin):
    def __init__(self):
        super(SubhaloRadialByMassPlugin,self).__init__()
        self.autofigname = 'subradbymass'
        self.logmassbins = np.array([4,5,6,7,8,9,10,11])
        self.numbins     = 7
        self.labels = [r'$'+str(self.logmassbins[i])+r'<logM<'+str(self.logmassbins[i+1])+r'$' for i in range(self.numbins)]
        self.bincolors = ['b','r','g','y','c','m','k']
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,legendon=True,**kwargs):
        assert lx != None
        subid,dist,submass,submgrav,subrvir,directsub = data
        ii = submgrav/submass > 0.0
        dist = dist[ii]; submass = submass[ii]; submgrav = submgrav[ii]; subrvir = subrvir[ii]

        rarr = np.logspace(-2,3,25)
        Varr = self.get_Varr(rarr)
        mvir,rvir,vvir=haloutils.load_haloprops(hpath)
        
        for i in range(self.numbins):
            binmin = 10.**self.logmassbins[i]
            binmax = 10.**self.logmassbins[i+1]
            ii = (submgrav > binmin) & (submgrav <= binmax)
            h_r, x_r = np.histogram(dist[ii], bins=np.concatenate(([0],rarr)))
            n_of_r = h_r/Varr
            n_mean = np.sum(h_r)/(4*np.pi/3. * rvir**3)
            n_ratio= n_of_r/n_mean
            if normtohost:
                rarr = rarr/rvir
            try: minii = np.min(np.where(n_ratio > 0)[0])
            except ValueError: continue
            ax.plot(rarr[minii:],n_ratio[minii:],color=self.bincolors[i],label=self.labels[i],**kwargs)
        if legendon: ax.legend(loc='lower left',fontsize='xx-small')
class IntegrableSubhaloRadialPlugin(SubhaloRadialPlugin):
    def __init__(self):
        super(IntegrableSubhaloRadialPlugin,self).__init__()
        self.ymin = 0;   self.ymax = 0.3
        self.n_ymin = 0; self.n_ymax = 0.3
        self.ylabel = r'$df/dlogr$'
        self.ylog=False
        self.autofigname = 'integrablesubrad'
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        subid,dist,submass,submgrav,subrvir,directsub = data
        ii = submgrav/submass > 0.0
        dist = dist[ii]; submass = submass[ii]; submgrav = submgrav[ii]; subrvir = subrvir[ii]

        rarr = np.logspace(-2,3,25)
        Varr = self.get_Varr(rarr)
        h_r, x_r = np.histogram(dist, bins=np.concatenate(([0],rarr)))
        n_of_r = h_r/Varr
        n_plot = n_of_r*rarr
        n_plot = n_plot/np.sum(n_plot)
        minii = np.min(np.where(n_plot > 0)[0])
        if normtohost:
            rarr = rarr/rvir
        if lx != None:
            ax.plot(rarr[minii:],n_plot[minii:],color=self.colordict[lx],**kwargs)
        else:
            ax.plot(rarr[minii:],n_plot[minii:],**kwargs)
class IntegrableSubhaloRadialByMassPlugin(SubhaloRadialByMassPlugin):
    def __init__(self):
        super(IntegrableSubhaloRadialByMassPlugin,self).__init__()
        self.ymin = 0;   self.ymax = 1.1
        self.n_ymin = 0; self.n_ymax = 1.1
        self.ylabel = r'$df/dlogr$'
        self.ylog=False
        self.autofigname = 'integrablesubradbymass'
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,legendon=True,**kwargs):
        assert lx != None
        subid,dist,submass,submgrav,subrvir,directsub = data
        ii = submgrav/submass > 0.0
        dist = dist[ii]; submass = submass[ii]; submgrav = submgrav[ii]; subrvir = subrvir[ii]

        rarr = np.logspace(-2,3,25)
        Varr = self.get_Varr(rarr)
        mvir,rvir,vvir=haloutils.load_haloprops(hpath)
        
        for i in range(self.numbins):
            binmin = 10.**self.logmassbins[i]
            binmax = 10.**self.logmassbins[i+1]
            ii = (submgrav > binmin) & (submgrav <= binmax)
            h_r, x_r = np.histogram(dist[ii], bins=np.concatenate(([0],rarr)))
            n_of_r = h_r/Varr
            n_plot = n_of_r*rarr
            n_plot = n_plot/np.sum(n_plot)
            try: minii = np.min(np.where(n_plot > 0)[0])
            except ValueError: continue
            if normtohost:
                rarr = rarr/rvir
            ax.plot(rarr[minii:],n_plot[minii:],color=self.bincolors[i],label=self.labels[i],**kwargs)
        if legendon: ax.legend(loc='lower left',fontsize='xx-small')
    
class SubhaloRadialSubmassFracPlugin(MultiPlugin):
    def __init__(self,rmin=10**0,rmax=10**3):
        allplug = ProfilePlugin()
        subplug = SubhaloRadialPlugin()
        super(SubhaloRadialSubmassFracPlugin,self).__init__([allplug,subplug])
        
        self.xmin = rmin; self.xmax = rmax
        self.ymin = 10**-4; self.ymax = 10**0
        self.xlabel = r'$r\ (kpc)$'
        self.ylabel = r'$f_{m,\ \rm{subs}}$'
        #self.ylabel = r'$\rm{Mass\ fraction\ in\ sub}$'
        self.xlog = True; self.ylog = True
        self.autofigname = 'subradialsubmassfrac'
    def _plot(self,hpath,datalist,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        if normtohost:
            raise NotImplementedError
        alldata = datalist[0]
        subdata = datalist[1]
        r,mltr,p03r,rvir,r200c,pNFW,pEIN = alldata        
        rbin = np.concatenate(([0],r))
        subid,dist,submass,submgrav,subrvir,directsub = subdata
        ii = np.logical_and(submgrav/submass > 0.0, directsub) #NEW in 2
        dist = dist[ii]; submass = submass[ii]; submgrav = submgrav[ii]; subrvir = subrvir[ii]

        dist /= 1000. #kpc to Mpc
        m_r, x_r = np.histogram(dist, weights = submgrav, bins=rbin)
        mltrsub = np.cumsum(m_r)

        r = r*1000
        iigood = mltr>0
        mfrac = mltrsub/mltr
        assert np.all(mfrac[iigood] >= 0),"Max: %f, Min: %f" % (np.max(mfrac[iigood]),np.min(mfrac[iigood]))
        eps = 1000*haloutils.load_soft(hpath)
        ii1 = r >= eps
        if lx != None:
            color = self.colordict[lx]
            ax.plot(r[ii1], mfrac[ii1], color=color, **kwargs)
        else:
            ax.plot(r[ii1], mfrac[ii1], **kwargs)

class SubhaloRadialMassPlugin(ProfilePlugin):
    def __init__(self,rmin=10**-2,rmax=10**3):
        print "WARNING: this plugin probably doesn't give useful results right now due to using all (bound+unbound) particles"
        super(SubhaloRadialMassPlugin,self).__init__()
        self.filename='subradialmass.npz'

        self.xmin = rmin; self.xmax = rmax
        self.ymin = 10**4; self.ymax = 10**13 #Msun
        self.xlabel = r'$r\ (kpc)$'
        self.ylabel = r'$Substructure\ enclosed\ mass$'
        self.xlog = True; self.ylog = True
        self.autofigname = 'subradialmass'
    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No Rockstar")
        snap = haloutils.get_numsnaps(hpath)-1
        snapstr = str(snap).zfill(3)
        snapfile = hpath+'/outputs/snapdir_'+snapstr+'/snap_'+snapstr
        header = rsg.snapshot_header(snapfile+'.0')
        zoomid = haloutils.load_zoomid(hpath)
        rscat = haloutils.load_rscat(hpath,snap)
        subs = self.get_rssubs(rscat,zoomid)
        subids = np.array(subs['id'])
        rarr = self.get_rarr()
        rarr,mltrarr,haloparts = self.compute_subpart_profile(rarr,hpath,rscat,subids,snap,header)
        np.savez(self.get_outfname(hpath),rarr=rarr,mltr=mltrarr,haloparts=haloparts)
    def compute_one_profile(self,rarr,hpath,rscat,rsid,snap,header,calcp03r=True,calcr200=True):
        raise NotImplementedError
    def compute_subpart_profile(self,rarr,hpath,rscat,subids,snap,header,calcp03r=False,calcr200=False):
        zoomid = haloutils.load_zoomid(hpath)
        subparts = np.array([])
        for subid in subids:
            subparts = np.concatenate((subparts,rscat.get_all_particles_from_halo(subid)))
            #subparts = np.union1d(subparts,rscat.get_all_particles_from_halo(subid))
        subparts = subparts.astype(int)
        subparts = np.unique(subparts)
        halopos = np.array(rscat.ix[zoomid][['posX','posY','posZ']])
        try:
            #subparts = np.sort(subparts)
            partpos = haloutils.load_partblock(hpath,snap,"POS ",parttype=1,ids=subparts)
        except IndexError as e:
            print e
            raise RuntimeError("Contamination in halo")
        dr = np.sort(self.distance(partpos,halopos))/header.hubble #Mpc
        mltrarr,p03rmin,r200c = self.calc_mltr_radii(rarr,dr,header,subparts,
                                                     calcp03r=calcp03r,calcr200=calcr200)
        return rarr,mltrarr,subparts #all in physical units
    def _read(self,hpath):
        thisfilename = self.get_filename(hpath)
        data = np.load(thisfilename)
        r = data['rarr']
        mltr = data['mltr']
        return r,mltr
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        if normtohost:
            raise NotImplementedError
        r,mltr = data
        r = r*1000. # kpc
        eps = 1000*haloutils.load_soft(hpath)
        ii1 = r >= eps
        if lx != None:
            color = self.colordict[lx]
            ax.plot(r[ii1],mltr[ii1],color=color,**kwargs)
        else:
            ax.plot(r[ii1],mltr[ii1],**kwargs)

class SubhaloRadialMassFracPlugin(MultiPlugin):
    def __init__(self,rmin=10**-2,rmax=10**3):
        print "WARNING: this is not the right plugin for substructure mass (use SubhaloRadialSubmassFracPlugin)"
        allplug = ProfilePlugin()
        subplug = SubhaloRadialMassPlugin()
        super(SubhaloRadialMassFracPlugin,self).__init__([allplug,subplug])
        
        self.xmin = rmin; self.xmax = rmax
        self.ymin = 10**-4; self.ymax = 10**0 #Msun?
        self.xlabel = r'$r\ (kpc)$'
        self.ylabel = r'$f_{m,\ \rm{subs}}$'
        self.xlog = True; self.ylog = True
        self.autofigname = 'subradialmassfrac'
    def _plot(self,hpath,datalist,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        if normtohost:
            raise NotImplementedError
        alldata = datalist[0]
        subdata = datalist[1]
        r,mltr,p03r,rvir,r200c = alldata        
        rsub,mltrsub = subdata
        assert np.sum(np.abs(r-rsub)) < 10**-9
        r = r*1000
        iigood = mltr>0
        mfrac = mltrsub/mltr
        assert np.all((mfrac[iigood] <= 1) & (mfrac[iigood] >= 0)),"Max: %f, Min: %f" % (np.max(mfrac),np.min(mfrac))
        eps = 1000*haloutils.load_soft(hpath)
        ii1 = r >= eps
        if lx != None:
            color = self.colordict[lx]
            ax.plot(r[ii1], mfrac[ii1], color=color, lw=1, **kwargs)
        else:
            ax.plot(r[ii1], mfrac[ii1], lw=1, **kwargs)
