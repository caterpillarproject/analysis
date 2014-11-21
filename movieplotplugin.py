import numpy as np
import pylab as plt
import os,subprocess
import haloutils
from analysisplugin import *
import pynbody as pnb
import mergertrees.MTCatalogue as MTC
import gc,sys
from brendanlib.grifflib import makecolormap
import matplotlib

class MoviePlotPlugin(AnalysisPluginBase):
    def __init__(self):
        self.filename = ''
        self.profileplug = ProfilePlugin()
        self.SHMFplug = SHMFPlugin()
        self.Nvmaxplug = NvmaxPlugin()

    def plot_movie_plots(self,hid,lx):
        hpath = haloutils.get_hpath_lx(hid,lx)
        vizdir   = '/spacebase/data/alexji/analysis/visualizations/movies/'+haloutils.hidstr(hid)+'_LX'+str(lx)
        datapath = hpath+'/'+OUTPUTFOLDERNAME
        massaccr = np.load(datapath+'/massaccr.npy')
        profarr  = np.load(datapath+'/profilemovie.npy')
        SHMFarr  = np.load(datapath+'/SHMFmovie.npy')
        Nvmaxarr = np.load(datapath+'/Nvmaxmovie.npy')
        cmap = makecolormap()
        vmin=2.*10**6;vmax=7.*10**8

        rarr  = np.logspace(-5,0,50)
        Mplot = 10.**((self.SHMFplug.histrange[:-1]+self.SHMFplug.histrange[1:])/2.0)
        vplot = self.Nvmaxplug.vmaxbins[1:]

        scales = massaccr[0,:]
        massaccrarr = massaccr[1,:]
        snaparr = np.arange(256-len(scales),256)[::-1]
        # go backwards through snaps since the merger tree goes backwards
        for i,snap in enumerate(snaparr):
            snapstr = str(snap).zfill(3)
            print snapstr
            fig,ax0,ax1,ax2,ax3,ax4 = self._make_figure()
            im = np.load(vizdir+'/im'+snapstr+'.npy')
            try:
                im[np.where(im==0)] = abs(im[np.where(abs(im != 0))]).min()
            except ValueError:
                raise ValueError, "Failed to make a sensible logarithmic image. This probably means there are no particles in the view."
            #import pdb
            #pdb.set_trace()

            if (im < 0).any():
                linthresh = np.nanmax(abs(im)) / 1000.
                norm = matplotlib.colors.SymLogNorm(
                    linthresh, vmin=vmin, vmax=vmax)
            else:
                norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
            #linthresh = np.nanmax(abs(im)) / 1000.
            #norm = matplotlib.colors.SymLogNorm(
            #    linthresh, vmin=vmin, vmax=vmax)
            ax0.imshow(im[::-1,:].view(np.ndarray),cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
            ax0.set_xticklabels([]); ax0.set_yticklabels([])

            prof = profarr[:,i]
            SHMF = SHMFarr[:,i]
            Nvmax= Nvmaxarr[:,i]
            hassubhalos = (Nvmax[0]>0)

            ax1.plot(rarr*1000.,(rarr)**2*prof)
            ax1.set_xscale('log'); ax1.set_yscale('log')
            ax1.set_xlim((10**-2,10**3)); ax1.set_ylim((10**-1.5,10**2.5))
            ax1.set_xlabel('r [kpc]'); ax1.set_ylabel(r'$r^2 \rho$ [$10^{10} M_\odot$ Mpc$^{-1}$]')

            if hassubhalos: ax2.plot(Mplot,SHMF)
            ax2.set_xscale('log'); ax2.set_yscale('log')
            ax2.set_xlim((10**4.5,10**10.6)); ax2.set_ylim((10**-10,10**-1.0))
            ax2.set_xlabel(r'M [$h^{-1} M_\odot$]'); ax2.set_ylabel(r'$dn/dM$')

            ax3.plot(scales[i:],massaccrarr[i:]/(10**12*.6711))
            ax3.set_xlim((0,1)); ax3.set_ylim((0,1))
            ax3.set_xlabel('a'); ax3.set_ylabel(r'M [$10^{12} M_\odot$]')

            if hassubhalos: ax4.plot(vplot,Nvmax)
            ax4.set_xscale('log'); ax4.set_yscale('log')
            ax4.set_xlim((.3,100)); ax4.set_ylim((1,10**4.5))
            ax4.set_xlabel('Vmax [km/s]'); ax4.set_ylabel('N(>Vmax)')
            fig.savefig(vizdir+'/movieplot'+snapstr+'.png')
            plt.close('all')
            
    def _make_figure(self,buf=.07,figsize=(12,6)):
        assert figsize[0] == 2*figsize[1]
        bufx = buf/2.; bufy = buf
        fig = plt.figure(figsize=figsize)
        ax0 = fig.add_axes((bufx,bufy,.5-2*bufx,(1-2*bufy)))
        ax1 = fig.add_axes((0.5+bufx,0.5+bufy,.25-2*bufx,.5-2*bufy))
        ax2 = fig.add_axes((.75+bufx,0.5+bufy,.25-2*bufx,.5-2*bufy))
        ax3 = fig.add_axes((0.5+bufx,0.0+bufy,.25-2*bufx,.5-2*bufy))
        ax4 = fig.add_axes((.75+bufx,0.0+bufy,.25-2*bufx,.5-2*bufy))
        for ax in [ax0,ax1,ax2,ax3,ax4]:
            ax.tick_params(axis='both',which='major',labelsize=8)
        return fig,ax0,ax1,ax2,ax3,ax4

    def calc_movie_plots(self,hid,lx):
        hpath = haloutils.get_hpath_lx(hid,lx)
        mtc = MTC.MTCatalogue(hpath+'/halos/trees',version=4,numHosts=1)
        mt = mtc[0]
        zoomid = haloutils.load_zoomid(hpath)
        if mt.rockstar_id != zoomid:
            raise ValueError()
        mb = mt.getMainBranch()
        numsnaps = len(mb)
        
        ## Density Profile
        numrbins = len(np.logspace(-5,0,50))
        ## SHMF
        SHMFbins = self.SHMFplug.histrange
        ## N(>Vmax)
        vmaxbins = self.Nvmaxplug.vmaxbins
        ## Mass Accr History
        massaccr = np.zeros((2,numsnaps))
        massaccr[0,:] = np.array(mb['scale'])
        massaccr[1,:] = np.array(mb['mvir'])
        np.save(hpath+'/'+OUTPUTFOLDERNAME+'/massaccr',massaccr)

        densityprofarr = np.zeros((numrbins,numsnaps))
        SHMFarr = np.zeros((len(SHMFbins)-1,numsnaps))
        Nvmaxarr = np.zeros((len(vmaxbins)-1,numsnaps))
        
        outdir = 'visualizations/movies/'+haloutils.hidstr(hid)+"_LX"+str(lx)
        for irow,row in enumerate(mb):
            snap = row['snap']
            print snap
            rscat = haloutils.load_rscat(hpath,snap)
            zoomid = row['origid']
            subs = self.get_rssubs(rscat,zoomid,True,None)
            ## Density Profile
            try:
                rarr,rhoarr,p03rmin,halorvir,r200c,halomass = self.profileplug.compute_one_profile(hpath,zoomid=zoomid,snap=snap)
                densityprofarr[:,irow] = rhoarr
            except RuntimeError:
                print "Profile error at snap "+str(snap)
            ## SHMF
            try:
                subm = np.array(subs['mvir'])
                x,shmf = self.SHMFplug.MassFunc_dNdM(subm,SHMFbins)
                SHMFarr[:,irow] = shmf
            except Exception as e:
                print e
                print "SHMF error at snap "+str(snap)                
            ## N(>Vmax)
            try:
                subvmax = np.array(subs['vmax'])
                Nvmaxarr[:,irow] = self.Nvmaxplug.calcNvmax(subvmax)
            except Exception as e:
                print e
                print "N(>Vmax) error at snap "+str(snap)                
            gc.collect()
        np.save(hpath+'/'+OUTPUTFOLDERNAME+'/profilemovie',densityprofarr)
        np.save(hpath+'/'+OUTPUTFOLDERNAME+'/SHMFmovie',SHMFarr)
        np.save(hpath+'/'+OUTPUTFOLDERNAME+'/Nvmaxmovie',Nvmaxarr)

if __name__=="__main__":
    plug = MoviePlotPlugin()
    #plug.calc_movie_plots(95289,12)
    plug.plot_movie_plots(95289,12)
