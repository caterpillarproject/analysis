import numpy as np
import pylab as plt
from caterpillaranalysis import *
from caterpillarplot import *
from subradplugin import *

import haloutils

def old():
    params = {'axes.labelsize': 22,
              'text.fontsize': 22,
              'legend.fontsize': 14,
              'xtick.labelsize': 18,
              'ytick.labelsize': 18,
              #'backend': 'ps',
              'text.usetex': True}
    
    plt.rcParams.update(params)

    fig,axarr = plt.subplots(3,2,figsize=(16,24))
    fig.subplots_adjust(hspace=.3,wspace=.25)
    axNvmax = axarr[0,0]
    axSHMF  = axarr[0,1]
    axLinMassAccr = axarr[1,0]
    axSubMassFrac = axarr[1,1]
    axSubRadProf  = axarr[2,0]
    axHostProf = axarr[2,1]

    plugNvmax = NvmaxPlugin()
    plugSHMF  = SHMFPlugin()
    plugLinMassAccr = LinearMassAccrPlugin()
    plugSubMassFrac = SubhaloRadialSubmassFracPlugin()
    plugSubRadProf  = SubhaloRadialPlugin()
    plugHostProf = ProfilePlugin()

    lx=14; lw=2
    paper_stackplot(lx,plugNvmax,ax=axNvmax,lw=lw)
    paper_stackplot(lx,plugSHMF,ax=axSHMF,lw=lw)
    paper_stackplot(lx,plugLinMassAccr,ax=axLinMassAccr,lw=lw)
    paper_stackplot(lx,plugSubMassFrac,ax=axSubMassFrac,lw=lw)
    paper_stackplot(lx,plugSubRadProf,ax=axSubRadProf,lw=lw)
    paper_stackplot(lx,plugHostProf,ax=axHostProf)

    colordict = haloutils.get_colors_for_halos()
    labeldict = haloutils.hid2name
    haloids = get_haloidlist(1)
    colors = [colordict[hid] for hid in haloids]
    labels = [labeldict[hid] for hid in haloids]
    axNvmax.legend(axNvmax.lines,labels,loc='lower left',ncol=2,frameon=False)  # ,fontsize='xx-small'

    plt.savefig('paper_panels.png')
    plt.savefig('paper_panels.eps')
    plt.show()

if __name__=="__main__":
    
    params = {'axes.labelsize': 22,
              'text.fontsize': 22,
              'legend.fontsize': 14,
              'xtick.labelsize': 18,
              'ytick.labelsize': 18,
              #'backend': 'ps',
              'text.usetex': True}
    
    plt.rcParams.update(params)

    plugNvmax = NvmaxPlugin()
    plugSHMF  = SHMFPlugin()
    plugLinMassAccr = LinearMassAccrPlugin()
    plugSubMassFrac = SubhaloRadialSubmassFracPlugin()
    plugSubRadProf  = SubhaloRadialPlugin()
    plugHostProf = ProfilePlugin()

    lx=14; lw=3; ext='.png'
    paper_stackplot(lx,plugNvmax,lw=lw,figfilename='p_Nvmax'+ext,legendloc='lower left')
    paper_stackplot(lx,plugSHMF,lw=lw,figfilename='p_SHMF'+ext,legendloc='upper right')
    #paper_stackplot(lx,plugLinMassAccr,lw=lw)
    paper_stackplot(lx,plugSubMassFrac,lw=lw,figfilename='p_SubRadMassFrac'+ext,legendloc='upper left')
    paper_stackplot(lx,plugSubRadProf,lw=lw,figfilename='p_SubRadNum'+ext,legendloc='lower left')
    paper_stackplot(lx,plugHostProf,figfilename='p_HostProf'+ext)

    plt.close('all')
