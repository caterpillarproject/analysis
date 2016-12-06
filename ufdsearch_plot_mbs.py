import glob
import haloutils
import os, sys, time
import cPickle as pickle
from six import iteritems
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd

import MTaddition as mtadd
import MTanalysis3 as mta
AE = mta.AllExtantData()
E = mtadd.ExtantDataReionization()

## Greg's abundance matching
sys.path.insert(0,"/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code")
import abundance_matching as greg_am
## Added infall_times.py to my local directory

if __name__=="__main__":
    with open("UFDSEARCHTMP/H1725139_MBs_7585.p","r") as f:
        ajmbs,gdmbs = pickle.load(f)

    ajkwds = {'color':'k','alpha':0.2,'lw':1.5}
    gdkwds = {'color':'r','alpha':0.2,'lw':1.0}

    fig,axes = plt.subplots(1,2,figsize=(12,6))
    ax = axes[0]
    for mb in ajmbs:
        ax.plot(1./mb['scale']-1,np.log10(mb['mvir']), **ajkwds)
    for mb in gdmbs:
        ax.plot(1./mb['scale']-1,np.log10(mb['mvir']), **gdkwds)
    ax.set_ylabel("log Mvir (mgrav)")

    ax = axes[1]
    vmax_ach = 9.48535156
    vmax_filt= 23.54248047
    for mb in ajmbs:
        ax.plot(1./mb['scale']-1,mb['vmax'], **ajkwds)
    for mb in gdmbs:
        ax.plot(1./mb['scale']-1,mb['vmax'], **gdkwds)
    ax.set_ylabel("Vmax")
    ax.plot([0,20],[vmax_ach,vmax_ach],'b:')
    ax.plot([0,20],[vmax_filt,vmax_filt],'b:')
    
    for ax in axes:
        ax.set_xlabel("z")
        ax.set_xlim((0,20))
        ylim = ax.get_ylim()
        ax.plot([12,12],ylim,'b:')
        ax.plot([8,8],ylim,'g:')
        ax.set_ylim(ylim)

    fig.savefig("fig_ufdsearch_vs_greg.png",bbox_inches='tight')

    plt.show()
    
