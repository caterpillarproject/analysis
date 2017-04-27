import numpy as np
import haloutils
import time,sys
import pandas as pd
import cPickle as pickle
import matplotlib.pyplot as plt

import sys
sys.path.append("./greg_dwarfs")
import DwarfMethods as dm

def get_hostpos(hpath, snap):
    from caterpillaranalysis import MassAccrPlugin
    plug = MassAccrPlugin()
    tab = plug.read(hpath)
    snaplist = tab['snap']
    ii = snaplist == snap
    if np.sum(ii) != 1: 
        raise ValueError("{0} snap {1} does not have a valid main branch rsid ({2} indices match snap)".format(get_foldername(hpath),snap,np.sum(ii)))
    return np.array([tab[ii]['x'],tab[ii]['y'],tab[ii]['z']]).reshape(-1)

def plot_dr(hpath, snap, ax):
    rscat = haloutils.load_rscat(hpath, snap, rmaxcut=False)
    pos = rscat[['posX','posY','posZ']]
    dr1 = np.sqrt(np.sum((pos-50.)**2, axis=1))

    hostpos = get_hostpos(hpath, snap)
    dr2 = np.sqrt(np.sum((pos-hostpos)**2, axis=1))
    
    dr = dr2
    
    # Volume 1: maximum radius
    drmax = np.max(dr)
    Vmax = 4.*np.pi/3 * drmax**3

    # Volume 2: dN/dr
    rbins = np.arange(0,drmax+0.02,0.02)
    rmid = (rbins[1:]+rbins[:-1])/2.
    h, _ = np.histogram(dr2, bins=rbins)
    
    h1, _ = np.histogram(dr1, bins=rbins)
    h2, _ = np.histogram(dr2, bins=rbins)
    ax.plot(rmid, h1, drawstyle='steps-mid', label='(50,50,50)')
    ax.plot(rmid, h2, drawstyle='steps-mid', label='host progenitor pos')
    ax.legend(loc='upper left')
    ax.set_xlabel('r (h^-1 Mpc)')
    ax.set_ylabel('dN/dr')
    
if __name__=='__main__':
    hpaths = dm.get_hpaths(field=False, lx=14)
    hpath = hpaths[5]
    snap = 67 #z_r = 8.346
    
    fig, axes = plt.subplots(6,5,figsize=(20,20))
    for hpath,ax in zip(hpaths[:30],np.ravel(axes)):
        plot_dr(hpath, snap, ax)

    #Ntot = h.sum()
    #plt.plot(rmid, Ntot*(rmid/drmax)**3)
    plt.show()


