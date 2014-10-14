import numpy as np
import pylab as plt
import haloutils
import asciitable
import os

def plot_halo_contamdist(hid,ax,pzindex,**kwargs):
    thisindex = pzindex[pzindex['parentid']==hid]
    lxlist = thisindex['LX']
    mindist = thisindex['min2']
    ax.plot(lxlist,mindist,'o-',label=haloutils.hidstr(hid),**kwargs)

if __name__=="__main__":
    pzindex = haloutils.get_parent_zoom_index()
    hidlist = [95289]
    fig,ax = plt.subplots()
    for hid in hidlist:
        plot_halo_contamdist(hid,ax,pzindex)
    ax.plot([10.5,14.5],[1.0,1.0],'k:')
    ax.set_xlim([10.5,14.5])
    ax.set_xticks([11,12,13,14])
    ax.legend(loc='best')
    ax.set_ylabel('distance to parttype2 [Mpc/h]')
    ax.set_xlabel('LX')
    fig.savefig('contamdistres.png',bbox_inches='tight')
    plt.close("all")
