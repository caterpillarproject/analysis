import numpy as np
import pylab as plt
import matplotlib

import haloutils
from caterpillaranalysis import *
from caterpillarplot import *

if __name__=="__main__":
    sheet = 1
    plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
    myparams = {'text.usetex' : True,
                'font.size' : 16,
                'font.family' : 'lmodern',
                'text.latex.unicode': True,
                'axes.linewidth': 2}
    plt.rcParams.update(myparams)

    ## Converge LX14 density profiles with Einasto fit and alpha
    R2Profile  = R2ProfilePlugin()
    fig = convergeplot(sheet,R2Profile,whichlx=[14],plotEIN=True,labelalpha=True,usehaloname=True)
    axlist = fig.axes
    for ax in axlist:
        ax.tick_params(which='major',axis='both', color='k', length=6,width=2)
        ax.tick_params(which='minor',axis='both', color='k', length=3,width=2)
    fig.savefig('brendan_hostdensity.png',bbox_inches='tight')

    ## Stacked SHMF
    SHMF = SHMFPlugin()
    haloids = haloutils.cid2hid.values()
    fig,ax = plt.subplots(figsize=(8,8))
    xvals = np.logspace(-7.5,2.5)
    yvals = xvals ** -1.9 * 10**-1.7
    ax.plot(xvals,yvals,'k:')
    stackplot(haloids,14,SHMF,ax=ax,normtohost=True,lw=2,color='k',alpha=.3)
    ax.tick_params(which='major',axis='both', color='k', length=6,width=2)
    ax.tick_params(which='minor',axis='both', color='k', length=3,width=2)
    fig.savefig('brendan_stackSHMF.png',bbox_inches='tight')

    ## Converge Nvmax
    Nvmax = NvmaxPlugin()
    fig = convergeplot(sheet,Nvmax,usehaloname=True)
    axlist = fig.axes
    for ax in axlist:
        ax.tick_params(which='major',axis='both', color='k', length=6,width=2)
        ax.tick_params(which='minor',axis='both', color='k', length=3,width=2)
    fig.savefig('brendan_Nvmax.png',bbox_inches='tight')
