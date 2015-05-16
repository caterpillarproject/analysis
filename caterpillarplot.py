import haloutils
import numpy as np
import pylab as plt
import seaborn.apionly as sns

plt.rcParams.update({'text.latex.preamble': r"\usepackage{amsmath}"})
#plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
#plt.rcParams.update({'text.usetex': True,
#                     'font.family': 'lmodern',
#                     'text.latex.unicode': True,
#                     'axes.linewidth': 2})

def get_haloidlist(sheet):
    if sheet==1:
        start = 1; stop = 13
        #haloidlist = [1631506,264569,  1725139,
        #              447649, 5320,    581141,
        #              94687,  1130025, 1387186,
        #              581180, 1725372, 1354437]
    elif sheet==2:
        start = 13; stop = 25
        #haloidlist = [1725272, 1195448, 1599988,
        #              796175,  388476,  1079897, 
        #              94638,   95289,   1232164,
        #              1422331, 196589,  1268839]
    else:
        raise ValueError("Sheet = {0} not valid".format(sheet))
        
    haloidlist = []
    for i in range(start,stop):
        haloidlist.append(haloutils.cid2hid[i])

    assert len(haloidlist) == 12
    return haloidlist

def sheetplot(sheetnum,plug,whichlx=[14],figfilename=None,figsize=None,**kwargs):
    """ Wrapper for convergeplot with whichlx=[14] """
    return convergeplot(sheetnum,plug,whichlx=whichlx,figfilename=figfilename,figsize=figsize,**kwargs)

def convergeplot(sheetnum,plug,whichlx=[14,13,12,11],figfilename=None,figsize=None,**kwargs):
    haloidlist = get_haloidlist(sheetnum)
    if figsize==None: figsize=(9,11)
    fig,axarr = plt.subplots(4,3,figsize=figsize,
                             sharex=True,sharey=True)
    plt.subplots_adjust(wspace=0,hspace=0)
    axlist = axarr.reshape(-1)
    for i,(hid,ax) in enumerate(zip(haloidlist,axlist)):
        plug.lxplot(hid,ax,whichlx=whichlx,**kwargs)
    for ax in np.ravel(axarr[0:2,:]):
        ax.set_xlabel('')
    for ax in np.ravel(axarr[:,1:3]):
        ax.set_ylabel('')
    if figfilename != None:
        fig.savefig(figfilename,bbox_inches='tight')
    else:
        plt.show()
    return fig

def stackplot(haloids,lx,plug,figfilename=None,ax=None,autocolor=None,**kwargs):
    if autocolor != None:
        colors = get_color_palette(autocolor,len(haloids))

    if ax == None: 
        fig,ax = plt.subplots()
        plotfig = True
    else:
        assert figfilename==None,'Cannot specify both ax and figfilename'
        plotfig = False
    for i,hid in enumerate(haloids):
        hpath = haloutils.get_hpath_lx(hid,lx)
        if autocolor == None:
            plug.plot(hpath,ax,**kwargs)
        else:
            plug.plot(hpath,ax,color=colors[i],**kwargs)
    if not plotfig: return
    if figfilename != None:
        fig.savefig(figfilename,bbox_inches='tight')
    else:
        plt.show()
    return fig

def paper_stackplot(lx,plug,figfilename=None,ax=None,legendloc=None,**kwargs):
    colordict = haloutils.get_colors_for_halos()
    haloids = get_haloidlist(1)
    labeldict = haloutils.hid2name
    if ax == None: 
        fig,ax = plt.subplots(figsize=(8,8))
        plotfig = True
    else:
        assert figfilename==None,'Cannot specify both ax and figfilename'
        plotfig = False
    for hid in haloids:
        hpath = haloutils.get_hpath_lx(hid,lx)
        plug.plot(hpath,ax,color=colordict[hid],label=labeldict[hid],**kwargs)
    if legendloc != None:
        ax.legend(loc=legendloc,ncol=2,frameon=False)
    if not plotfig: return
    if figfilename != None:
        fig.savefig(figfilename,bbox_inches='tight')
    else:
        plt.show()
    return fig

def animated_stackplot(haloids,lx,plug,figfilename=None):
    # TODO
    # Make one figure with all transparent lines
    # Make N figures with one halo id highlighted using solid line and hid label
    raise NotImplementedError

def haloplot(hid,lx,pluglist,savefig=False,savepath=None,pdf=False,eps=False,normtohost=False,**kwargs):
    hpath = haloutils.get_hpath_lx(hid,lx)
    hidstr = haloutils.hidstr(hid)
    ictype,lx,nv = haloutils.get_zoom_params(hpath)
    figlist = []
    for plug in pluglist:
        fig,ax = plt.subplots()
        plug.plot(hpath,ax,normtohost=normtohost,**kwargs)
        figlist.append(fig)
        if savefig:
            assert plug.autofigname != ''
            # TODO automatic savepath inside the actual halo directory?
            figfilename = './'
            figfilename += hidstr+'_LX'+str(lx)+'_'+plug.autofigname
            if normtohost: figfilename += '_norm'
            if pdf:   figfilename += '.pdf'
            elif eps: figfilename += '.eps'
            else:     figfilename += '.png'
            fig.savefig(figfilename,bbox_inches='tight')
    if not savefig:
        plt.show()
    return figlist

def plot_5x5(plug,lx=14,figfilename=None,**kwargs):
    fig,axarr = plt.subplots(5,5,figsize=(12,12))
    for i,ax in enumerate(np.ravel(axarr)):
        hid = haloutils.cid2hid[i+1]
        hpath = haloutils.get_hpath_lx(hid,lx)
        plug.plot(hpath,ax,**kwargs)
        if hpath==None: continue
        plug.label_plot(hpath,ax,label='catnum')
    for i in range(4):
        for j in range(1,5):
            axarr[i,j].set_xlabel('')
            axarr[i,j].set_ylabel('')
    for i in range(4):
        axarr[i,0].set_xlabel('')
    for j in range(1,5):
        axarr[4,j].set_ylabel('')
    if figfilename != None:
        fig.savefig(figfilename,bbox_inches='tight')
    else:
        plt.show()
    return fig
