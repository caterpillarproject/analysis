import haloutils
import numpy as np
import pylab as plt

def get_haloidlist(sheet):
    if sheet==1:
        haloidlist = [95289,  1195448, 1725139,
                      447649, 5320,    581141,
                      94687,  1130025, 1387186,
                      581180, 1725372,  264569]
    #elif sheet==2:
    #    haloidlist = [1599902, 1354437, 1631506, 
    #                  1232164, 1422331, 768257, 
    #                  649861, 1725272, 196589,
    #                  1268839. 1268839, 1268839]
    else:
        exit("Invalid sheet number")
    assert len(haloidlist) == 12
    return haloidlist

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

def stackplot(haloids,lx,plug,figfilename=None,ax=None,**kwargs):
    if ax == None: fig,ax = plt.subplots()
    else: assert figfilename==None,'Cannot specify both ax and figfilename'
    for hid in haloids:
        hpath = haloutils.get_hpath_lx(hid,lx)
        plug.plot(hpath,ax,**kwargs)
    if ax == None:
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

def haloplot(hid,lx,pluglist,savefig=False,savepath=None,pdf=False,**kwargs):
    hpath = haloutils.get_hpath_lx(hid,lx)
    hidstr = haloutils.hidstr(hid)
    ictype,lx,nv = haloutils.get_zoom_params(hpath)
    figlist = []
    for plug in pluglist:
        fig,ax = plt.subplots()
        plug.plot(hpath,ax,**kwargs)
        figlist.append(fig)
        if savefig:
            assert plug.autofigname != ''
            # TODO automatic savepath inside the actual halo directory?
            figfilename = './'
            figfilename += hidstr+'_LX'+str(lx)+'_'+plug.autofigname
            if pdf: figfilename += '.pdf'
            else:   figfilename += '.png'
            fig.savefig(figfilename,bbox_inches='tight')
    if not savefig:
        plt.show()
    return figlist

