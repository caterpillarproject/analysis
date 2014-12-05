import haloutils
from caterpillaranalysis import *

def get_haloidlist(sheet):
    if sheet==1:
        haloidlist = [95289,  1195448, 1725139,
                      447649, 5320,    581141,
                      94687,  1130025, 1387186,
                      581180, 649861,  264569]
    #elif sheet==2:
    #    haloidlist = [649861,1599902, 1354437,
    #                  1631506, 1232164, 264569, 
    #                  581180, 1422331, 1725372] #1725272, 
    else:
        exit("Invalid sheet number")
    assert len(haloidlist) == 12
    return haloidlist

def convergeplot(sheetnum,plug,figfilename=None,figsize=None,**kwargs):
    haloidlist = get_haloidlist(sheetnum)
    if figsize==None: figsize=(9,11)
    fig,axarr = plt.subplots(4,3,figsize=figsize,
                             sharex=True,sharey=True)
    plt.subplots_adjust(wspace=0,hspace=0)
    axlist = axarr.reshape(-1)
    for i,(hid,ax) in enumerate(zip(haloidlist,axlist)):
        plug.lxplot(hid,ax,**kwargs)
    for ax in np.ravel(axarr[0:2,:]):
        ax.set_xlabel('')
    for ax in np.ravel(axarr[:,1:3]):
        ax.set_ylabel('')
    if figfilename != None:
        fig.savefig(figfilename,bbox_inches='tight')
    else:
        plt.show()

if __name__=="__main__":
    Nvmax = NvmaxPlugin()
    convergeplot(1,Nvmax,lw=2)
