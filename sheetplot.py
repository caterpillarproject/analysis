import pylab as plt
import numpy as np
import os,asciitable
import haloutils

def get_haloidlist(sheet):
    if sheet==1:
        haloidlist = [1194083, 1292049, 1327707,
                      231858,  1725139, 230667,
                      1269360, 649524, 706754]
    elif sheet==2:
        haloidlist = [299792, 1232333, 1476079,
                      1233057, 1327707, 1327707,
                      1327707, 1327707, 1327707]
    elif sheet==3:
        haloidlist = [1327707, 1327707, 1327707,
                      1327707, 1327707, 1327707,
                      1327707, 1327707, 1327707]
    elif sheet==10:
        haloidlist = [1327707, 1327707, 1327707,
                      1327707, 1327707, 1327707,
                      1327707, 1327707, 1327707]
    else:
        exit("Invalid sheet number")
    assert len(haloidlist) == 9
    return haloidlist

def sheetplot(sheetnum,Reader,Plotter,figsize=None,save=True,pdf=False):
    #TODO take a list of Plotter objects for the same reader
    figfilename = Plotter.fprefix+'_s'+str(sheetnum)+'_'+Plotter.fpostfix
    if pdf: figfilename += '.pdf'
    else:   figfilename += '.png'
    haloidlist = get_haloidlist(sheetnum)
    if figsize==None: figsize=(15,15)
    fig,axarr = plt.subplots(3,3,figsize=figsize,
                             sharex=True,sharey=True)
    plt.subplots_adjust(wspace=0,hspace=0)
    axlist = axarr.reshape(-1)
    for i,(hid,ax) in enumerate(zip(haloidlist,axlist)):
        data = Reader(hid)
        Plotter(ax,data,fignum=i)
    for ax in np.ravel(axarr[0:2,:]):
        ax.set_xlabel('')
    for ax in np.ravel(axarr[:,1:3]):
        ax.set_ylabel('')
    
    if save:
        fig.savefig(figfilename,bbox_inches='tight')
    else:
        plt.show()

class ReaderBase(object):
    def get_lxlist(self,hid):
        outlist = []
        for lx in [14,13,12,11]:
            if os.path.exists(self.get_filename(hid,lx)): outlist.append(lx)
        return outlist
    def get_filename(self,hid,lx):
        return haloutils.get_hpath(hid,self.ictype,lx,self.nv)+'/'+self.filename

class PlotterBase(object):
    colordict = {11:'blue',12:'red',13:'green',14:'purple'}
