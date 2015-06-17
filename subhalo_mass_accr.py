import numpy as np
import pylab as plt
import os, sys, subprocess, time
import cPickle as pickle
from seaborn.apionly import rugplot

import haloutils
from caterpillaranalysis import PluginBase
from MTanalysis2 import getInfall

class SubMassAccrPlugin(PluginBase):
    def __init__(self):
        super(SubMassAccrPlugin,self).__init__()
        self.filename='submassaccr.p'
        self.xmin = 0; self.xmax = 1
        self.ymin = 10**5; self.ymax = 10**11
        self.xlog = False; self.ylog = True
        self.xlabel = 'scale'
        self.ylabel = 'mass'
    def _analyze(self,hpath):
        if not haloutils.check_mergertree_exists(hpath):
            raise IOError("No Merger Tree")
        # Load mtc of host halo and all its subs at z=0
        mtc = haloutils.load_zoom_mtc(hpath,indexbyrsid=True)
        hostid = haloutils.load_zoomid(hpath)
        #hostmb = mtc[hostid].getMainBranch()

        header = haloutils.get_halo_header(hpath)

        numsnaps = haloutils.get_numsnaps(hpath)
        lastsnap = haloutils.get_lastsnap(hpath)
        numsubs = len(mtc.Trees)-1
        massaccrarr = np.zeros((numsubs,numsnaps))
        idarr = np.zeros(numsubs)
        massarr = np.zeros(numsubs)
        vmaxarr = np.zeros(numsubs)
        mpeakarr = np.zeros(numsubs)
        mpeaksnaparr = np.zeros(numsubs)
        vpeakarr = np.zeros(numsubs)
        vpeaksnaparr = np.zeros(numsubs)
        
        foundhost=0
        for i,(rsid,mt) in enumerate(mtc.Trees.iteritems()):
            if rsid==hostid: foundhost=1; continue
            ix = i-foundhost
            
            mb = mt.getMainBranch()
            massaccrarr[ix,mb['snap']] = mb['mvir']
            idarr[ix] = rsid
            massarr[ix] = mb['mvir'][0]
            vmaxarr[ix] = mb['vmax'][0]
            
            ix_mpeak = np.argmax(mb['mvir'])
            mpeakarr[ix] = mb['mvir'][ix_mpeak]
            mpeaksnaparr[ix] = mb['snap'][ix_mpeak]
            
            ix_vpeak = np.argmax(mb['vmax'])
            vpeakarr[ix] = mb['vmax'][ix_vpeak]
            vpeaksnaparr[ix] = mb['snap'][ix_vpeak]

        # sort by mpeak
        iisort = np.argsort(mpeakarr)[::-1]
        idarr = idarr[iisort]
        massaccrarr = massaccrarr[iisort,:]/header.hubble
        vmaxarr = vmaxarr[iisort]
        mpeakarr = mpeakarr[iisort]/header.hubble
        mpeaksnaparr = mpeaksnaparr[iisort]
        vpeakarr = vpeakarr[iisort]
        vpeaksnaparr = vpeaksnaparr[iisort]
            
        with open(self.get_outfname(hpath),'w') as f:
            pickle.dump([idarr,massaccrarr,vmaxarr,mpeakarr,mpeaksnaparr,vpeakarr,vpeaksnaparr],f)
    def _read(self,hpath):
        with open(self.get_outfname(hpath),'r') as f:
            data = pickle.load(f)
        idarr,massaccrarr,vmaxarr,mpeakarr,mpeaksnaparr,vpeakarr,vpeaksnaparr = data
        idarr = idarr.astype(int)
        mpeaksnaparr = mpeaksnaparr.astype(int)
        vpeaksnaparr = vpeaksnaparr.astype(int)
        return idarr,massaccrarr,vmaxarr,mpeakarr,mpeaksnaparr,vpeakarr,vpeaksnaparr
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,alpha=0.1,
              logMpeakcut=None,whatrug='mpeak',**kwargs):
        assert whatrug in [None,'mpeak','vpeak'] #TODO 'infall'
        idarr,massaccrarr,vmaxarr,mpeakarr,mpeaksnaparr,vpeakarr,vpeaksnaparr = data

        
        numsnaps = haloutils.get_numsnaps(hpath)
        snaparr = np.arange(numsnaps)
        scalearr = haloutils.get_scale_snap(hpath,snaparr)
        assert massaccrarr.shape[1] == len(scalearr)
        
        if logMpeakcut != None:
            ii = np.log10(mpeakarr) > logMpeakcut
            idarr = idarr[ii]
            massaccrarr = massaccrarr[ii,:]
            vmaxarr = vmaxarr[ii]
            mpeakarr = mpeakarr[ii]
            mpeaksnaparr = mpeaksnaparr[ii]
            vpeakarr = vpeakarr[ii]
            vpeaksnaparr = vpeaksnaparr[ii]
            
        if 'color' in kwargs:
            color = kwargs.pop('color')
        else:
            if lx != None:
                color = self.colordict[lx]
            else:
                color = 'k'

        for massrow in massaccrarr:
            ax.plot(scalearr,massrow,color=color,alpha=alpha,**kwargs)

        if whatrug != None:
            if whatrug=='mpeak':
                rugdat = scalearr[mpeaksnaparr]
            elif whatrug=='vpeak':
                rugdat = scalearr[vpeaksnaparr]
            rugplot(rugdat,ax=ax,color=color,alpha=alpha,height=self.ymin*10**.5)
