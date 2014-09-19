import numpy as np
import pylab as plt
import os,asciitable
import haloutils,sheetplot

from profilefit import NFWprofile,fitNFW

class ProfileReader(sheetplot.ReaderBase):
    def __init__(self,whichtype,ictype="BB",nv=4):
        self.whichtype = whichtype
        if whichtype == 0: 
            self.filename = 'rs-halo-profile.dat'
            plotextra = 'rockstarnosubs'
        elif whichtype == 1: 
            self.filename = 'rs-halo-profile-allpart.dat'
            plotextra = 'rockstarall'
        elif whichtype == 2: 
            self.filename = 'subf-halo-profile.dat'
            plotextra = 'subfind'
        elif whichtype == 3: 
            self.filename = 'subf-halo-profile-radius.dat'
            plotextra = 'subfindradius'
        else: raise ValueError("whichtype must be 0,1,2,3")
        assert ictype.upper() in ["BB","BE","BC"]
        self.ictype = ictype
        self.nv = nv

    def __call__(self,hid):
        ictype = self.ictype; nv = self.nv
        lxlist = self.get_lxlist(hid)
        rlist = []; rholist = []; p03rlist = []
        rvirlist = []; r200clist = []
        for lx in lxlist:
            thisfilename = self.get_filename(hid,lx)
            data = np.array(asciitable.read(thisfilename,delimiter=" ",data_start=1))
            rlist.append(1000.*data['col1'])
            rholist.append(data['col2'])
            f = open(thisfilename,'r')
            p03r,rvir,r200c,halomass = f.readline().split(" ")
            p03r = 1000.*float(p03r); rvir = float(rvir); r200c = float(r200c)
            f.close()
            p03rlist.append(p03r); rvirlist.append(rvir)
            r200clist.append(r200c)
        return hid,lxlist,rlist,rholist,p03rlist,rvirlist,r200clist

class ProfilePlotter(sheetplot.PlotterBase):
    def __init__(self,whichtype):
        self.whichtype = whichtype
        if whichtype == 0: 
            self.filename = 'rs-halo-profile.dat'
            plotextra = 'rockstarnosubs'
        elif whichtype == 1: 
            self.filename = 'rs-halo-profile-allpart.dat'
            plotextra = 'rockstarall'
        elif whichtype == 2: 
            self.filename = 'subf-halo-profile.dat'
            plotextra = 'subfind'
        elif whichtype == 3: 
            self.filename = 'subf-halo-profile-radius.dat'
            plotextra = 'subfindradius'
        else: raise ValueError("whichtype must be 0,1,2,3")
        self.fprefix = 'rhor2'
        self.fpostfix= plotextra
    def __call__(self,ax,data,fignum=None,**kwargs):
        hid,lxlist,rlist,rholist,p03rlist,rvirlist,r200clist = data
        ymin = 10**-1.5; ymax = 10**2.5
        ax.set_xscale('log'); ax.set_yscale('log')
        for i,lx in enumerate(lxlist):
            color = self.colordict[lx]
            r = rlist[i]; rho = rholist[i]
            p03r = p03rlist[i]; rvir = rvirlist[i]; r200c = r200clist[i]
            ax.plot(r,(r/1000.)**2 * rho, color=color, **kwargs)
            ax.plot([p03r,p03r],[ymin,ymax],color=color,ls='--',**kwargs)
            ax.plot([rvir,rvir],[ymin,ymax],'k-.')
            ax.plot([r200c,r200c],[ymin,ymax],'k:')

            good_r = (r > p03r) & (r < rvir)
            rs,rhos = fitNFW(r[good_r],rho[good_r],x0=[20,6],bounds=[(1,300),(4,8)])
            ax.plot(r,(r/1000.)**2 * NFWprofile(r,rs,rhos),ls=':',lw='0.5',**kwargs)

            yexponent = 1.3 + 0.2*(14-lx)
            ax.text(10**-1.8,10**yexponent,r"LX%i $r_s=$%3.2f kpc" % (lx,rs),
                    color=color,fontsize='x-small')

        ax.set_xlabel(r'r [$h^{-1}$ kpc]')
        ax.set_ylabel(r'$r^2 \rho(r)$ [$10^{10} M_\odot$ Mpc$^{-1}$]')
        ax.set_xlim([10**-2,10**3])
        ax.set_ylim([ymin,ymax])
        plotlabel = haloutils.hidstr(hid)
        ax.text(10**-1.8,10**2.2,plotlabel,color='black',fontsize='small')
        
