import numpy as np
import pylab as plt
import os,subprocess,sys
import asciitable
import pynbody as pnb
import readsnapshots.readsnapHDF5_greg as rsg
from profiles.densityprofile import densityprofilesorted,getr200sorted
import haloutils

class PluginBase(object):
    """
    When extending this class, make sure to define the following variables in __init__:
    Data: filename
    Plotting: xmin, xmax, ymin, ymax, xlog, ylog, xlabel, ylabel
    Figure name (for haloplot): autofigname

    Define the following methods:
    _analyze(self,hpath)
        compute relevant quantities, save in hpath/OUTPUTFOLDERNAME
    _read(self,hpath)
        read the data computed/saved with _analyze(), return data
    _plot(self,hpath,data,ax,lx=None,**kwargs)
        take data from _read() and plot in ax.
        lx input is used when stacking multiple LX's on same plot (see convergeplot)
        **kwargs should be used for plotting functions
    """
    colordict = {11:'b',12:'r',13:'g',14:'m'}
    OUTPUTFOLDERNAME = 'analysis'
    def __init__(self):
        self.filename=None
        self.allhalos=True
        self.radius=None
        self.verbose=False
        self.xmin=None
        self.xmax=None
        self.ymin=None
        self.ymax=None
        self.xlog=None
        self.ylog=None
        self.xlabel=''
        self.ylabel=''
        self.autofigname=None
    def get_outfname(self,hpath):
        subprocess.call("mkdir -p "+hpath+'/'+self.OUTPUTFOLDERNAME,shell=True)
        return hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename
    def get_filename(self,hpath):
        assert hpath != None
        assert self.filename != None
        fname = hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename
        if not os.path.exists(fname): raise IOError
        return fname
    def file_exists(self,hpath):
        if hpath==None: return False
        try:
            fname = self.get_filename(hpath)
            return True
        except IOError:
            return False
    def _analyze(hpath):
        raise NotImplementedError
    def analyze(self,hpath,recalc=False):
        if hpath==None: return
        if recalc: self._analyze(hpath)
        else:
            if self.file_exists(hpath):
                if self.verbose: print "Already analyzed: "+haloutils.get_foldername(hpath)
            else: self._analyze(hpath)
    def get_rssubs(self,rscat,zoomid):
        if self.allhalos:
            return rscat.get_all_subhalos_within_halo(zoomid,radius=self.radius)
        else:
            return rscat.get_subhalos_within_halo(zoomid,radius=self.radius)

    def _read(self,hpath):
        raise NotImplementedError
    def read(self,hpath,autocalc=True,recalc=False):
        if hpath==None: return None
        if self.file_exists(hpath):
            return self._read(hpath)
        elif autocalc:
            try:
                print "Automatically analyzing "+haloutils.get_foldername(hpath)+"..."
                self.analyze(hpath,recalc=recalc)
                print "Done!"
                return self._read(hpath)
            except:
                print "Automatic analysis failed..."
                print sys.exc_info()
                return None
        else:
            return None

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        raise NotImplementedError
    def plot(self,hpath,ax,lx=None,labelon=False,recalc=False,**kwargs):
        data = self.read(hpath,recalc=recalc)
        if data==None:
            self.format_plot(ax)
            return
        self._plot(hpath,data,ax,lx=lx,labelon=labelon,**kwargs)
        self.format_plot(ax)
    def lxplot(self,hid,ax,**kwargs):
        lxlist = haloutils.get_lxlist(hid)
        hpaths = haloutils.get_lxlist(hid,gethpaths=True)
        for lx,hpath in zip(lxlist,hpaths):
            self.plot(hpath,ax,lx=lx,labelon=True,**kwargs)

    def customplot(self,ax,*args,**kwargs):
        raise NotImplementedError
    def format_plot(self,ax):
        assert self.xmin != None and self.xmax != None
        assert self.ymin != None and self.ymax != None
        assert self.xlog != None and self.ylog != None
        ax.set_xlim((self.xmin,self.xmax))
        ax.set_ylim((self.ymin,self.ymax))
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)
        if self.xlog: ax.set_xscale('log')
        if self.ylog: ax.set_yscale('log')

class NvmaxPlugin(PluginBase):
    def __init__(self,vmin=0.3,vmax=100,Nmin=1,Nmax=10**4.5):
        super(NvmaxPlugin,self).__init__()
        self.filename='Nvmax.dat'
        self.logvmin = -1.
        self.logvmax = 3.
        self.dlogv = 0.05
        self.vmaxbins = 10.**np.arange(self.logvmin,self.logvmax+self.dlogv,self.dlogv)
        self.xmin = vmin; self.xmax = vmax
        self.ymin = Nmin;  self.ymax = Nmax
        self.xlabel = r'$V_{\rm max}$ [km/s]'
        self.ylabel = r'N($>V_{\rm max}$)'
        self.xlog = True; self.ylog = True
        self.autofigname='Nvmax'
    def calcNvmax(self,vmax):
        h,x = np.histogram(vmax,bins=self.vmaxbins)
        return np.cumsum(h[::-1])[::-1]
    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
        numsnaps = haloutils.get_numsnaps(hpath)
        rscat = haloutils.load_rscat(hpath,numsnaps-1)
        zoomid = haloutils.load_zoomid(hpath)
        eps = 1000*haloutils.load_soft(hpath)
        
        subs = self.get_rssubs(rscat,zoomid)
        svmax = np.array(subs['vmax'])
        srmax = np.array(subs['rvmax'])
        svmaxp = svmax * np.sqrt(1+(eps/srmax)**2)

        Nvmax  = self.calcNvmax(svmax)
        Nvmaxp = self.calcNvmax(svmaxp)

        try:
            scat = haloutils.load_scat(hpath)
            bestgroup = 0
            ssvmax = scat.sub_vmax[0:scat.group_nsubs[0]]
            ssrmax = scat.sub_vmaxrad[0:scat.group_nsubs[0]]
            ssvmaxp = ssvmax*np.sqrt(1+((eps/1000.)/ssrmax)**2)
            sNvmax  = self.calcNvmax(ssvmax)
            sNvmaxp = self.calcNvmax(ssvmaxp)
        except IOError: #No Subfind
            ssvmax = 0
            sNvmax = np.zeros(len(Nvmax))
            sNvmaxp = np.zeros(len(Nvmax))

        with open(self.get_outfname(hpath),'w') as f:
            f.write(str(np.min(svmax))+" "+str(np.min(ssvmax))+'\n')
            for v,N,sN,Np,sNp in zip(self.vmaxbins[1:],Nvmax,sNvmax,Nvmaxp,sNvmaxp):
                f.write(str(v)+" "+str(N)+" "+str(sN)+" "+str(Np)+" "+str(sNp)+'\n')
    def _read(self,hpath):
        thisfilename = self.get_filename(hpath)
        data = asciitable.read(thisfilename,delimiter=' ',data_start=1)
        v  = data['col1']
        N  = data['col2']
        sN = data['col3']
        Np = data['col4']
        sNp= data['col5']
        with open(thisfilename,'r') as f:
            split = f.readline().split(" ")
            minv = float(split[0])
            sminv= float(split[1])
        return v,N,minv,sN,sminv,Np,sNp
    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        v,N,minv,sN,sminv,Np,sNp = data
        ii = np.where(v >= minv)
        if lx != None:
            ax.plot(v[ii],N[ii],color=self.colordict[lx],**kwargs)
        else:
            ax.plot(v[ii],N[ii],**kwargs)
        if labelon:
            plotlabel = haloutils.hidstr(haloutils.get_parent_hid(hpath))
            ax.text(self.xmin * 10**0.1, self.ymax * 10**-0.5, plotlabel,color='black',fontsize='medium')

class SHMFPlugin(PluginBase):
    def __init__(self,Mmin=10**4.5,Mmax=10**10.6,ymin=10**-10,ymax=10**-1.0):
        super(SHMFPlugin,self).__init__()
        self.filename='SHMF.dat'
        self.histrange = np.arange(4.0,10.5,0.2)

        self.xmin = Mmin; self.xmax = Mmax
        self.ymin = ymin;  self.ymax = ymax
        self.xlabel = r'$M_{\rm sub} [h^{-1} M_\odot]$'
        self.ylabel = r'$dn/dM_{\rm sub} [h/M_\odot]$'
        self.xlog = True; self.ylog = True
        self.autofigname = 'SHMF'
    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
        numsnaps = haloutils.get_numsnaps(hpath)
        rscat = haloutils.load_rscat(hpath,numsnaps-1)
        zoomid = haloutils.load_zoomid(hpath)
        
        subs = self.get_rssubs(rscat,zoomid)
        subM = np.array(subs['mvir'])
        x,y = self.MassFunc_dNdM(subM,self.histrange)

        try:
            scat = haloutils.load_scat(hpath)
            bestgroup = 0
            ssubM = scat.sub_mass[0:scat.group_nsubs[0]]*10**10
            sx,sy = self.MassFunc_dNdM(ssubM,self.histrange)
        except IOError: #No Subfind
            sx = np.zeros(len(x))
            sy = np.zeros(len(y))

        with open(self.get_outfname(hpath),'w') as f:
            for a,b,sa,sb in zip(x,y,sx,sy):
                f.write(str(a)+' '+str(b)+' '+str(sa)+' '+str(sb)+'\n')
    def MassFunc_dNdM(self,masses,histrange):
        """
        Adapted from Greg's MassFunctions code
        """
        numbins = len(histrange) - 1
        hist, r_array = np.histogram(np.log10(masses), bins=histrange)
        x_array = self._getMidpoints(r_array)
        dM = 10.**r_array[1:]-10.**r_array[0:numbins] #Mass size of bins in non-log space
        dNdM = hist/dM
        return 10**x_array, dNdM
    def _getMidpoints(self,bins):
        spacing = bins[1:]-bins[:-1]
        return bins[:-1]+spacing/2.0

    def _read(self,hpath):
        thisfilename = self.get_filename(hpath)
        data = asciitable.read(thisfilename,delimiter=' ')
        return data['col1'],data['col2'],data['col3'],data['col4']
    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        x,y,sx,sy = data
        if lx != None:
            ax.plot(x,y,color=self.colordict[lx],**kwargs)
        else:
            ax.plot(x,y,**kwargs)
        if labelon:
            plotlabel = haloutils.hidstr(haloutils.get_parent_hid(hpath))
            ax.text(self.xmin*10**0.2,self.ymax*10**-.75,plotlabel,color='black',fontsize='medium')

class ProfilePlugin(PluginBase):
    def __init__(self,rmin=10**-2,rmax=10**3,ymin=10**-1.5,ymax=10**2.5):
        super(ProfilePlugin,self).__init__()
        self.filename='rsprofile.dat'

        self.xmin = rmin; self.xmax = rmax
        self.ymin = ymin;  self.ymax = ymax
        self.xlabel = r'r [$h^{-1}$ kpc]'
        self.ylabel = r'$r^2 \rho(r)$ [$10^{10} M_\odot$ Mpc$^{-1}$]'
        self.xlog = True; self.ylog = True
        self.autofigname = 'rhor2'
    def _analyze(self,hpath):
        rarr, rhoarr, p03rmin, halorvir, r200c, halomass = self.compute_one_profile(hpath)
        with open(self.get_outfname(hpath),'w') as f:
            f.write(str(p03rmin)+" "+str(halorvir)+" "+str(r200c)+" "+str(halomass)+"\n")
            for r,rho in zip(rarr[2:],rhoarr[1:]):
                f.write(str(r)+" "+str(rho)+"\n")
    def get_rarr(self):
        return np.logspace(-5,0,50)
    def compute_one_profile(self,hpath,zoomid=-1,rarr=-1,snap=255):
        snapstr = str(snap).zfill(3)
        snapfile = hpath+'/outputs/snapdir_'+snapstr+'/snap_'+snapstr
        header = rsg.snapshot_header(snapfile+'.0')
        if rarr==-1: rarr = self.get_rarr()

        rscat = haloutils.load_rscat(hpath,snap)
        haloid = haloutils.get_parent_hid(hpath)
        ictype,lx,nv = haloutils.get_zoom_params(hpath)
        if snap==255 and zoomid==-1:
            zoomid = haloutils.load_zoomid(hpath)
        else: assert zoomid != -1,'Need to pass in zoomid for snap != 255'
        haloparts = rscat.get_all_particles_from_halo(zoomid)
        halopos = np.array(rscat.ix[zoomid][['posX','posY','posZ']])
        halorvir = float(rscat.ix[zoomid]['rvir']) #h^-1 kpc; rvir is 18pi^2
        halomass = rscat.ix[zoomid]['mvir']/header.hubble

        try:
            haloparts = np.sort(haloparts)
            partpos = haloutils.load_partblock(hpath,snap,"POS ",parttype=1,ids=haloparts)
            rhoarr,p03rmin = densityprofilesorted(rarr,partpos,header,haloparts,halopos,power03=True)
            r200c = getr200sorted(haloparts,partpos,header,halopos)
        except IndexError as e:
            raise RuntimeError("Contamination in halo")
        return rarr, rhoarr, p03rmin, halorvir, r200c, halomass

    def _read(self,hpath):
        thisfilename = self.get_filename(hpath)
        data = np.array(asciitable.read(thisfilename,delimiter=" ",data_start=1))
        r = 1000.*data['col1']
        rho = data['col2']
        f = open(thisfilename,'r')
        p03r,rvir,r200c,halomass = f.readline().split(" ")
        p03r = 1000.*float(p03r); rvir = float(rvir); r200c = float(r200c)
        return r,rho,p03r,rvir,r200c

    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        r,rho,p03r,rvir,r200c = data
        eps = 1000*haloutils.load_soft(hpath)
        ii1 = r >= eps
        ii2 = r >= p03r
        if lx != None:
            color = self.colordict[lx]
            ax.plot(r[ii1], (r[ii1]/1000.)**2 * rho[ii1], color=color, lw=1, **kwargs)
            ax.plot(r[ii2], (r[ii2]/1000.)**2 * rho[ii2], color=color, lw=3, **kwargs)
        else:
            ax.plot(r[ii1], (r[ii1]/1000.)**2 * rho[ii1], lw=1, **kwargs)
            ax.plot(r[ii2], (r[ii2]/1000.)**2 * rho[ii2], lw=3, **kwargs)
        if labelon:
            plotlabel = haloutils.hidstr(haloutils.get_parent_hid(hpath))
            ax.text(self.xmin*10**0.2,self.ymax*10**-.3,plotlabel,color='black',fontsize='medium')

class MassAccrPlugin(PluginBase):
    def __init__(self,Mmin=10**4.5,Mmax=10**10.6,ymin=10**-10,ymax=10**-1.0):
        super(MassAccrPlugin,self).__init__()
        self.filename='massaccr.dat'

        self.xmin = 0; self.xmax = 1
        self.ymin = 10**6;  self.ymax = 10**13
        self.xlabel = r'scale'
        self.ylabel = r'$M$ [$M_\odot$]'
        self.xlog = False; self.ylog = True
        self.autofigname = 'massaccr'
    def _analyze(self,hpath):
        if not haloutils.check_mergertree_exists(hpath,autoconvert=True):
            raise IOError("No Merger Tree")
        zoomid = haloutils.load_zoomid(hpath)
        rscat = haloutils.load_rscat(hpath,haloutils.get_numsnaps(hpath)-1)
        mtc = haloutils.load_mtc(hpath,haloids=[zoomid])
        mt = mtc[0]
        mb = mt.getMainBranch()
        scale = mb['scale'][::-1]
        snap = mb['snap'][::-1]
        mvir = mb['mvir'][::-1]/rscat.h0
        sammvir = mb['sam_mvir'][::-1]/rscat.h0
        vmax = mb['vmax'][::-1]
        TU = mb['T/|U|'][::-1]
        scaleMM = mb['scale_of_last_MM'][::-1]
        x = mb['posX'][::-1]
        y = mb['posY'][::-1]
        z = mb['posZ'][::-1]
        spin = mb['spin'][::-1]
        spinbullock = mb['spin_bullock'][::-1]
        asciitable.write({'scale':scale,'snap':snap,
                          'mvir':mvir,'sam_mvir':sammvir,
                          'vmax':vmax,'T/|U|':TU,
                          'scale_of_last_MM':scaleMM,
                          'x':x,'y':y,'z':z,
                          'spin':spin,'spin_bullock':spinbullock},
                         self.get_outfname(hpath),
                         names=['scale','snap','mvir','sam_mvir','vmax',
                                'T/|U|','scale_of_last_MM','x','y','z',
                                'spin','spin_bullock'])
    def _read(self,hpath):
        thisfilename = self.get_filename(hpath)
        tab = asciitable.read(thisfilename,header_start=0)
        return tab
    def _plot(self,hpath,data,ax,lx=None,labelon=False,**kwargs):
        tab = data
        x = tab['scale']
        y = tab['mvir']
        if lx != None:
            ax.plot(x,y,color=self.colordict[lx],**kwargs)
        else:
            ax.plot(x,y,**kwargs)
        if labelon:
            plotlabel = haloutils.hidstr(haloutils.get_parent_hid(hpath))
            ax.text(self.xmin+0.1,self.ymax*10**-.5,plotlabel,color='black',fontsize='medium')
