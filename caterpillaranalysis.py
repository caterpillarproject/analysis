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
    When extending this class, make sure to define the following variables:
    Data: filename
    Plotting: xmin, xmax, ymin, ymax, xlog, ylog, xlabel, ylabel

    Define the following methods:
    _analyze()
    _read()
    _plot()
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
    def get_filename(self,hpath):
        assert self.filename != None
        fname = hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename
        if not os.path.exists(fname): raise IOError
        return fname
    def file_exists(self,hpath):
        try:
            fname = self.get_filename(hpath)
            return True
        except IOError:
            return False
    def get_outfname(self,hpath):
        subprocess.call("mkdir -p "+hpath+'/'+self.OUTPUTFOLDERNAME,shell=True)
        return hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename
    def _analyze(hpath):
        raise NotImplementedError
    def analyze(self,hpath,recalc=False):
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

    def _plot(self,hpath,data,ax,lx=None,**kwargs):
        raise NotImplementedError
    def plot(self,hpath,ax,lx=None,**kwargs):
        data = self.read(hpath)
        if data==None:
            self.format_plot(ax)
            return
        self._plot(hpath,data,ax,lx=lx,**kwargs)
        self.format_plot(ax)
    def lxplot(self,hid,ax,**kwargs):
        lxlist = haloutils.get_lxlist(hid)
        hpaths = haloutils.get_lxlist(hid,gethpaths=True)
        for lx,hpath in zip(lxlist,hpaths):
            self.plot(hpath,ax,lx=lx,**kwargs)

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
    def _plot(self,hpath,data,ax,lx=None,**kwargs):
        v,N,minv,sN,sminv,Np,sNp = data
        ii = np.where(v >= minv)
        if lx != None:
            ax.plot(v[ii],N[ii],color=self.colordict[lx],**kwargs)
        else:
            ax.plot(v[ii],N[ii],**kwargs)
        plotlabel = haloutils.hidstr(haloutils.get_parent_hid(hpath))
        ax.text(self.xmin * 10**0.1, self.ymax * 10**-0.5, plotlabel,color='black',fontsize='medium')
