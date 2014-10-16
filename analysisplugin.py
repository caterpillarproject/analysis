import numpy as np
import haloutils
import os,subprocess

OUTPUTFOLDERNAME = 'analysis'

class AnalysisPluginBase(object):
    def __init__(self):
        raise NotImplementedError
    def __call__(self,hpath,**kwargs):
        raise NotImplementedError
    def get_outfname(self,hpath):
        subprocess.call("mkdir -p "+hpath+'/'+OUTPUTFOLDERNAME,shell=True)
        return hpath+'/'+OUTPUTFOLDERNAME+'/'+self.filename
    def get_rssubs(self,rscat,zoomid,allhalos,radius):
        if allhalos:
            return rscat.get_all_subhalos_within_halo(zoomid,radius=radius)
        else:
            return rscat.get_subhalos_within_halo(zoomid,radius=radius)

class NvmaxPlugin(AnalysisPluginBase):
    def __init__(self,logvmin=-1.,logvmax=3.,dlogv=0.05,
                 filename='Nvmax.dat',allhalos=True,radius=None,
                 verbose=False):
        self.vmaxbins = 10.**np.arange(logvmin,logvmax+dlogv,dlogv)
        self.filename = filename
        self.allhalos = allhalos
        self.radius = radius
        self.verbose = verbose
    def __call__(self,hpath):
        numsnaps = haloutils.get_numsnaps(hpath)
        rscat = haloutils.load_rscat(hpath,numsnaps-1)
        zoomid = haloutils.load_zoomid(hpath)
        eps = 1000*haloutils.load_soft(hpath)
        
        subs = self.get_rssubs(rscat,zoomid,self.allhalos,self.radius)
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
    def calcNvmax(self,vmax):
        h,x = np.histogram(vmax,bins=self.vmaxbins)
        return np.cumsum(h[::-1])[::-1]

class SHMFPlugin(AnalysisPluginBase):
    def __init__(self,logMmin=4.0,logMmax=10.5,dlogM=.2,
                 filename='SHMF.dat',allhalos=True,radius=None,
                 verbose=False):
        self.histrange = np.arange(logMmin,logMmax,0.2)
        self.filename = filename
        self.allhalos = allhalos
        self.radius = radius
        self.verbose = verbose
    def __call__(self,hpath):
        numsnaps = haloutils.get_numsnaps(hpath)
        rscat = haloutils.load_rscat(hpath,numsnaps-1)
        zoomid = haloutils.load_zoomid(hpath)
        
        subs = self.get_rssubs(rscat,zoomid,self.allhalos,self.radius)
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
