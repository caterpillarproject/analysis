import numpy as np
import pylab as plt
import warnings

from numpy import random

class AbundMatch(object):
    def __init__(self):
        raise NotImplementedError
    def get_Mstar(self,Mhalo,z=0.):
        raise NotImplementedError
    def plot(self,z=0.,ax=None,plotratio=False,**kwargs):
        Mhaloarr = np.logspace(6,15,101)
        Mstararr = self.get_Mstar(Mhaloarr,z=z)
        ratio = Mstararr/Mhaloarr
        if ax==None:
            fig,ax = plt.subplots(figsize=(8,8))
            ax.set_xscale('log')
            ax.set_yscale('log')
        if plotratio:
            ax.plot(Mhaloarr,ratio,**kwargs)
        else:
            ax.plot(Mhaloarr,Mstararr,**kwargs)
        return ax
    def __getitem__(self,key):
        raise NotImplementedError

class Behroozi13AbundMatch(AbundMatch):
    def __init__(self):
        self.M10 = 11.514
        self.M1a = -1.793
        self.M1z = -0.251
        self.eps0 = -1.777
        self.epsa = -0.006
        self.epsz = -0.000
        self.epsa2= -0.119
        self.alpha0 = -1.412
        self.alphaa = 0.731
        self.delta0 = 3.508
        self.deltaa = 2.608
        self.deltaz = -0.251
        self.gamma0 = 0.316
        self.gammaa = 1.319
        self.gammaz = 0.279
        self.xi0 = 0.218
        self.xia = -0.023
        self.MICL0 = 12.515
        self.MICLa = -2.503
        
    def __getitem__(self,key):
        raise NotImplementedError
    def get_logMstar(self,Mhalo,z):
        logM1,logeps,alpha,delta,gamma,xi,logMICL = self.scale_params(1./(1.+z))
        def _f(x):
            return -np.log10(10**(alpha*x)+1) + delta*(np.log10(1+np.exp(x)))**gamma / (1+np.exp(10**-x))
        return logM1 + logeps + _f(np.log10(Mhalo)-logM1)-_f(0)
    def get_Mstar(self,Mhalo,z=0.):
        return 10.**self.get_logMstar(Mhalo,z)
    def get_Mstar_scatter(self,Mhalo,z=0.):
        xi = self._xi(1./(1.+z))
        try:
            N = len(Mhalo)
        except:
            N = 1
        scatter = random.normal(scale=xi,size=N)
        return self.get_Mstar(Mhalo,z=z)*10**scatter

    def scale_params(self,a):
        nu = np.exp(-4.*a*a)
        z = 1./a - 1
        logM1 = self.M10 + nu*(self.M1a*(a-1.)+self.M1z*z)
        logeps = self.eps0 + nu*(self.epsa*(a-1.)+self.epsz*z) + self.epsa2*(a-1.)
        alpha = self.alpha0 + nu*(self.alphaa*(a-1.))
        delta = self.delta0 + nu*(self.deltaa*(a-1.)+self.deltaz*z)
        gamma = self.gamma0 + nu*(self.gammaa*(a-1.)+self.gammaz*z)
        xi = self.xi0 + self.xia*(a-1.)
        logMICL = self.MICL0 + self.MICLa*(a-1.)
        return logM1,logeps,alpha,delta,gamma,xi,logMICL
    
    def _nu(self,a):
        return np.exp(-4.*a*a)
    def _logM1(self,a):
        return self.M10 + self._nu(a)*(self.M1a*(a-1.)+self.M1z*(1./a-1.))
    def _logeps(self,a):
        return self.eps0 + self._nu(a)*(self.epsa*(a-1.)+self.epsz*(1./a-1.)) + self.epsa2*(a-1.)
    def _alpha(self,a):
        return self.alpha0 + self._nu(a)*(self.alphaa*(a-1.))
    def _delta(self,a):
        return self.delta0 + self._nu(a)*(self.deltaa*(a-1.)+self.deltaz*(1./a-1.))
    def _gamma(self,a):
        return self.gamma0 + self._nu(a)*(self.gammaa*(a-1.)+self.gammaz*(1./a-1.))
    def _xi(self,a):
        return self.xi0 + self.xia*(a-1.)

class GK14AbundMatch(Behroozi13AbundMatch):
    def __init__(self):
        super(GK14AbundMatch,self).__init__()
        self.alpha0 = -1.92
    def get_logMstar(self,Mhalo,z):
        if z != 0: warnings.warn("GK14 has only calibrated their AM with z=0")
        return super(GK14AbundMatch,self).get_logMstar(Mhalo,z)

class Moster13AbundMatch(AbundMatch):
    def __init__(self):
        self.M10 = 11.590
        self.M11 = 1.195
        self.N10 = 0.0351
        self.N11 = -0.0247
        self.beta10 = 1.376
        self.beta11 = 0.153
        self.gamma10 = 0.608
        self.gamma11 = 0.173

        self.scatter = 0.15
    def __getitem__(self,key):
        raise NotImplementedError
    def get_Mratio(self,Mhalo,z):
        logM1,N,beta,gamma = self.scale_params(1./(1.+z))
        Mrat = Mhalo/10**logM1
        return 2*N/(Mrat**(-beta) + Mrat**(gamma))
    def get_logMstar(self,Mhalo,z):
        return np.log10(self.get_Mstar(Mhalo,z))
    def get_Mstar(self,Mhalo,z=0.):
        return self.get_Mratio(Mhalo,z)*Mhalo
    def scale_params(self,a):
        logM1 = self.M10 + self.M11*(1-a)
        N = self.N10 + self.N11*(1-a)
        beta = self.beta10 + self.beta11*(1-a)
        gamma = self.gamma10 + self.gamma11*(1-a)
        return logM1,N,beta,gamma

    def get_Mstar_scatter(self,Mhalo,z=0.):
        try:
            N = len(Mhalo)
        except:
            N = 1
        scatter = random.normal(scale=self.scatter,size=N)
        return self.get_Mstar(Mhalo,z=z)*10**scatter
