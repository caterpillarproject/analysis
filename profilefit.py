import numpy as np
from scipy import integrate
from scipy.optimize import fmin_l_bfgs_b,curve_fit

def NFWprofile(r,rs,rhos):
    x = r/rs
    return rhos/(x*(1+x)**2)
def M99profile(r,rM,rhoM):
    x = r/rM
    return rhoM/(x**1.5 *(1+x)**1.5)
def EINprofile(r,r2,rho2,alpha):
    x = r/r2
    return rho2 * np.exp(-2/alpha * (x**alpha - 1))

def logNFWprofile(r,rs,logrhos):
    x = r/rs
    return logrhos-np.log10(x*(1+x)**2)
def logM99profile(r,rM,logrhoM):
    x = r/rM
    return logrhoM-np.log10(x**1.5 *(1+x)**1.5)
def logEINprofile(r,r2,logrho2,alpha):
    x = r/r2
    return logrho2+np.log10(np.e)*(-2/alpha * (x**alpha - 1))

def _ufunclike(f,x):
    return np.array(map(f,np.ravel(x)))
def NFWmltr(r,rs,rhos):
    def f(r):
        return 4*np.pi*integrate.quad(lambda x: x*x*NFWprofile(x,rs,rhos),0,r)[0]
    return _ufunclike(f,r)
def NFWmltr_analytic(r,rs,rhos):
    def F(t):
        return np.log(1+t)-t/(1.+t)
    def f(r):
        return 4*np.pi*rhos*rs**3*F(r/rs)
    return _ufunclike(f,r)
def M99mltr(r,rM,rhoM):
    return NotImplementedError
def EINmltr(r,r2,rho2,alpha):
    def f(r):
        return 4*np.pi*integrate.quad(lambda x: x*x*EINprofile(x,r2,rho2,alpha),0,r)[0]
    return _ufunclike(f,r)

def calc_rhoarr(rbin,dr,mpart):
    Varr = 4*np.pi/3 * (rbin[1:]**3-rbin[:-1]**3)
    h,x = np.histogram(dr,bins=rbin)
    Marr = h*mpart
    return Marr/Varr
def _Q2(y1,y2):
    assert len(y1)==len(y2)
    return np.sum((y1-y2)**2)/len(y1)
def fitNFW(rarr,rhoarr,p0,verbose=False,retQ2=False):
    rmid = 10**((np.log10(rarr[1:])+np.log10(rarr[:-1]))/2.)
    assert len(rmid) == len(rhoarr)
    logrho = np.log10(rhoarr)
    pNFW = curve_fit(logNFWprofile,rmid,logrho,p0=p0)[0]
    Q2 = _Q2(logrho,logNFWprofile(rmid,pNFW[0],pNFW[1]))
    if verbose: print "NFW Fit value:",pNFW,Q2
    if retQ2: return pNFW[0],10**pNFW[1],Q2
    return pNFW[0],10**pNFW[1]
def fitM99(rarr,rhoarr,p0,verbose=False,retQ2=False):
    rmid = 10**((np.log10(rarr[1:])+np.log10(rarr[:-1]))/2.)
    assert len(rmid) == len(rhoarr)
    logrho = np.log10(rhoarr)
    pM99 = curve_fit(logM99profile,rmid,logrho,p0=p0)[0]
    Q2 = _Q2(logrho,logM99profile(rmid,pM99[0],pM99[1]))
    if verbose: print "M99 Fit value:",pM99,Q2
    if retQ2: return pM99[0],10**pM99[1],Q2
    return pM99[0],10**pM99[1]
def fitEIN(rarr,rhoarr,p0,verbose=False,retQ2=False):
    rmid = 10**((np.log10(rarr[1:])+np.log10(rarr[:-1]))/2.)
    assert len(rmid) == len(rhoarr)
    logrho = np.log10(rhoarr)
    pEIN = curve_fit(logEINprofile,rmid,logrho,p0=p0,maxfev=1000000)[0]
    Q2 = _Q2(logrho,logEINprofile(rmid,pEIN[0],pEIN[1],pEIN[2]))
    if verbose: print "EIN Fit value:",pEIN,Q2
    if retQ2: return pEIN[0],10**pEIN[1],pEIN[2],Q2
    return pEIN[0],10**pEIN[1],pEIN[2]

#def fitNFW(rarr,rhoarr,x0=[.05,6.5],bounds=[(.001,3),(5,8)],verbose=False):
#    nbins = len(rarr)
#    logrho = np.log10(rhoarr)
#    def Q2(x):
#        rs,logrhos = x
#        logrhomodel = logNFWprofile(rarr,rs,logrhos)
#        return np.sum((logrho-logrhomodel)**2)/nbins
#    x,f,d = fmin_l_bfgs_b(Q2,x0,approx_grad=True,bounds=bounds)
#    if verbose: print "NFW Fit value:",x,Q2(x)
#    return x[0],10**x[1]
    
#def fitM99(rarr,rhoarr,x0=[.05,6.5],bounds=[(.001,3),(5,8)],verbose=False):
#    nbins = len(rarr)
#    logrho = np.log10(rhoarr)
#    def Q2(x):
#        rs,logrhos = x
#        logrhomodel = logM99profile(rarr,rs,logrhos)
#        return np.sum((logrho-logrhomodel)**2)/nbins
#    x,f,d = fmin_l_bfgs_b(Q2,x0,approx_grad=True,bounds=bounds)
#    if verbose: print "M99 Fit value:",x,Q2(x)
#    return x[0],10**x[1]
    
#def fitEinasto(rarr,rhoarr,x0=[.05,6.5,.17],bounds=[(.001,3),(5,8),(.01,.3)],verbose=False):
#    nbins = len(rarr)
#    logrho = np.log10(rhoarr)
#    def Q2(x):
#        r2,logrho2,alpha = x
#        logrhomodel = logEINprofile(rarr,r2,logrho2,alpha)
#        return np.sum((logrho-logrhomodel)**2)/nbins
#    x,f,d = fmin_l_bfgs_b(Q2,x0,approx_grad=True,bounds=bounds)
#    if verbose: print "EIN Fit value:",x,Q2(x)
#    return x[0],10**x[1],x[2]
