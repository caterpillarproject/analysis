import numpy as np
from scipy.optimize import fmin_l_bfgs_b

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
    #return np.log10(EINprofile(r,r2,10**logrho2,alpha))

def fitNFW(rarr,rhoarr,x0=[.05,6.5],bounds=[(.001,3),(5,8)]):
    nbins = len(rarr)
    logrho = np.log10(rhoarr)
    def Q2(x):
        rs,logrhos = x
        logrhomodel = logNFWprofile(rarr,rs,logrhos)
        return np.sum((logrho-logrhomodel)**2)/nbins
    x,f,d = fmin_l_bfgs_b(Q2,x0,approx_grad=True,bounds=bounds)
    print "NFW Fit value:",x,Q2(x)
    return x[0],10**x[1]
    
def fitM99(rarr,rhoarr,x0=[.05,6.5],bounds=[(.001,3),(5,8)]):
    nbins = len(rarr)
    logrho = np.log10(rhoarr)
    def Q2(x):
        rs,logrhos = x
        logrhomodel = logM99profile(rarr,rs,logrhos)
        return np.sum((logrho-logrhomodel)**2)/nbins
    x,f,d = fmin_l_bfgs_b(Q2,x0,approx_grad=True,bounds=bounds)
    print "M99 Fit value:",x,Q2(x)
    return x[0],10**x[1]
    
def fitEinasto(rarr,rhoarr,x0=[.05,6.5,.17],bounds=[(.001,3),(5,8),(.01,.3)]):
    nbins = len(rarr)
    logrho = np.log10(rhoarr)
    def Q2(x):
        r2,logrho2,alpha = x
        logrhomodel = logEINprofile(rarr,r2,logrho2,alpha)
        return np.sum((logrho-logrhomodel)**2)/nbins
    x,f,d = fmin_l_bfgs_b(Q2,x0,approx_grad=True,bounds=bounds)
    print "EIN Fit value:",x,Q2(x)
    return x[0],10**x[1],x[2]
