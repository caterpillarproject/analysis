import matplotlib.pyplot as plt
import numpy as np
import methods
import readhalos.RSDataReaderv2 as rsr2
import AnalyticMF as MF
from fitting import *

def MassFunc_dNdM_Volume(masses, numbins, boxsize):
    """
    Produce mass function data for N(m). Fractional number density of halos in
    mass range (m,m+dm). Integrates to total number of halos/volume.
    @param masses: list of all masses in Msun/h
    @param numbins: number of bins to use in histogram
    @param boxsize: Size of box in MPc/h.
    @return: x-axis in mass, y-axis in N(m) 
    xlabel: '$M [M_\odot /h)]$'
    ylabel: '$N(M) [Num Halos/Mpc^3]$'
    """
    logmasses = np.log10(masses)
    hist, r_array = np.histogram(logmasses, numbins)
    dlogM = r_array[1]-r_array[0]
    x_array = r_array[1:] - .5*dlogM
    dM = 10.**r_array[1:]-10.**r_array[0:numbins] #Mass size of bins in non-log space.
    volume = np.float(boxsize**3) # in MPc^3
    return 10**x_array, hist/volume/dM

def MassFunc_dNdlog10M(masses, numbins, boxsize, histrange=False):
    """
    Produce mass function data for dn(m)/dlog10(M). n(m) is number of
    halos with mass < m, divided by box volume.
    mass range (m,m+dm). Integrates to total number of halos.
    @param masses: list of all masses in Msun/h
    @param numbins: number of bins to use in histogram
    @param boxsize: Size of box in MPc
    @return: x-axis in mass [Msun/h], y-axis in dn/dlog10(M) 
    xlabel: '$M [M_\odot /h)]$'
    ylabel: '$dn/d\log_{10}(M) [num halos/Mpc^3]$'
    """
    logmasses = np.log10(masses)
    ### Make n(M) number of halos with mass < M.
    if type(histrange) == bool:
        hist, r_array = np.histogram(logmasses, numbins)
    else:
        hist, r_array = np.histogram(logmasses, bins = histrange)

    dlogM = r_array[1]-r_array[0]
    x_array = r_array[1:] - .5*dlogM
    dM = 10.**r_array[1:]-10.**r_array[0:numbins] #Mass size of bins in non-log space
    dn_dlogM = np.log(10)*hist*10**x_array/dM
    volume = np.float(boxsize**3) # in MPc^3
    return 10**x_array, dn_dlogM/volume

def MassFunc_dNdM(masses, histrange):
    """
    Produce mass function data for dN/dM. dN is number of halos in mass
    interval dM. Not normalized to the volume. Best for SHMF
 
    @param masses: list of all masses in Msun
    @param histrange: binning range in logspace. ex: 
    numbins = 30
    histrange = np.arange(8.5,11.51,3./numbins)
    @return: [x-axis in log10(mass), y-axis in dN/dM, xlabel, ylabel] 
    xlabel: '$\log_{10}(M [M_\odot)])$'
    ylabel: '$\log_{10}(dN/dM) [M_\odot^{-1}]$'
    """
    numbins = len(histrange) - 1
    hist, r_array = np.histogram(np.log10(masses), bins=histrange)
    x_array = methods.getMidpoints(r_array)
    dM = 10.**r_array[1:]-10.**r_array[0:numbins] #Mass size of bins in non-log space
    print hist, 'num subs per interval'
    dNdM = hist/dM
    return 10**x_array, dNdM

def PlotMassFunc_Host(cat):
    """
    Given RSHalo catalogue of simulation, plot host halo mass function
    """
    masses_host = cat.get_hosts()['mvir']
    x_array, y_array = MassFunc_dNdlog10M(masses_host, 30, cat.box_size)
    #m_array = np.linspace(x_array[0], x_array[-1],numpoints)
    #mass_array = 10**m_array
    Predicted = MF.AnalyticMF(h=cat.h0, Omega_m=cat.Om,sigma8=0.83,n=0.96)
    STdndlogM = Predicted.STdNdLogM(x_array, z=1./cat.scale-1)

    ax = plt.axes()
    ax.plot(x_array, y_array, linewidth=2)
    ax.plot(x_array, STdndlogM, linestyle=':',linewidth=2, label='S&T')

    plt.legend()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('M $\mathrm{[M_{\odot}h^{-1}]}$')
    ax.set_ylabel('$dn/dlog_{10}(\mathrm{M})$ $[\mathrm{N_{halos}/Mpc^3 \ h^{-3}}]$')
    plt.savefig('ParentHostMassFunction')
    plt.show()

def Fit(x_axis,y_axis,error):
    """
    x_axis, y_axis in normal units.
    """
    def func(x, a, b):
        return a*x + b
    avar=-0.1
    bvar=2.5
    p0 = np.array([avar , bvar])

    mask = np.isnan(error)
    error[mask] = 0.0
    
    if error.min() == 0.0:
        popt, punc, rchi2, dof = general_fit(func, np.log10(x_axis), np.log10(y_axis), p0)
    else:
        popt, punc, rchi2, dof = general_fit(func, np.log10(x_axis), np.log10(y_axis), p0,np.log10(error))
    #print " a:", '{:.2f}'.format(popt[0]),"+-",'{:.3f}'.format(punc[0]),", b:",'{:.2f}'.format(popt[1]),"+-",'{:.3f}'.format(punc[1])
    print popt[0], 'slope', popt[1], 'intercept'
    return 10.**func(np.log10(x_axis),popt[0],popt[1])


def PlotMassFunc_SHMF(cat, haloID):
    """
    Given RSHalo Catalogue and halo ID, plot its SHMF
    """
    numbins = 12#18
    histrange = np.arange(7.0,10.6,3.6/numbins)

    subs = cat.get_subhalos_from_halo(haloID)
    masses_sub = np.array(subs['mvir'])
    x_axis, y_axis =  MassFunc_dNdM(masses_sub, histrange)
    yfit = Fit(x_axis,y_axis,np.zeros(len(x_axis)))
    ax = plt.axes()
    ax.plot(x_axis, y_axis, label='SHMF')
    ax.plot(x_axis, yfit, label='best fit')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.legend()
    plt.savefig('SHMF_%s' %(haloID))
    plt.show()
    

# parent simulation catalogue
#cat = rsr2.RSDataReader("/bigbang/data/AnnaGroup/caterpillar/parent/RockstarData/", 63, digits=2, AllParticles = False) 
#id = 103791





def vmaxfunction3(vmax,histrange):
    """
    Produce mass function data for dN/dV. dN is number of halos in mass
    interval dVmax.
    @param vmax: list of all vmaxes in km/sec
    @param numbins: number of bins to use in histogram
    @param histrange: binning range in logspace
    @return: [x-axis in log10(vmax), y-axis in dN/dVmax] 
    """
    hist, r_array = np.histogram(np.log10(vmax), bins=histrange)
    x_array = methods.getMidpoints(r_array)
    dV = 10.**r_array[1:]-10.**r_array[0:len(histrange)-1] #Mass size of bins in non-log space.
    dNdV = hist/dV
    return [x_array, dNdV]

def vmaxfunction3(nu,histrange):
    hist, r_array = np.histogram(np.log10(nu), bins=histrange)
    x_array = methods.getMidpoints(r_array)
    return [x_array, np.cumsum(hist[::-1])[::-1]]
