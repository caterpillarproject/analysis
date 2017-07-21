import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import haloutils as htils
from scipy import stats
import methods
import MTanalysis3 as mta
import DwarfMethods as dm
import MTanalysis_field as mtaf

import statsmodels.api as sm
from PlotParams import *


#data = np.random.rand(25)
#x = np.arange(0,1,.1)
#cumhist = [np.sum(data>x[i]) for i in range(len(x))]
#c2 = np.cumsum(np.histogram(data,bins=x)[0])[::-1]



def getMidpoints(bins):
    """                                                                                                             
    Given a range of cut-off values, find the midpoints.                                                            
    @return: array of length len(bins)-1                                                                            
    """
    spacing = bins[1:]-bins[:-1]
    return bins[:-1]+spacing/2.0


def shvmax_function(cat, haloid, target_halo_vmax=151,factor=1):
    """                           
    @ return: x-axis array of M_sub. y-axis array of dN/dM_sub   
    """
    bins_per_dex = 12
    min_vmax=np.log10(10.); max_vmax=np.log10(cat.ix[haloid]['vmax']*.8) 
    histrange = 10**np.arange(min_vmax,max_vmax+.01,1./bins_per_dex)
    radius=factor*cat.ix[haloid]['rvir']
    vmaxes_sub = np.array(cat.get_subhalos_within_halo(haloid, radius)['vmax'])
    cumhist = np.array([np.sum(vmaxes_sub>histrange[i]) for i in range(len(histrange))])
    rescaled = cumhist*(target_halo_vmax)/(cat.ix[haloid]['vmax'])
    return histrange, rescaled


def peak_vmax_function(hpath,cat,haloid,target_halo_vmax=151, version='peak',lx=14):
    if lx==13:
        AE = mta.AllExtantData()
        dataE = AE.read(hpath)
    else:
        dataE = dm.get_extant_data(hpath,False)
    peak_vmaxes = dataE['peak_vmax']/cat.h0
    bins_per_dex = 12
    min_vmax=np.log10(16); max_vmax=np.log10(cat.ix[haloid]['vmax']*0.8 )
    if lx==13:
        min_mass=np.log10(10)
    histrange = 10**np.arange(min_vmax,max_vmax+.01,1./bins_per_dex)
    cumhist = np.array([np.sum(peak_vmaxes>histrange[i]) for i in range(len(histrange))])
    print (target_halo_vmax)/(cat.ix[haloid]['vmax']), 'ratio if vmax of host'
    rescaled = cumhist*(target_halo_vmax)/(cat.ix[haloid]['vmax'])
    return histrange, rescaled


# use this function to plot SHMF for host halo, peak mass function for host halo
def plotFits(version='peak',lx=14,target_halo_vmax=151,factor=1):
    hpaths = dm.get_hpaths(field=False, lx=lx)[0:15]  #htils.get_all_halo_paths_lx(lx)
    plt.figure()
    ax = plt.subplot(111)
    y_values_full = []
    for hpath, color in zip(hpaths[0:],colorlist[0:len(hpaths)]):
        snap_z0 = htils.get_numsnaps(hpath)-1
        cat=htils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID = htils.load_zoomid(hpath)
        if version=='normal':
            x_axis, y_axis = shvmax_function(cat,hostID, target_halo_vmax, factor)
        if version[0:4]=='peak':
            x_axis, y_axis = peak_vmax_function(hpath,cat,hostID,target_halo_vmax,version,lx)
        y_values_full.append(y_axis)
        print 'done with cat', htils.hpath_catnum(hpath)
        ## Plot all the data                                                                 
        ax.plot(x_axis, y_axis, linewidth=1.0,color=color)


    ###
    # I want to average the values of each column. if NaN, ignore it. only average over what you can.
    # To keep the number of bins per dex the same, I have to use variable length y_axis.
    ###
    y_values_full = np.array(y_values_full)
    lengths = np.array([len(row) for row in y_values_full])
    maxlength = np.max(lengths)
    mean_y_axis = np.zeros(maxlength)
    for i in range(maxlength):
        print i
        mask = lengths > i
        mean_y_axis[i] = np.mean([a[i] for a in y_values_full[mask]])
        #std_err_yaxis = np.std(y_values_full,axis=0)
    lowbin = 0
    highbin = len(mean_y_axis)-2
    iszero = np.where(mean_y_axis<=0)[0]
    if len(iszero)>0:
        highbin = min(highbin, iszero[0])

    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x_axis)[lowbin:highbin],np.log10(mean_y_axis)[lowbin:highbin])
    # FIXED INTERCEPT FIT
    slope_fixed = -3.30
    fixintercept = np.mean(np.log10(mean_y_axis)[lowbin:highbin] - slope_fixed*np.log10(x_axis)[lowbin:highbin])
    
    print 'alpha = slope = ', slope
    K = 10**intercept/(target_halo_vmax)
    Kfix = 10**fixintercept/(target_halo_vmax)
    print 'K = ', K, 'Kfixed', Kfix
    print std_err, 'std err'
    plt.figtext(.15,.22, '$dn/dV_{sub} = %f \, V_{sub}^{%f}\,  V_{host}$' %(Kfix, slope_fixed))
    plt.figtext(.15,.15, '$dn/dV_{sub} = %f \, V_{sub}^{%f}\,  V_{host}$' %(K, slope))
    ax.plot(x_axis[lowbin:],  10**(intercept)*x_axis[lowbin:]**slope, linewidth=2.0, linestyle='--',color='black') 
    ax.plot(x_axis[lowbin:],  10**(fixintercept)*x_axis[lowbin:]**slope_fixed, linewidth=2.0, linestyle='--',color='red') 

    # but I don't want this value to change with number of samples
    # I want the width of the distribution.
    #std_err_yaxis/mean_y_axis
    #np.log10(std_err_yaxis)/np.log10(mean_y_axis)

    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.xlim((8,130))
    plt.title("Vmax Function normalized to a 151 km/s host")
    plt.xlabel('$V_{max}[km/s]$')
    plt.ylabel('$dN/dV_{sub} [(km/s)^{-1}] $')
    matplotlib.rcParams.update({'font.size': 15})

    ext=''
    if lx==13:
        ext='_'+str(lx)

    if version=='normal':
        plt.figtext(.15,.25, 'Vmax definition at z=0')	
        plt.savefig('hostfigs/HostHaloVmaxNgreater'+ext)

    if version=='peak':
        plt.figtext(.15,.25, 'Vmax definition at peak vmax')
        plt.savefig('hostfigs/HostHaloPeakVmaxNgreater'+ext)

    plt.close()


plotFits(lx=14, version='peak')
plotFits(lx=14, version='normal')
plotFits(lx=13, version='peak')
plotFits(lx=13, version='normal')




