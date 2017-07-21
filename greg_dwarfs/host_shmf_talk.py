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


def getMidpoints(bins):
    """                                                                                                             
    Given a range of cut-off values, find the midpoints.                                                            
    @return: array of length len(bins)-1                                                                            
    """
    spacing = bins[1:]-bins[:-1]
    return bins[:-1]+spacing/2.0


def massfunction3(masses, histrange):
    """                                        
    Produce mass function data for dN/dM. dN is number of halos in mass                     
    interval dM. numbins needs to be len(histrange)-1. Outdated, but too many places to fix.                        
    @param masses: list of all masses in Msun
    @param histrange: binning range in logspace    
    @return: [x-axis in log10(mass), y-axis in dN/dM, xlabel, ylabel]
    """
    hist, r_array = np.histogram(np.log10(masses), bins=histrange)
    #print hist, 'histogram'
    x_array = getMidpoints(r_array)
    dM = 10.**r_array[1:]-10.**r_array[0:-1] #Mass size of bins in non-log space.
    dNdM = hist/dM
    return [x_array, dNdM, '$\log_{10}(M [M_\odot)])$', '$\log_{10}(dN/dM) [M_\odot^{-1}]$']


# given cat, haloid, find shmf
# find all subhalos from host halos in a mass range given by low and high                
# host come from catalogue cat, and SHMF normalized to a host                            
# of mass target_halo_mass            

def shmf(cat, haloid, target_halo_mass=1e12,factor=1):
    """                                                                                  
    @ return: x-axis array of M_sub. y-axis array of dN/dM_sub   
    """
    bins_per_dex = 5
    min_mass=7.0; max_mass=np.log10(cat.ix[haloid]['mgrav']/cat.h0)-1.5  #max_mass=10.0
    histrange = np.arange(min_mass,max_mass+.01,1./bins_per_dex) #determines subhalo mass 
    radius=factor*cat.ix[haloid]['rvir']  #  This is good without dividing by h0.
    masses_sub = np.array(cat.get_subhalos_within_halo(haloid, radius)['mgrav']/cat.h0)
    [x_array_sub, y_sub, xlabel, ylabel] = massfunction3(masses_sub, histrange)
    dNdM = y_sub*(target_halo_mass)/(cat.ix[haloid]['mgrav']/cat.h0)
    return 10**x_array_sub, dNdM


def get_mass_string(v):
    if v=='peak':
        return 'max_mass'
    if v=='peak_full':
        return 'max_mass_full'
    if v=='peak_half':
        return 'max_mass_half'
    if v=='peak_third':
        return 'max_mass_third'
    if v=='peak_fourth':
        return 'max_mass_fourth'
    
def peak_mass_function(hpath,cat,haloid,target_halo_mass=1e12, version='peak',lx=14):
    if lx==13:
        AE = mta.AllExtantData()
        dataE = AE.read(hpath)
    else:
        dataE = dm.get_extant_data(hpath,False)
    mass_string = get_mass_string(version)
    peak_masses = dataE[mass_string]/cat.h0

    bins_per_dex = 5
    min_mass=7.5; max_mass=np.log10(cat.ix[haloid]['mgrav']/cat.h0)-1.5 # peak mass must be 7.5, or else dataE doesn't hav enough values
    if lx==13:
        min_mass=8
    histrange = np.arange(min_mass,max_mass+.01,1./bins_per_dex)
    [x_array_sub, y_sub, xlabel, ylabel] = massfunction3(peak_masses, histrange)
    dNdM = y_sub*(target_halo_mass)/(cat.ix[haloid]['mgrav']/cat.h0)
    return 10**x_array_sub, dNdM


##### START MODIFYING FIELD FUNCTIONS TO HOST HALO #######
# add to the catalogs - m200 infall
# m350 max mass
# m200 peak
# - means re-running all of the extant+destroyed catalogs?
# or can I do it faster, and append data?
# create function to keep appending data to it all
# only need to re-run the extant catalogs


# this is the z=0 M 200 subhalo mass function
def m200_function(cat, haloid, target_halo_mass=1e12):
    """                                                                                  
    @ return: x-axis array of M_sub. y-axis array of dN/dM_sub   
    """
    bins_per_dex = 5
    min_mass=7.0;max_mass=np.log10(cat.ix[haloid]['altm2']/cat.h0)-1.5
    histrange = np.arange(min_mass,max_mass+.01,1./bins_per_dex) #determines subhalo mass range
    masses_sub = np.array(cat.get_subhalos_within_halo(haloid)['altm2']/cat.h0)
    [x_array_sub, y_sub, xlabel, ylabel] = massfunction3(masses_sub, histrange)
    dNdM = y_sub*(target_halo_mass)/(cat.ix[haloid]['altm2']/cat.h0)
    return 10**x_array_sub, dNdM


# this is for moster
# subhalos must be taken at infall mass 200
# host field halos taken as m200 at z=0
def infall_mass200_function(hpath,cat,haloid,target_halo_mass=1e12):
    dataE = dm.get_extant_data(hpath,False)
    infall_m200 = dataE['infall_mass200']/cat.h0
    bins_per_dex = 5
    min_mass=7.5;max_mass=np.log10(cat.ix[haloid]['altm2']/cat.h0)-1.5
    histrange = np.arange(min_mass,max_mass+.01,1./bins_per_dex) #determines subhalo mass range
    [x_array_sub, y_sub, xlabel, ylabel] = massfunction3(infall_m200, histrange)
    dNdM = y_sub*(target_halo_mass)/(cat.ix[haloid]['altm2']/cat.h0)
    return 10**x_array_sub, dNdM


# use this for Brook
# subhalos at max mass 350
# host halos at z=0 M350
def max_mass350_function(hpath,cat,haloid,target_halo_mass=1e12):
    # first get host halo mass in 350 crit.     # I need the m350 of all host halos!!
    hostM350 = mtaf.HostHaloM350()
    m350data = hostM350.read(hpath)
    hostmassNFW = float(m350data['m350NFW_host']/cat.h0)
    print np.log10(hostmassNFW)
    
    dataE = dm.get_extant_data(hpath,False)  
    maxmass = dataE['max_mass350']/cat.h0
    maxmassNFW = dataE['max_mass350NFW']/cat.h0

    badmask = maxmassNFW < 0    
    peak_masses = dataE['max_mass']/cat.h0
    ratios = maxmassNFW[~badmask]/peak_masses[~badmask]
    meanratio = np.mean(ratios) 
    maxmassNFW[badmask] = meanratio*peak_masses[badmask]
    

    bins_per_dex = 5
    min_mass=7.5; max_mass=np.log10(cat.ix[haloid]['mgrav']/cat.h0)-1.5
    histrange = np.arange(min_mass,max_mass+.01,1./bins_per_dex) #determines subhalo mass range
    [x_array_sub, y_sub, xlabel, ylabel] = massfunction3(maxmassNFW, histrange)
    dNdM = y_sub*(target_halo_mass)/(hostmassNFW)
    return 10**x_array_sub, dNdM

###### END CHANGES ######



# use this function to plot SHMF for host halo, peak mass function for host halo
def plotFits(version='peak',lx=14,target_halo_mass=1e12,factor=1):
    hpaths = dm.get_hpaths(field=False, lx=lx)  #htils.get_all_halo_paths_lx(lx)
    plt.figure()
    ax = plt.subplot(111)
    y_values_full = []; x_full=[]
    for hpath, color in zip(hpaths[0:],colorlist[0:len(hpaths)]):
        snap_z0 = htils.get_numsnaps(hpath)-1
        cat=htils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID = htils.load_zoomid(hpath)
        if version=='normal':
            x_axis, y_axis = shmf(cat,hostID, target_halo_mass, factor)
        if version[0:4]=='peak':
            x_axis, y_axis = peak_mass_function(hpath,cat,hostID,target_halo_mass,version,lx)
        if version=='m200':
            x_axis, y_axis = m200_function(cat, hostID,target_halo_mass)
        if version=='m200_infall':
            x_axis, y_axis = infall_mass200_function(hpath,cat,hostID, target_halo_mass)
        if version=='m350_max':
            x_axis, y_axis = max_mass350_function(hpath,cat,hostID, target_halo_mass)
            
        y_values_full.append(y_axis)
        if len(x_axis)>len(x_full):
            x_full = x_axis
        print 'done with cat', htils.hpath_catnum(hpath)
        """
        lowbin = 0
        highbin = len(x_axis)-2
        slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x_axis)[lowbin:highbin],np.log10(y_axis)[lowbin:highbin])
        print 'alpha = slope = ', slope
        print intercept, 'intercept'
        print 'K = ',10**intercept/(10**12)
        print std_err, 'std err'
        """
        ## Plot all the data 
        mask0 = y_axis == 0
        ax.plot(x_axis[~mask0], y_axis[~mask0], linewidth=1.0,color=color)
        #ax.plot(x_axis[lowbin:],  10**(intercept)*x_axis[lowbin:]**slope, linewidth=1.0, linestyle='--') 


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

    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x_full)[lowbin:highbin],np.log10(mean_y_axis)[lowbin:highbin])
    # FIXED INTERCEPT FIT
    slope_fixed = -1.90
    fixintercept = np.mean(np.log10(mean_y_axis)[lowbin:highbin] - slope_fixed*np.log10(x_full)[lowbin:highbin])
    
    print 'alpha = slope = ', slope
    K = 10**intercept/(target_halo_mass)
    Kfix = 10**fixintercept/(target_halo_mass)
    print 'K = ', K, 'Kfixed', Kfix
    print std_err, 'std err'
    #plt.figtext(.15,.22, '$dn/dM_{sub} = %f \, M_{sub}^{%f}\,  M_{host}$' %(Kfix, slope_fixed))
    plt.figtext(.19,.2, '$\mathrm{dn/dM_{sub} = %.4f \, M_{sub}^{%.2f}\,  M_{host}}$' %(K, slope),fontsize=22)
    #ax.plot(x_full[lowbin:],  10**(intercept)*x_full[lowbin:]**slope, linewidth=2.0, linestyle='--',color='black') 
    ax.plot(x_full[lowbin:],  10**(fixintercept)*x_full[lowbin:]**slope_fixed, linewidth=4.0, linestyle='--',color='red') 

    # but I don't want this value to change with number of samples
    # I want the width of the distribution.
    #std_err_yaxis/mean_y_axis
    #np.log10(std_err_yaxis)/np.log10(mean_y_axis)

    ax.set_xscale('log')
    ax.set_yscale('log')        
    #plt.legend()
    #plt.title("SHMF normalized to a 1e12 host")
    plt.xlabel('$M_{\mathrm{sub}}^{\mathrm{peak}} [\mathrm{M_\odot /h}]$', fontsize=25)
    plt.ylabel('$\mathrm{dN/dM_{sub} [\mathrm{M_\odot^{-1}}]}$', fontsize=25)
    matplotlib.rcParams.update({'font.size': 20})
    ax.set_ylim((10**-11, 10**-4.3))
    ax.set_xlim((10**7.5, 10**11))

    plt.gcf().subplots_adjust(bottom=0.17,left=0.18, top=0.96,right=0.94)
    print 'adjusted topp'
    ext='TALK'
    if lx==13:
        ext='_'+str(lx)

    if version=='normal':
        plt.figtext(.15,.25, 'Mgrav definition at z=0')	
        plt.savefig('hostfigs/HostHaloSHMF_%d-%d_rvir' %(int(factor), round(10*(factor-int(factor)))) + ext )

    if version=='peak':
        #plt.figtext(.15,.25, 'Mgrav definition at max mass')
        plt.savefig('hostfigs/HostHaloPeakSHMF'+ext)


    if version=='peak_full':
        plt.figtext(.15,.25, 'Max mass full')
        plt.savefig('hostfigs/HostHaloPeakSHMF_full'+ext)
    if version=='peak_half':
        plt.figtext(.15,.25, 'Max mass half')
        plt.savefig('hostfigs/HostHaloPeakSHMF_half'+ext)
    if version=='peak_third':
        plt.figtext(.15,.25, 'Max mass third')
        plt.savefig('hostfigs/HostHaloPeakSHMF_third'+ext)
    if version=='peak_fourth':
        plt.figtext(.15,.25, 'Max mass fourth')
        plt.savefig('hostfigs/HostHaloPeakSHMF_fourth'+ext)

    if version =='m200':
        plt.figtext(.15,.25, 'M200 definition at z=0')
	plt.savefig('hostfigs/HostHaloM200z0'+ext)

    if version =='m200_infall':
        plt.figtext(.15,.25, 'M200 definition at infall mass')
	plt.savefig('hostfigs/HostHaloM200_infall'+ext)

    if version =='m350_max':
        plt.figtext(.15,.25, 'M350 definition at max mass')
	plt.savefig('hostfigs/HostHaloM350NFW_max'+ext)
    plt.close()




plotFits(version='peak')
