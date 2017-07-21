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
    # xarray should be the same in all cases for me to average things
    """
    bins_per_dex = 5
    min_mass=7.0; max_mass=np.log10(cat.ix[haloid]['mgrav']/cat.h0)-1.5  #max_mass=10.0
    histrange = np.arange(min_mass,max_mass+.01,1./bins_per_dex) #determines subhalo mass 
    radius=factor*cat.ix[haloid]['rvir']  #  This is good without dividing by h0.
    masses_sub = np.array(cat.get_subhalos_within_halo(haloid, radius)['mgrav']/cat.h0)
    [x_array_sub, y_sub, xlabel, ylabel] = massfunction3(masses_sub, histrange)
    dNdM = y_sub*(target_halo_mass)/(cat.ix[haloid]['mgrav']/cat.h0)
    return 10**x_array_sub, dNdM


def shea_shmf(cat,haloid, field = False):
    bins_per_dex = 5
    min_mass = -4.5; max_mass = -1.5
    if field:
        min_mass = -3; max_mass = -1.3
    xdata = np.logspace(min_mass, max_mass, bins_per_dex*(max_mass - min_mass))
    hostmass = float(cat.ix[haloid]['mgrav']/cat.h0)
    histrange = np.log10(xdata * hostmass)
    masses_sub = np.array(cat.get_subhalos_within_halo(haloid)['mgrav']/cat.h0)
    [x_array_sub, y_sub, xlabel, ylabel] = massfunction3(masses_sub, histrange)
    return 10**x_array_sub/hostmass , y_sub


def make_fit_shea(y_axis, x_axis):
    y_values_full = np.array(y_axis)
    lengths = np.array([len(row) for row in y_values_full])
    maxlength = np.max(lengths)
    mean_y_axis = np.zeros(maxlength)
    for i in range(maxlength):
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
    slope_fixed = -1.9
    fixintercept = np.mean(np.log10(mean_y_axis)[lowbin:highbin] - slope_fixed*np.log10(x_axis)[lowbin:highbin])
    
    K = 10**intercept
    Kfix = 10**fixintercept
    return K, slope, Kfix, slope_fixed



def make_fit(y_values_full, target_halo_mass, x_axis):
    ###
    # I want to average the values of each column. if NaN, ignore it. only average over what you can.
    # To keep the number of bins per dex the same, I have to use variable length y_axis.
    ###
    y_values_full = np.array(y_values_full)
    lengths = np.array([len(row) for row in y_values_full])
    maxlength = np.max(lengths)
    mean_y_axis = np.zeros(maxlength)
    for i in range(maxlength):
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
    slope_fixed = -1.90
    fixintercept = np.mean(np.log10(mean_y_axis)[lowbin:highbin] - slope_fixed*np.log10(x_axis)[lowbin:highbin])
    
    K = 10**intercept/(target_halo_mass)
    Kfix = 10**fixintercept/(target_halo_mass)
    return K, slope, Kfix, slope_fixed


def get_field_halos(cat,hpath,hostID,mlow=10,mhigh=11.5):
    # get halos beyond host halo virial radius, but less than 
    # contamination radius
    hosts = cat.get_hosts()
    dists = dm.distance(cat.ix[hostID][['posX','posY','posZ']], hosts[['posX','posY','posZ']])
    contam_dist = dm.get_contam_dist(hpath)
    mask = dists < contam_dist
    mass_mask = (10**mlow < hosts['mgrav']/cat.h0) &(10**mhigh > hosts['mgrav']/cat.h0)
    return hosts[mask & mass_mask]


def field_shmf(cat, haloid, target_halo_mass=1e12, factor=1):
    """           
    @ return: x-axis array of M_sub. y-axis array of dN/dM_sub   
    """
    bins_per_dex = 5
    min_mass=7.5;max_mass=np.log10(cat.ix[haloid]['mgrav']/cat.h0)-1.5
    histrange = np.arange(min_mass,max_mass+.01,1./bins_per_dex) #determines subhalo mass range
    radius=factor*cat.ix[haloid]['rvir']
    masses_sub = np.array(cat.get_subhalos_within_halo(haloid,radius)['mgrav']/cat.h0)
    [x_array_sub, y_sub, xlabel, ylabel] = massfunction3(masses_sub, histrange)
    dNdM = y_sub*(target_halo_mass)/(cat.ix[haloid]['mgrav']/cat.h0)
    return 10**x_array_sub, dNdM

# use this function to plot SHMF for host halo, peak mass function for host halo
def plotFits(lx=14,target_halo_mass=1e12,factor=1):
    hpaths = dm.get_hpaths(field=False, lx=lx) #htils.get_all_halo_paths_lx(lx)
    f, (ax, ax2) = plt.subplots(2, 1)
    y_values_full = []; bh=0; mh=0; sh=0
    y_shea_full = []

    small=[];medium=[];large=[]; field=[]
    small_shea=[];medium_shea=[];large_shea=[]; field_shea=[]
    

    x_field_full = []; x_full = []
    

    for hpath, color in zip(hpaths[0:],colorlist[0:len(hpaths)]):
        snap_z0 = htils.get_numsnaps(hpath)-1
        cat=htils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID = htils.load_zoomid(hpath)

        # FIELD HALOS
        field_halos = get_field_halos(cat,hpath,hostID,mlow=10,mhigh=11.5)
        for i in range(len(field_halos)):
            field_halo=field_halos[i:i+1]
            x_axis_field, y_axis_field = field_shmf(cat,int(field_halo['id']), target_halo_mass, factor)
            xshea_field, yshea_field = shea_shmf(cat,int(field_halo['id']), field=True)
            field_shea.append(yshea_field)

            field.append(y_axis_field)
            if len(x_axis_field)>len(x_field_full):
                x_field_full = x_axis_field
        # END FIELD HALOS

        x_axis, y_axis = shmf(cat,hostID, target_halo_mass, factor)
        if len(x_axis)>len(x_full):
            x_full = x_axis
        y_values_full.append(y_axis)
        print 'done with cat', htils.hpath_catnum(hpath)
        ## Plot all the data    
        hostmass = float(cat.ix[hostID]['mgrav']/cat.h0)

        xshea, yshea = shea_shmf(cat,hostID)
        y_shea_full.append(yshea)

        if hostmass < 1.3e12:
            color='blue'
            small.append(y_axis)
            small_shea.append(yshea)
            sh+=1
            ax.plot(x_axis, .05*y_axis, linewidth=0.8,color=color)
        elif hostmass > 1.3e12 and hostmass < 1.6e12:
            color='grey'
            medium.append(y_axis)
            medium_shea.append(yshea)
            mh+=1
        else:
            color='red'
            large.append(y_axis)
            large_shea.append(yshea)
            bh+=1
        #ax.plot(x_axis, y_axis, linewidth=0.8,color=color)
    print sh, mh, bh, 'small medium big'

    #full
    K, slope, Kfix, slope_fixed = make_fit(y_values_full, target_halo_mass, x_full)
    ax.plot(x_axis,  Kfix*x_axis**slope_fixed * target_halo_mass, linewidth=2.0, linestyle='--',label='full')
    print 'full values K, slope:', Kfix, slope_fixed, K ,slope

    K, slope, Kfix, slope_fixed = make_fit_shea(y_shea_full, xshea)
    ax2.plot(xshea,  Kfix*xshea**slope_fixed, linewidth=2.0, linestyle='--',label='full')
    print 'Shea full values K, slope:', Kfix, slope_fixed, K ,slope

    # small
    K, slope, Kfix, slope_fixed = make_fit(small, target_halo_mass, x_full)
    ax.plot(x_axis,  Kfix*x_axis**slope_fixed * target_halo_mass, linewidth=2.0, linestyle='--',label='small')
    print 'small values K, slope:', Kfix, slope_fixed, K ,slope

    K, slope, Kfix, slope_fixed = make_fit_shea(small_shea, xshea)
    ax2.plot(xshea,  Kfix*xshea**slope_fixed, linewidth=2.0, linestyle='--',label='small')
    print 'Shea small values K, slope:', Kfix, slope_fixed, K ,slope


    # medium
    K, slope, Kfix, slope_fixed = make_fit(medium, target_halo_mass, x_full)
    ax.plot(x_axis,  Kfix*x_axis**slope_fixed * target_halo_mass, linewidth=2.0, linestyle='--',label='medium')
    print 'medium values K, slope:', Kfix, slope_fixed, K ,slope

    K, slope, Kfix, slope_fixed = make_fit_shea(medium_shea, xshea)
    ax2.plot(xshea,  Kfix*xshea**slope_fixed, linewidth=2.0, linestyle='--',label='medium')
    print 'Shea medium values K, slope:', Kfix, slope_fixed, K ,slope

    # large
    K, slope, Kfix, slope_fixed = make_fit(large, target_halo_mass, x_full)
    ax.plot(x_axis,  Kfix*x_axis**slope_fixed * target_halo_mass, linewidth=2.0, linestyle='--', label='large')
    print 'large values K, slope:', Kfix, slope_fixed, K ,slope

    K, slope, Kfix, slope_fixed = make_fit_shea(large_shea, xshea)
    ax2.plot(xshea,  Kfix*xshea**slope_fixed, linewidth=2.0, linestyle='--',label='large')
    print 'Shea large values K, slope:', Kfix, slope_fixed, K ,slope

    # field
    K, slope, Kfix, slope_fixed = make_fit(field, target_halo_mass, x_field_full)
    ax.plot(x_axis_field,  Kfix*x_axis_field**slope_fixed * target_halo_mass, linewidth=2.0, linestyle='--', label='field')
    print 'field values K, slope:', Kfix, slope_fixed, K ,slope

    K, slope, Kfix, slope_fixed = make_fit_shea(field_shea, xshea_field)
    ax2.plot(xshea_field,  Kfix*xshea_field**slope_fixed, linewidth=2.0, linestyle='--',label='field')
    print 'Shea field values K, slope:', Kfix, slope_fixed, K ,slope


    #plt.figtext(.15,.22, '$dn/dM_{sub} = %f \, M_{sub}^{%f}\,  M_{host}$' %(Kfix, slope_fixed))
    #plt.figtext(.15,.15, '$dn/dM_{sub} = %f \, M_{sub}^{%f}\,  M_{host}$' %(K, slope))

    #ax.plot(x_axis[lowbin:],  10**(fixintercept)*x_axis[lowbin:]**slope_fixed, linewidth=2.0, linestyle='--',color='red')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_yscale('log')


    plt.legend()
    plt.title("SHMF normalized to a 1e12 host")
    ax.set_xlabel('$M_{sub}[M_\odot /h]$')
    ax2.set_xlabel('$M_{sub} / M_{host}$')
    ax.set_ylabel('$dN/dM_{sub} [M_\odot^{-1}] $')
    ax2.set_ylabel('$dN/dM_{sub} [M_\odot^{-1}] $')
    matplotlib.rcParams.update({'font.size': 15})
    
    plt.figtext(.15,.25, 'Mgrav definition at z=0')
    plt.savefig('hostfigs/SHMF_Functional_Form')


plotFits()





"""


full values K, slope: 0.00119063900062 -1.9 0.000708473495181 -1.87408546014
field values K, slope: 0.00121676436775 -1.9 0.00446402108224 -1.96641408868
small values K, slope: 0.00129391962638 -1.9 0.000491090599445 -1.85050050829
medium values K, slope: 0.0011784583331 -1.9 0.000641663802632 -1.86930154755
large values K, slope: 0.00107283189486 -1.9 0.00149939179926 -1.91671074193




# X values here must be Msub / Mhost
Shea full values K, slope: 1.41933147462e-14 -1.9 1.82718765087e-14 -1.86587132337
Shea field values K, slope: 6.56672291297e-13 -1.9 5.25451238094e-13 -1.94046053144
Shea small values K, slope: 2.18001300895e-14 -1.9 3.13617808098e-14 -1.85086260401
Shea medium values K, slope: 1.37775947571e-14 -1.9 2.04936502686e-14 -1.84635016011
Shea large values K, slope: 9.85130448142e-15 -1.9 6.87983239731e-15 -1.94850716391




Final confirmation: use my host halo mass function to plot N(>vmax) for a 1e10.5 mass host
use Shea's mass function for a 1e10.5 mass host
plot against the real data from field halos (integrate my mass function values from above)

I believe shea's function will be very off, but mine will be fairly close




field with sub-subs: K = .001313, slope = -1.90
Without:             K = .001217
.926 is the ratio

"""


# plot-elvis-mycode. look at dn/dm vs msub/mhost


