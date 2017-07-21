import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import haloutils as htils
from scipy import stats
import methods
import MTanalysis_field as mtaf
import DwarfMethods as dm


colorlist = ['black','red','green','blue','orange','yellow','purple','grey','brown','lawngreen','steelblue','peru','khaki','indigo','silver','maroon','sage','midnightblue','orangered','gold','deeppink','gainsboro','salmon','seagreen','slateblue','sandybrown','goldenrod','pink','cyan','darkslategrey']


def massfunction3(masses,histrange):
    """                                        
    Produce mass function data for dN/dM. dN is number of halos in mass 
    interval dM. numbins needs to be len(histrange)-1. Outdated, but too many places to fix.
    @param masses: list of all masses in Msun
    @param histrange: binning range in logspace
    @return: [x-axis in log10(mass), y-axis in dN/dM, xlabel, ylabel]
    """
    hist, r_array = np.histogram(np.log10(masses), bins=histrange)
    #print hist, 'histogram'
    x_array = dm.getMidpoints(r_array)
    dM = 10.**r_array[1:]-10.**r_array[0:-1] #Mass size of bins in non-log space.
    dNdM = hist/dM
    return [x_array, dNdM, '$\log_{10}(M [M_\odot)])$', '$\log_{10}(dN/dM) [M_\odot^{-1}]$']


# given cat, haloid, find shmf
# find all subhalos from host halos in a mass range given by low and high                
# host come from catalogue cat, and SHMF normalized to a host                            
# of mass target_halo_mass            

# mgrav of subhalos and host at z=0. Use for GarrisonKimmel4
def shmf(cat, haloid, target_halo_mass=1e12, factor=1):
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

def m200_function(cat, haloid, target_halo_mass=1e12):
    """                                                                                  
    @ return: x-axis array of M_sub. y-axis array of dN/dM_sub   
    """
    bins_per_dex = 5
    min_mass=7.5;max_mass=np.log10(cat.ix[haloid]['altm2']/cat.h0)-1.5
    histrange = np.arange(min_mass,max_mass+.01,1./bins_per_dex) #determines subhalo mass range
    masses_sub = np.array(cat.get_subhalos_within_halo(haloid)['altm2']/cat.h0)
    [x_array_sub, y_sub, xlabel, ylabel] = massfunction3(masses_sub, histrange)
    dNdM = y_sub*(target_halo_mass)/(cat.ix[haloid]['altm2']/cat.h0)
    return 10**x_array_sub, dNdM

# mgrav of subhalos at max mass, mgrav of hosts at z=0. Use for GarrisonKimmel16
def max_mass_function(hpath,cat,haloid,target_halo_mass=1e12):
    data = dm.get_extant_data(hpath, field=True)
    mask = data['field_rsid']==haloid
    data = data[mask]

    max_masses = data['max_mass']/cat.h0
    bins_per_dex = 5
    min_mass=7.5;max_mass=np.log10(cat.ix[haloid]['mgrav']/cat.h0)-1.5
    histrange = np.arange(min_mass,max_mass+.01,1./bins_per_dex) #determines subhalo mass range
    [x_array_sub, y_sub, xlabel, ylabel] = massfunction3(max_masses, histrange)
    dNdM = y_sub*(target_halo_mass)/(cat.ix[haloid]['mgrav']/cat.h0)
    return 10**x_array_sub, dNdM


# this is for moster
# subhalos must be taken at infall mass 200
# host field halos taken as m200 at z=0
def infall_mass200_function(hpath,cat,haloid,target_halo_mass=1e12):
    data = dm.get_extant_data(hpath, field=True)
    mask = data['field_rsid']==haloid
    data = data[mask]

    infall_m200 = data['infall_mass200']/cat.h0
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
    # first get host halo mass in 350 crit
    data = dm.get_extant_data(hpath, field=True)
    mask = data['field_rsid']==haloid
    data = data[mask]
    #hostmass = np.array(data['host_m350z0'])[0]/cat.h0
    hostmassNFW = np.array(data['host_m350NFWz0'])[0]/cat.h0

    maxmass = data['max_mass350']/cat.h0
    maxmassNFW = data['max_mass350NFW']/cat.h0

    # sort out the NFW mass values that are < 0
    badmask = maxmassNFW < 0    
    peak_masses = data['max_mass']/cat.h0
    meanratio = np.mean(  maxmassNFW[~badmask]/peak_masses[~badmask])
    maxmassNFW[badmask] = meanratio*peak_masses[badmask]

    bins_per_dex = 5
    min_mass=7.5;max_mass=np.log10(hostmassNFW)-1.5
    histrange = np.arange(min_mass,max_mass+.01,1./bins_per_dex) #determines subhalo mass range
    [x_array_sub, y_sub, xlabel, ylabel] = massfunction3(maxmassNFW, histrange)
    dNdM = y_sub*(target_halo_mass)/(hostmassNFW)
    return 10**x_array_sub, dNdM


def get_field_halos(cat,hpath,hostID,mlow=10,mhigh=11.5):
    # get halos beyond host halo virial radius, but less than 
    # contamination radius
    hosts = cat.get_hosts()
    dists = dm.distance(cat.ix[hostID][['posX','posY','posZ']], hosts[['posX','posY','posZ']])
    contam_dist = dm.get_contam_dist(hpath)
    mask = dists < contam_dist
    mass_mask = (10**mlow < np.array(hosts['mgrav']/cat.h0)) &(10**mhigh > np.array(hosts['mgrav']/cat.h0))
    return hosts[mask & mass_mask]


# Plot shmf of field halos, and peak mass function of field halos
def plotFits(version='max',lx=14,target_halo_mass=1e11, factor=1):
    hpaths = dm.get_hpaths(field=True)  #htils.get_all_halo_paths_lx(lx)[0:30]
    plt.figure()
    ax = plt.subplot(111)
    y_values_full = []
    x_axis_full = []
    for hpath, color in zip(hpaths[0:],colorlist[0:len(hpaths)]):
        snap_z0 = htils.get_numsnaps(hpath)-1
        cat=htils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID = htils.load_zoomid(hpath)
        field_halos = get_field_halos(cat,hpath,hostID,mlow=10,mhigh=11.5)  #mlow = np.log10(5*10**10), mhigh = np.log10(5*10**11))
        for i in range(len(field_halos)):
            field_halo=field_halos[i:i+1]
            if version=='normal':
                x_axis, y_axis = shmf(cat,int(field_halo['id']), target_halo_mass, factor)
            if version=='max':
                x_axis, y_axis = max_mass_function(hpath,cat,int(field_halo['id']),target_halo_mass)
	    if version=='m200':
		x_axis, y_axis = m200_function(cat,int(field_halo['id']),target_halo_mass)
	    if version=='m200_infall':
		x_axis, y_axis = infall_mass200_function(hpath,cat,int(field_halo['id']),target_halo_mass)
            if version=='m350_max':
                x_axis, y_axis = max_mass350_function(hpath,cat,int(field_halo['id']),target_halo_mass)

            y_values_full.append(y_axis)
            if len(x_axis)>len(x_axis_full):
                x_axis_full = x_axis
            ## Plot all the data         
            ax.plot(x_axis, y_axis, linewidth=0.8,color=color)
        print 'done with cat', htils.hpath_catnum(hpath)

    y_values_full = np.array(y_values_full)
    print len(y_values_full), 'number of field halos in total'
    ###
    # I want to average the values of each column. if NaN, ignore it. only average over what you can.
    ###
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

    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x_axis_full)[lowbin:highbin],np.log10(mean_y_axis)[lowbin:highbin])
    # FIXED INTERCEPT FIT
    slope_fixed = -1.90
    fixintercept = np.mean(np.log10(mean_y_axis)[lowbin:highbin] - slope_fixed*np.log10(x_axis_full)[lowbin:highbin])

    print 'alpha = slope = ', slope
    K = 10**intercept/(target_halo_mass)
    Kfix = 10**fixintercept/(target_halo_mass)
    print 'K = ', K, 'Kfixed', Kfix
    print std_err, 'std err'
    #"std_err" is 0.13228756555322888 is the standard error of the slope coefficient, not of the fitting error.


    plt.figtext(.15,.22, '$dn/dM_{sub} = %f \, M_{sub}^{%f}\,  M_{host}$' %(Kfix, slope_fixed))
    plt.figtext(.15,.15, '$dn/dM_{sub} = %f \, M_{sub}^{%f}\,  M_{host}$' %(K, slope))

    ax.plot(x_axis_full, mean_y_axis, linewidth=2.5,linestyle='-',color='blue',label='field halo mean data')
    ax.plot(x_axis_full[lowbin:],  10**(intercept)*x_axis_full[lowbin:]**slope, linewidth=2.0, linestyle='--',color='black',label='field halo fit')
    
    ax.plot(x_axis_full[lowbin:],  10**(fixintercept)*x_axis_full[lowbin:]**slope_fixed, linewidth=2.0, linestyle='--',color='red',label='forced slope') 

    ax.set_xscale('log')
    ax.set_yscale('log')  
    plt.legend(loc='upper right')
    plt.title("SHMF normalized to a 1e11 host")
    plt.xlabel('$M_{sub}[M_\odot /h]$')
    plt.ylabel('$dN/dM_{sub} [M_\odot^{-1}] $')
    matplotlib.rcParams.update({'font.size': 15})


    extra = '_LMCsized'
    if version =='normal':
        #ax.plot(x_axis_full[lowbin:],  0.00102*target_halo_mass*x_axis_full[lowbin:]**-1.895046, linewidth=2.0, linestyle='--',color='red',label='host halo fit')
	plt.figtext(.15,.25, 'Mgrav definition at z=0')	
        plt.savefig('fieldfigs/FieldHaloSHMF_%d-%d_rvir_LMCsized' %(int(factor), round(10*(factor-int(factor)))))
    if version =='max':
        ax.plot(x_axis_full[lowbin:],  0.002734*target_halo_mass*x_axis_full[lowbin:]**-1.883824, linewidth=2.0, linestyle='--',color='red',label='host halo fit')
        plt.figtext(.15,.25, 'Mgrav definition at max mass')
        plt.savefig('fieldfigs/FieldHaloMaxSHMF'+extra)
    if version =='m200':
        plt.figtext(.15,.25, 'M200 definition at z=0')
	plt.savefig('fieldfigs/FieldHaloM200z0'+extra)
    if version =='m200_infall':
        plt.figtext(.15,.25, 'M200 definition at infall mass')
	plt.savefig('fieldfigs/FieldHaloM200_infall'+extra)
    if version =='m350_max':
        plt.figtext(.15,.25, 'M350 definition at max mass')
	plt.savefig('fieldfigs/FieldHaloM350NFW_max'+extra)
    plt.close()




def jackknife_shmf(version='max',lx=14,target_halo_mass=1e11, factor=1, exclude=None):
    hpaths = dm.get_hpaths(field=True)  #htils.get_all_halo_paths_lx(lx)[0:30]
    if exclude is not None:
        hpaths = np.append(hpaths[0:exclude], hpaths[exclude+1:])
    y_values_full = []
    x_axis_full = []
    for hpath, color in zip(hpaths[0:],colorlist[0:len(hpaths)]):
        snap_z0 = htils.get_numsnaps(hpath)-1
        cat=htils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID = htils.load_zoomid(hpath)
        field_halos = get_field_halos(cat,hpath,hostID,mlow=10,mhigh=11.5)
        for i in range(len(field_halos)):
            field_halo=field_halos[i:i+1]
            if version=='normal':
                x_axis, y_axis = shmf(cat,int(field_halo['id']), target_halo_mass, factor)
            if version=='max':
                x_axis, y_axis = max_mass_function(hpath,cat,int(field_halo['id']),target_halo_mass)
	    if version=='m200':
		x_axis, y_axis = m200_function(cat,int(field_halo['id']),target_halo_mass)
	    if version=='m200_infall':
		x_axis, y_axis = infall_mass200_function(hpath,cat,int(field_halo['id']),target_halo_mass)
            if version=='m350_max':
                x_axis, y_axis = max_mass350_function(hpath,cat,int(field_halo['id']),target_halo_mass)

            y_values_full.append(y_axis)
            if len(x_axis)>len(x_axis_full):
                x_axis_full = x_axis
        print 'done with cat', htils.hpath_catnum(hpath)

    y_values_full = np.array(y_values_full)
    print len(y_values_full), 'number of field halos in total'
    ###
    # I want to average the values of each column. if NaN, ignore it. only average over what you can.
    ###
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

    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x_axis_full)[lowbin:highbin],np.log10(mean_y_axis)[lowbin:highbin])
    # FIXED INTERCEPT FIT
    slope_fixed = -1.88
    fixintercept = np.mean(np.log10(mean_y_axis)[lowbin:highbin] - slope_fixed*np.log10(x_axis_full)[lowbin:highbin])


    K = 10**intercept/(target_halo_mass)
    Kfix = 10**fixintercept/(target_halo_mass)
    print Kfix, 'Kfix for 1.88'

    slope_fixed = -1.75
    fixintercept = np.mean(np.log10(mean_y_axis)[lowbin:highbin] - slope_fixed*np.log10(x_axis_full)[lowbin:highbin])
    Kfix = 10**fixintercept/(target_halo_mass)
    print Kfix, 'Kfix for 1.75'

    return slope, K


"""
s_list=[]; k_list=[]
for i in range(32):
    slope,K = jackknife_shmf(lx=14, version='m200_infall',exclude=i)
    print slope, K, i
    s_list.append(slope); k_list.append(K)
print s_list
print k_list
"""



#plotFits(lx=14,version='normal')
#plotFits(lx=14,version='max')
#plotFits(lx=14,version='m200_infall')
#plotFits(lx=14,version='m350_max')
#plotFits(lx=14,version='m200')



"""
plotFits(lx=14,version='normal',factor=0.2)
plotFits(lx=14,version='normal',factor=0.3)
plotFits(lx=14,version='normal',factor=0.4)
plotFits(lx=14,version='normal',factor=0.5)
plotFits(lx=14,version='normal',factor=0.6)
"""
"""
plotFits(lx=14,version='normal',factor=0.7)
plotFits(lx=14,version='normal',factor=0.8)
plotFits(lx=14,version='normal',factor=0.9)
plotFits(lx=14,version='normal',factor=1.0)
plotFits(lx=14,version='normal',factor=1.1)
"""
"""
plotFits(lx=14,version='normal',factor=1.2)
plotFits(lx=14,version='normal',factor=1.3)
plotFits(lx=14,version='normal',factor=1.4)
plotFits(lx=14,version='normal',factor=1.5)
plotFits(lx=14,version='normal',factor=1.6)
"""

"""
plotFits(lx=14,version='normal',factor=1.7)
plotFits(lx=14,version='normal',factor=1.8)
plotFits(lx=14,version='normal',factor=1.9)
plotFits(lx=14,version='normal',factor=2.0)
"""
