import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import haloutils as htils
from scipy import stats
import methods
#import MTanalysis3 as mta
from PlotParams import *
from scipy import interpolate
#import MTanalysis_field as mtaf
import DwarfMethods as dm
from PlotParams import *


## this file used to get the infall time distributions for moster et all
# abundance matching


# want to get all infall times of current z=0 satellites, of field and host
# plot infall time distribution, averaged over all 30 halos.
# difference in field vs host infall time distribution??


# divide into low mass and high mass systems, does it change much?
# divide into MW vs field halo hosts, does it change much?
## ANSWER: it does not depend on mass significantly.

## Why a bi-modal distribution? 

# if it doesn't change much by subhalo peak mass, then I can randomly assign
# stellar masses by sampling SHMF, then infall distribution, then assigning
# stellar mass from abundance matching assuming a subhalo mass and infall time.
# TEST: does it change it much?
# 

# even when binning by lookback time and by infall scale, you get the bimodal distribution.
# something is weird, but it doesn't affect things that much, and I don't want to do it!
def get_infall_distr(hpath, field=True):
    dataE = dm.get_extant_data(hpath, field)
    dataE = dataE[dataE['infall_snap']>0]
    
    peak_masses = dataE['max_mass']/0.67
    # mask for lower mass subhalos
    #mask = np.log10(peak_masses)> 9.5
    #dataE = dataE[mask]
    infall_snaps = dataE['infall_snap']   #'infall_scale'
    bins = np.arange(0,322)
    hist, snaps = np.histogram(infall_snaps, bins=bins,normed=True)

    # tmp to get mean infall time
    #print infall_snaps, 'infall snaps'
    if len(infall_snaps)>0:
        infall_scales =  htils.get_scale_snap(hpath,infall_snaps)
        infall_times = dm.lookback(infall_scales)
    else:
        infall_times = []
    """
    # histogram of times
    bins = np.arange(0,14.001,14/75.)  # 75 time bins in Gyrs
    infall_scales =  htils.get_scale_snap(hpath,infall_snaps)
    infall_times = dm.lookback(infall_scales)
    hist_times, times = np.histogram(infall_times, bins=bins,normed=True)

    # histogram of scale factors
    bins = np.arange(0,1.001,1./75.)
    hist_scales, scales = np.histogram(infall_scales, bins=bins,normed=True)
    """
    # if snapshots per simulation vary, need to convert to scale factors here
    return hist, snaps[0:-1], infall_times #, hist_times, times[0:-1], hist_scales, scales[0:-1] # last snapshot is upper bound of hist bounds.
    

def load_infall_distr(field=True):
    path = '/nfs/blank/h4231/gdooley/Dropbox/DwarfsOfDwarfs/code/SHMF/infall_distr.npy'
    if field:
        path = '/nfs/blank/h4231/gdooley/Dropbox/DwarfsOfDwarfs/code/SHMF/infall_distr_field.npy'
    snaps, mean_distr = np.load(path)
    return snaps,mean_distr


# x-axis is the infall snap
# y-axis is number of counts
def infall_distribution_all(field=True, exclude=None):
    #hpaths = htils.get_all_halo_paths_lx(lx=14)[11:30]
    hpaths = dm.get_hpaths(field)
    if exclude is not None:
        hpaths = np.append(hpaths[0:exclude], hpaths[exclude+1:])

    y_values_full = []; itimes=[]
    distrs = np.ndarray([len(hpaths),321]); badrows = []
    for hpath, color,i in zip(hpaths[0:],colorlist[0:len(hpaths)],range(len(hpaths))):
        #print htils.hpath_catnum(hpath),' cat num', i, 'value of i'        
        hist, snaps, infall_times = get_infall_distr(hpath,field)
        distrs[i] = hist
        itimes.append(infall_times)
        if np.isnan(distrs[i][0]):
            badrows.append(i)
        
        #tmp=np.repeat(snaps,2)
        #plt.plot(np.append(tmp[1:],tmp[-1]+1), np.repeat(hist,2),linewidth=1.5)
        #plt.savefig('infall_dist_%d' %i)
        #plt.close()
    x_axis = snaps

    distrs = np.delete(distrs,badrows,0)
    itimes = np.delete(itimes,badrows,0)
    itimes = np.array([item for row in itimes for item in row])

    mean_distr = np.mean(distrs,axis=0)
    #print mean_distr, 'mean distr'
    mean_itime = np.mean(itimes)
    print mean_itime #, 'mean infall time. Field =', field




    """
    tmp=np.repeat(x_axis,2)
    # make it a histogram. do the doubling of the data thing
    plt.figure()
    ax = plt.subplot(111)
    ax.plot(np.append(tmp[1:],tmp[-1]+1), np.repeat(mean_distr,2),linewidth=1.5)
    # also want to save data to file for sampling it later
    fname = 'infall_distr'
    if field:
        fname = 'infall_distr_field'
    np.save(fname,np.array([x_axis, mean_distr]))
    plt.savefig(fname)
    plt.close()
    """


# use the infall distribution to generate a random list of infall times 
def random_infall_scale(tck, N=1, offset=0):
    if N==0:
        return np.array([])    
    # convert randnums to scale factors    
    randnums = np.random.rand(N)
    a = interpolate.splev(randnums,tck)  #given cumulative distr, what is scalefactor?
    #print zip(randnums, a)
 
    if offset != 0:
        # convert a to time
        times = dm.lookback(a)
        # add offset
        times += offset
        # convert back to time
        a = np.array([dm.time_to_scale(time) for time in times])
    return a


def get_infall_tck(field=True):
    snaps, mean_distr = load_infall_distr(field)
    cum_distr = np.cumsum(mean_distr)
    
    # convert snaps to scale factors
    elements, snaps = np.unique(cum_distr, return_index=True)
    hpath = htils.catnum_hpath(1,lx=14)
    scales = htils.get_scale_snap(hpath,snaps)
    
    # create interpolation from cum_distr to scale factor
    tck = interpolate.splrep(elements,scales, k=1)
    return tck
    


def plot_infall_distribution(field=True, exclude=None):
    plt.figure()
    bins = np.arange(0,13.801,13.8/22.)  # 75 time bins in Gyrs

    #hpaths = htils.get_all_halo_paths_lx(lx=14)[11:30]
    hpaths = dm.get_hpaths(field)
    field = False
    if exclude is not None:
        hpaths = np.append(hpaths[0:exclude], hpaths[exclude+1:])

    hist_sum = []; itimes=[]
    distrs = np.ndarray([len(hpaths),321]); badrows = []
    for hpath, color,i in zip(hpaths[0:],colorlist[0:len(hpaths)],range(len(hpaths))):
        #print htils.hpath_catnum(hpath),' cat num', i, 'value of i'        
        hist, snaps, infall_times = get_infall_distr(hpath,field)
        distrs[i] = hist #np.histogram(infall_times, bins=bins,normed=True) 
        itimes.append(infall_times)
        if np.isnan(distrs[i][0]):
            badrows.append(i)

    distrs = np.delete(distrs,badrows,0)
    itimes = np.delete(itimes,badrows,0)
    itimes = np.array([item for row in itimes for item in row])

    print len(infall_times), len(itimes)
    hist_times, times = np.histogram(itimes, bins=bins,normed=True)
    x_axis = np.repeat(times,2)  #dm.getMidpoints(times)
    y_axis = np.repeat(hist_times,2)  # hist_times
    plt.plot(x_axis[1:-1], y_axis,lw=linewidth,ls='-',color='blue', label='MW-sized halos')

    #x_axis = dm.getMidpoints(times)
    #y_axis = hist_times
    #plt.plot(x_axis, y_axis,lw=linewidth,ls='-',color='blue', label='MW-sized halos')



    field = True
    hpaths = dm.get_hpaths(field)
    if exclude is not None:
        hpaths = np.append(hpaths[0:exclude], hpaths[exclude+1:])

    y_values_full = []; itimes=[]
    distrs = np.ndarray([len(hpaths),321]); badrows = []
    for hpath, color,i in zip(hpaths[0:],colorlist[0:len(hpaths)],range(len(hpaths))):
        #print htils.hpath_catnum(hpath),' cat num', i, 'value of i'        
        hist, snaps, infall_times = get_infall_distr(hpath,field)
        distrs[i] = hist
        itimes.append(infall_times)
        if np.isnan(distrs[i][0]):
            badrows.append(i)

    distrs = np.delete(distrs,badrows,0)
    itimes = np.delete(itimes,badrows,0)
    itimes = np.array([item for row in itimes for item in row])



    print len(infall_times), len(itimes)
    hist_times, times = np.histogram(itimes, bins=bins,normed=True)
    x_axis = np.repeat(times,2)  #dm.getMidpoints(times)
    y_axis = np.repeat(hist_times,2)  # hist_times
    plt.plot(x_axis[1:-1], y_axis,lw=linewidth,ls='-',color='red', label='Field halos')
    #plt.plot(x_axis, y_axis,lw=linewidth,ls='-',color='red', label='Field halos')


    plt.xlabel('Infall time [Gyr from z=0]',fontsize=label_font)
    plt.ylabel('Normalized frequency',fontsize=label_font)
    plt.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.xlim((0,14))
    plt.legend(loc='upper left',frameon=False,fontsize=legend_size-2)
    plt.gcf().subplots_adjust(bottom=0.14)
    plt.savefig('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/Paper/Submission2/infall_time_distr.pdf')
    plt.close()



# use this for a jackknife approach to infall time uncertainty
#for i in range(33):
#    infall_distribution_all(True, i)



#infall_distribution_all(False)
#infall_distribution_all(True)

#scale_factors = random_infall_scale(1000)
#plt.hist(scale_factors,bins=np.arange(0,1.04,.04))
#plt.savefig('test_random_sample_scales')
#plt.close()




#0.089960104700241184, 0.0   # why does this turn into a scale factor of 0? it should be 0.27
