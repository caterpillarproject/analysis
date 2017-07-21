import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import statsmodels.nonparametric.smoothers_lowess as smoother
import DwarfMethods as dm
from PlotParams import *
from matplotlib import colors


# z = 9.33, 10.73, 12.1, 13.57
h0 = 0.6711
vmax_ach = 9.48535156 
vmax_filt = 23.54248047 

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    from math import factorial
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients                                                                               
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with                                                                   
    # values taken from the signal itself                                                                   
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')


def get_Barber_data():
    path = '/nfs/blank/h4231/gdooley/DwarfsOfDwarfs/code/for_Dooley/'
    classical = np.loadtxt(path+'classical.dat')
    luminous = np.loadtxt(path+'luminous.dat')
    allsats = np.loadtxt(path+'allSats.dat')
    # convert masses to log10 in Msun, not log10(msun/h)
    classical = np.log10(10**classical/.73)  # hubble value from the Barber paper is 0.73
    luminous = np.log10(10**luminous/.73)
    allsats = np.log10(10**allsats/.73)
    return classical, luminous,allsats
    

# generate log_mass array and a dark fraction array.
# then create an interpolation (linear) between the points
def get_lum_fraction_Barber(plot=False):
    classical, luminous,allsats = get_Barber_data()
    log_mass = np.arange(6.7,10.5,.1) # play with binnings
    hist_all, rvals = np.histogram(allsats,bins=log_mass)
    hist_lum, rvals = np.histogram(luminous,bins=log_mass)
    x_axis = dm.getMidpoints(rvals)
    fracs = np.array(hist_lum,dtype=float)/np.array(hist_all)

    altfracs = np.array(fracs)
    #altfracs[-12] = fracs[-11]
    #altfracs[-11] = 1
    altfracs[-11]=1.00
    altfracs[-12]=0.95
    smoothed_fracs= smoother.lowess(fracs,x_axis,frac=0.25)[:,-1]
    sg_fracs = savitzky_golay(altfracs,9,5)
    tck = interpolate.splrep(x_axis, sg_fracs,k=1)

    n = np.where(tck[1]>1)[0][0]
    tck[1][n:] = 1.0


    if plot:
        plt.plot(x_axis, smoothed_fracs,color='pink',ls='-',lw=2)
        plt.plot(x_axis, sg_fracs,color='green',ls='-',lw=2)
        plt.plot(x_axis, interpolate.splev(x_axis,tck),color='red',ls='--',lw=2)
        plt.plot(x_axis, fracs)
        plt.ylabel('fraction luminous')
        plt.xlabel('$\log_{10}(M_{200}) M_\odot$')
        plt.savefig('lum_frac_barber')   
        plt.close()

    return tck, x_axis, fracs # smoothed_fracs

    """ # for plotting the luminous fraction
    plt.bar(rvals[:-1], hist_all,width=0.1)
    plt.bar(rvals[:-1], hist_lum,width=0.1,color='purple')
    plt.ylabel('Total Number of halos')
    plt.xlabel('$\log_{10}(M_{200}) M_\odot$')
    plt.yscale('log')
    plt.savefig('lum_frac_barber_hist')
    plt.close()
    """

# given halo mass, what fraction are luminous? Halo mass must be in terms of M200 infall.
# For all but moster, must therefore make a conversion from Mhalo to M200 infall.
def frac_luminous(M,tck):  
    if len(M)==0:
        return np.array([])
    values = interpolate.splev(np.log10(M),tck)
    # if >1 or < 0, doesn't matter. always hosts stars or never does when compared to random int
    mask0=values<0
    values[mask0]=0
    mask1=values>1
    values[mask1]=1
    return values

def get_lum_fraction_caterpillar(mdef,z=13,vmax_ach=vmax_ach, vmax_filt=vmax_filt, plot=False):
    if mdef=='m200':
        mstring='infall_mass200'
    if mdef=='mpeak':
        mstring='max_mass'
    if mdef=='m350':
        mstring='max_mass350NFW'
    vstring = 'vmax_'+str(z)
    tstring = 'tvir_'+str(z)
    m200string = 'm200_'+str(z)

    h0=.6711
    # want to collect infall 200 masses for all halos
    # and for just the luminous ones
    m200_all = np.array([])
    m200_lum9 = np.array([])
    hpaths = dm.get_hpaths(False)
    for hpath in hpaths:   #[0:25]:   # I don't remember why I did this
        data = dm.get_extant_data(hpath,False)
        m200_all = np.append(m200_all, np.array(data[mstring]/0.6711,dtype=np.float64))

        # now do the reionization checks for z=11
        mask9 = (data[vstring] > vmax_ach)
        #mask9 = (data[tstring] > 3000)
        #mask9 = (data[m200string]/h0 > 10**6.8)
        #print np.sum(mask), 'num bigger than 15'
        mask2 = data['peak_vmax'] >= vmax_filt  
        #print np.sum(mask&mask2), 'num also peaking higher than 50'
        m200_lum9 = np.append(m200_lum9,  np.array(data[mstring][mask9 | mask2]/0.6711,dtype=np.float64))

    log_mass = np.arange(7.5,10.5,.1) # play with binnings. I don't have any mass smaller than 7.5 recorded
    hist_all, rvals = np.histogram(np.log10(m200_all),bins=log_mass)
    hist_lum9, rvals = np.histogram(np.log10(m200_lum9),bins=log_mass)

    x_axis = dm.getMidpoints(rvals) # this is in log10 form
    fracs9 = np.array(hist_lum9,dtype=float)/np.array(hist_all)

    smoothed_fracs= smoother.lowess(fracs9,x_axis,frac=0.25)[:,-1]
    #smoothed_fracs = savitzky_golay(fracs9,9,5)
    tck = interpolate.splrep(x_axis, smoothed_fracs,k=1)

    if plot:
        tck_barb, xbarber, frac_barber = get_lum_fraction_Barber()
        plt.plot(10**xbarber, interpolate.splev(xbarber,tck_barb),ls='-',lw=linewidth,label='$M_{\mathrm{200}}^{\mathrm{infall}}$ (Barber)')
        plt.plot(10**xbarber, interpolate.splev(xbarber, tck), ls='-',lw=linewidth, label='Caterpillar')
        plt.legend(loc='upper left', fontsize=legend_size, frameon=False)
        plt.ylim((0,1))
        plt.xlim((10**6.9, 10**11))
        plt.xscale('log')
        plt.ylabel('$\mathrm{f_{lum}}$',fontsize=label_font )
        plt.xlabel('$\mathrm{M_{halo}} \, [\mathrm{M_\odot}]$',fontsize=label_font)
        plt.gcf().subplots_adjust(bottom=0.15)
        plt.tick_params(axis='both', which='major', labelsize=tick_size)
        plt.savefig('lum_frac_caterpillar_'+str(z))
        plt.close()

        # for plotting the luminous fraction
        #plt.bar(rvals[:-1], hist_all,width=0.1)
        #plt.bar(rvals[:-1], hist_lum11,width=0.1,color='purple')
        #plt.ylabel('Total Number of halos')
        #plt.xlabel('$\log_{10}(M_{200}) M_\odot$')
        #plt.yscale('log')
        #plt.savefig('lum_frac_caterpillar_hist')
        #plt.close()
    return tck, x_axis, fracs9

def plot_frac_lum_mdef():
    tck_barb, xbarber, frac_barber = get_lum_fraction_Barber()
    tck_m350, xm350, frac_m350 = get_lum_fraction_m350()
    tck_mpeak, xmpeak, frac_mpeak = get_lum_fraction_mpeak()
    plt.plot(10**xbarber, interpolate.splev(xbarber,tck_barb),ls='-',lw=linewidth,label='$M_{\mathrm{200}}^{\mathrm{infall}}$ (Barber)', color=colors.hex2color('#377eb8'))
    plt.plot(10**xm350, interpolate.splev(xm350,tck_m350),ls='-',lw=linewidth,label='$M_{\mathrm{350}}^{\mathrm{peak}}$', color=colors.hex2color('#e41a1c'))
    plt.plot(10**xmpeak, interpolate.splev(xmpeak,tck_mpeak),ls='-',lw=linewidth,label='$M_{\mathrm{vir}}^{\mathrm{peak}}$', color=colors.hex2color('#ff7f00'))
    plt.legend(loc='upper left', fontsize=legend_size, frameon=False)

    plt.ylim((0,1))
    plt.xlim((10**6.9, 10**11))
    plt.xscale('log')
    plt.ylabel('$\mathrm{f_{lum}}$',fontsize=label_font )
    plt.xlabel('$\mathrm{M_{halo}} \, [\mathrm{M_\odot}]$',fontsize=label_font)
    plt.gcf().subplots_adjust(bottom=0.15)

    plt.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.savefig('lum_frac_mdef.pdf')
    plt.close()




def plot_frac_lum_talk():
    # get_lum_fraction_mpeak() 
    tck_cat13, xcat13, fracs13 = get_lum_fraction_caterpillar(mdef='mpeak',z=13)
    tck_cat10, xcat10, fracs10 = get_lum_fraction_caterpillar(mdef='mpeak',z=10)
    tck_cat7, xcat7, fracs7 = get_lum_fraction_caterpillar(mdef='mpeak',z=7)
    xcat = np.arange(7.0,10.5,.1)
    
    plt.plot(10**xcat, interpolate.splev(xcat,tck_cat13),ls='-',lw=5,label='z=13.3',color=colors.hex2color('#377eb8'))
    plt.plot(10**xcat, interpolate.splev(xcat,tck_cat10),ls='-',lw=5,label='z=10.3', color=colors.hex2color('#e41a1c'))
    plt.plot(10**xcat, interpolate.splev(xcat,tck_cat7),ls='-',lw=5,label='z=7.3', color=colors.hex2color('#ff7f00'))
    plt.legend(loc='lower right', fontsize=25, frameon=False)

    plt.ylim((0,1))
    plt.xlim((10**6.9, 10**11))
    plt.xscale('log')
    plt.ylabel('$\mathrm{f_{lum}}$',fontsize=25)
    plt.xlabel('$\mathrm{M_{halo}} \, [\mathrm{M_\odot}]$',fontsize=25)
    plt.gcf().subplots_adjust(bottom=0.15,left=.13,top=0.95,right=0.93)
    print 'changed'
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.savefig('lum_frac_TALK')
    plt.close()




def plot_frac_luminous():
    tck_barb, xbarber, frac_barber = get_lum_fraction_Barber()
    tck_cat, xcat, fracs9,fracs11,fracs12,fracs13 = get_lum_fraction_caterpillar(mdef='m200')
    tck_cat_350, xcat_350, fracs9_350,fracs11_350,fracs12_350,fracs13_350 = get_lum_fraction_caterpillar(mdef='m350')
    tck_cat_peak, xcat_peak, peak_fracs9,peak_fracs11,peak_fracs12,peak_fracs13 = get_lum_fraction_caterpillar(mdef='mpeak')

    plt.plot(xbarber, interpolate.splev(xbarber,tck_barb),color='red',ls='--',lw=3)
    plt.plot(xbarber, frac_barber, label='Barber',lw=3)

    plt.plot(xcat, fracs13,lw=3,label='m200')
    plt.plot(xcat_350, fracs13_350,lw=3,label='m350')
    plt.plot(xcat_peak, peak_fracs13,lw=3,label='mpeak')

    #plt.plot(xcat, fracs9,lw=3,label='$z_{reion} = 9.33$')
    #plt.plot(xcat, fracs11,lw=3,label='$z_{reion} = 10.73$')
    #plt.plot(xcat, fracs12,lw=3,label='$z_{reion} = 12.1$')
    #plt.plot(xcat, fracs13,lw=3,label='$z_{reion} = 13.57$')
    plt.legend(loc='upper left', fontsize=legend_size, frameon=False)

    plt.ylabel('Fraction Luminous', fontsize=label_font)
    plt.xlabel('$\log_{10}(M_{\mathrm{halo}}) M_\odot$',fontsize=label_font)

    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig('lum_frac_comparison')
    plt.close()




# can take lum_frac_barber applied to all the m200 masses and generate random yes or no on them.
# get data on yes or no of the m350 values and plot the fractions.
# re-generate it 100 times and combine the values.
def barber_to_other():
    tck_barb, xbarber, frac_barber = get_lum_fraction_Barber()    

    lum_m350=[]; all_m350 = []
    lum_mpeak = []; all_mpeak = []
    hpaths = dm.get_hpaths(False)

    for hpath in hpaths:
        data = dm.get_extant_data(hpath,False)
        m200 = data['infall_mass200']/h0
        lum_chances = frac_luminous(m200, tck_barb)

        for i in range(20):
            randnums = np.random.rand(len(lum_chances))
            mask = randnums < lum_chances # chooses which halos will form stars
        
            lum_m350.append( np.array(data['max_mass350NFW'][mask]/h0))
            all_m350.append( np.array(data['max_mass350NFW']/h0))
            
            lum_mpeak.append( np.array(data['max_mass'][mask]/h0))
            all_mpeak.append( np.array(data['max_mass']/h0))

    lum_m350 = np.array([item for arr in lum_m350 for item in arr])
    all_m350 = np.array([item for arr in all_m350 for item in arr])
    lum_mpeak = np.array([item for arr in lum_mpeak for item in arr])
    all_mpeak = np.array([item for arr in all_mpeak for item in arr])

    np.save('lum_m350', lum_m350)
    np.save('all_m350', all_m350)
    np.save('lum_mpeak', lum_mpeak)
    np.save('all_mpeak', all_mpeak)



def get_m350_data():
    path = '/nfs/blank/h4231/gdooley/DwarfsOfDwarfs/code/'
    return np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/lum_m350.npy'), np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/dark_m350.npy')
    

# generate log_mass array and a dark fraction array.
# then create an interpolation (linear) between the points
def get_lum_fraction_m350():
    luminous, allsats = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/lum_m350.npy'), np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/all_m350.npy')
    allsats = allsats[allsats>0]
    luminous = luminous[luminous>0]  # had to add these masks in Jan 2017, do not know why I didn't need them before. Update to numpy histogram??

    log_mass = np.arange(6.7,10.5,.1) # play with binnings
    hist_all, rvals = np.histogram(np.log10(allsats),bins=log_mass)
    hist_lum, rvals = np.histogram(np.log10(luminous),bins=log_mass)
    x_axis = dm.getMidpoints(rvals)
    fracs = np.array(hist_lum,dtype=float)/np.array(hist_all)

    smoothed_fracs= smoother.lowess(fracs,x_axis,frac=0.25)
    tck = interpolate.splrep(smoothed_fracs[:,0], smoothed_fracs[:,-1],k=1)
    return tck, x_axis, fracs # smoothed_fracs


def get_lum_fraction_mpeak():
    luminous, allsats = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/lum_mpeak.npy'), np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/all_mpeak.npy')
    log_mass = np.arange(7.5,10.5,.1) # play with binnings
    hist_all, rvals = np.histogram(np.log10(allsats),bins=log_mass)
    hist_lum, rvals = np.histogram(np.log10(luminous),bins=log_mass)
    x_axis = dm.getMidpoints(rvals)
    fracs = np.array(hist_lum,dtype=float)/np.array(hist_all)

    smoothed_fracs= smoother.lowess(fracs,x_axis,frac=0.25)
    tck = interpolate.splrep(smoothed_fracs[:,0], smoothed_fracs[:,-1],k=1)
    return tck, x_axis, fracs # smoothed_fracs





#def frac_luminous(M,tck):  
#    if len(M)==0:
#        return np.array([])
#    return interpolate.splev(np.log10(M),tck)
    # if >1 or < 0, doesn't matter. always hosts stars or never does when compared to random int
    #mask0=values<0
    #values[mask0]=0
    #mask1=values>1
    #values[mask1]=1
    #return values


#def get_lum_frac_error(vmax_ach=9.91434884, vmax_filt= 22.73337007 , plot=True, z=13):  # for the lowess smoother, 0:25
#def get_lum_frac_error(vmax_ach= 10.20675659, vmax_filt=  22.16201782 , plot=True, z=13): # sav gol, 0:25
def get_lum_frac_error(vmax_ach=9.48535156 , vmax_filt=  23.54248047  , plot=False, z=8): # lowess, all
    vstring = 'vmax_'+str(z)
    mstring='infall_mass200'
    h0=.6711
    # want to collect infall 200 masses for all halos
    # and for just the luminous ones
    m200_all = np.array([])
    m200_lum13 = np.array([])
    hpaths = dm.get_hpaths(False)
    for hpath in hpaths:   #[0:25]   # something weird with number 30
        data = dm.get_extant_data(hpath,False)
        m200_all = np.append(m200_all, np.array(data[mstring]/0.6711,dtype=np.float64))

        # now do the reionization checks
        mask13 = (data[vstring] > vmax_ach)
        #mask13 = (data['tvir_13'] > 3000)
        #mask13 = (data['m200_13']/h0 > 10**6.8)
        mask2 = data['peak_vmax'] >= vmax_filt #25
        m200_lum13 = np.append(m200_lum13,  np.array(data[mstring][mask13 | mask2]/0.6711,dtype=np.float64))

    log_mass = np.arange(7.5,10.5,.1) # play with binnings
    hist_all, rvals = np.histogram(np.log10(m200_all),bins=log_mass)
    hist_lum13, rvals = np.histogram(np.log10(m200_lum13),bins=log_mass)

    x_axis = dm.getMidpoints(rvals) # this is in log10 form
    yvals = np.array(hist_lum13,dtype=float)/np.array(hist_all)

    smoothed_fracs= smoother.lowess(yvals,x_axis,frac=0.25)[:,-1]
    #smoothed_fracs = savitzky_golay(yvals,9,5)
    tck = interpolate.splrep(x_axis, smoothed_fracs,k=1)

    tck_barb, xbarber, frac_barber = get_lum_fraction_Barber()
    if plot:
        plt.plot(10**xbarber, interpolate.splev(xbarber,tck_barb),ls='-',lw=linewidth,label='$M_{\mathrm{200}}^{\mathrm{infall}}$ (Barber)')
        #print xbarber, 'xbarber'
        plt.plot(10**xbarber, interpolate.splev(xbarber, tck), ls='-',lw=linewidth, label='Caterpillar')
        plt.legend(loc='upper left', fontsize=legend_size, frameon=False)
        plt.ylim((0,1))
        plt.xlim((10**6.9, 10**11))
        plt.xscale('log')
        plt.ylabel('$\mathrm{f_{lum}}$',fontsize=label_font )
        plt.xlabel('$\mathrm{M_{halo}} \, [\mathrm{M_\odot}]$',fontsize=label_font)
        plt.gcf().subplots_adjust(bottom=0.15)
        plt.tick_params(axis='both', which='major', labelsize=tick_size)
        plt.savefig('lum_frac_caterpillar')
        plt.close()

    # compute the error
    error = np.sum((interpolate.splev(xbarber,tck_barb) - interpolate.splev(xbarber, tck))**2)
    #print error
    return error

def helper(vmaxes):
    return get_lum_frac_error(vmaxes[0], vmaxes[1], plot=False)







# plot diversity of models
def add_lum_frac_axis(ax, vmax_ach, vmax_filt, z, label): 
    vstring = 'vmax_'+str(z)
    mstring='infall_mass200'
    h0=.6711
    m200_all = np.array([])
    m200_lum13 = np.array([])
    hpaths = dm.get_hpaths(False)
    for hpath in hpaths:   #[0:25]   # something weird with number 30
        data = dm.get_extant_data(hpath,False)
        m200_all = np.append(m200_all, np.array(data[mstring]/0.6711,dtype=np.float64))

        # now do the reionization checks for z=11
        mask13 = (data[vstring] > vmax_ach)
        mask2 = data['peak_vmax'] >= vmax_filt
        m200_lum13 = np.append(m200_lum13,  np.array(data[mstring][mask13 | mask2]/0.6711,dtype=np.float64))
    log_mass = np.arange(7.5,10.5,.1) # play with binnings
    hist_all, rvals = np.histogram(np.log10(m200_all),bins=log_mass)
    hist_lum13, rvals = np.histogram(np.log10(m200_lum13),bins=log_mass)

    x_axis = dm.getMidpoints(rvals) # this is in log10 form
    yvals = np.array(hist_lum13,dtype=float)/np.array(hist_all)

    smoothed_fracs= smoother.lowess(yvals,x_axis,frac=0.25)[:,-1]
    #smoothed_fracs = savitzky_golay(yvals,9,5)
    tck = interpolate.splrep(x_axis, smoothed_fracs,k=1)
    ax.plot(10**x_axis, interpolate.splev(x_axis, tck), ls='-',lw=linewidth, label=label)


def plot_rion_vary():
    fig, ax = plt.subplots()
    #tck_barb, xbarber, frac_barber = get_lum_fraction_Barber()
    #ax.plot(10**xbarber, interpolate.splev(xbarber,tck_barb),ls='-',lw=linewidth,label='$M_{\mathrm{200}}^{\mathrm{infall}}$ (Barber)')
    add_lum_frac_axis(ax, vmax_ach, vmax_filt, 13, label='z=13.3')     
    add_lum_frac_axis(ax, vmax_ach, vmax_filt, 14, label='z=14.3')     
    add_lum_frac_axis(ax, vmax_ach, vmax_filt, 11, label='z=11.3')     
    add_lum_frac_axis(ax, vmax_ach, vmax_filt, 9, label='z=9.3')

    add_lum_frac_axis(ax, vmax_ach*.75, vmax_filt*.75, 13, label='vmax lower')
    add_lum_frac_axis(ax, vmax_ach*1.25, vmax_filt*1.25, 13, label='vmax bigger')

    #add_lum_frac_axis(ax, 5, 20, 13, label='$v_{\mathrm{max}}^{\mathrm{filt}} = 20 \, \mathrm{km/s}$, $v_{\mathrm{max}}^{\mathrm{pre}} = 5 \, \mathrm{km/s}$  ')
    #add_lum_frac_axis(ax, 15, 30, 13, label=' $v_{\mathrm{max}}^{\mathrm{filt}} = 30 \, \mathrm{km/s}$, $v_{\mathrm{max}}^{\mathrm{pre}} = 15 \, \mathrm{km/s}$  ')


    plt.legend(loc='upper left', fontsize=legend_size-1, frameon=False)
    plt.ylim((0,1))
    plt.xlim((10**7, 10**10))
    plt.xscale('log')
    plt.ylabel('$\mathrm{f_{lum}}$',fontsize=label_font )
    plt.xlabel('$\mathrm{M_{halo}} \, [\mathrm{M_\odot}]$',fontsize=label_font)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.savefig('lum_frac_vary')
    plt.close()





#get_lum_fraction_Barber(plot=True)



#print get_lum_frac_error(plot=True)


# Code to get best fit vmaxes:
"""
import scipy.optimize
r = scipy.optimize.fmin(helper, x0=np.array([10,25]), maxiter=100)
print r

"""

##### WITH THE OLD Z=13 and OLD LOWESS FIT
# z=13 better solutions
#9.90332031,  22.89550781   0.022858
# z=12
#10.85888672,  23.77319336   0.037238
#z=11
# 13.046875  ,  23.05664062   0.032135
###### END OLD FITS ##########

##### NEW FITS  ######
# z=13 better solutions
# 10.20675659  22.16201782     0.035851   # with savitzky golay
# 9.91434884  22.73337007   0.018449    # with lowess

# 9.48535156  23.54248047    0.031688   # with lowess and all halos
###### END NEW FITS  ######



## Figure 2 #####
#plot_frac_lum_mdef():


