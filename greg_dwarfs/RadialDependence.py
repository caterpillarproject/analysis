import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import haloutils as htils
from scipy import stats
from PlotParams import *
from scipy.optimize import curve_fit
from scipy import interpolate
import DwarfMethods as dm
from scipy.integrate import dblquad
import field_shmf as fshmf


def plot_radial_fits():
    radius2 = np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2])
    radius = np.array([0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2])

    Khost = np.array([0.0, 0.000008, 0.000047,0.000129,0.00023,0.000362,0.000504,0.000655,0.000811,0.000971,0.001122,0.001269,0.001421,0.001549,0.001681,0.001802,0.001905,0.002003,0.002095,0.00218,0.002257])

    Kfield = np.array([0.0, 0.000008, 0.00006,0.000144,0.00025,0.000386,0.000511,0.00078,0.000961,0.001139,0.001314,0.001526,0.001677,0.001805,0.001972,0.002118,0.002277,0.002452,0.002603,0.002761,0.003])

    Khost2 = np.array([0, 0.000008, 0.000045, 0.000139, 0.000245, 0.000384, 0.000521,0.000672,0.000815,0.000992,0.001133,0.001282,0.001427,0.001545,0.001662,0.001784,0.001879,0.001962,0.002068,0.002148,0.002228])  # no rmax cut

    volume = 4/3.*np.pi*radius**3

    # fit a line to the data
    slope, intercept, r_value, p_value, std_err = stats.linregress(radius,Khost)
    slopeV, interceptV, r_valueV, p_valueV, std_errV = stats.linregress(np.log10(volume),np.log10(Khost))
    slopeF, interceptF, r_valueF, p_valueF, std_errF = stats.linregress(radius,Kfield)
    slopeFV, interceptFV, r_valueFV, p_valueFV, std_errFV = stats.linregress(np.log10(volume),np.log10(Kfield))
    #print slope, intercept, r_value, p_value, std_err, 'host slope, intercept, r, p, std error'
    #print slopeF, interceptF, r_valueF, p_valueF, std_errF, 'field slope, intercept, r, p, std error'
    #print slopeV, interceptV, r_valueV, p_valueV, std_errV, 'host by volume slope, intercept, r, p, std error'

    #slope, a,b,c = np.linalg.lstsq(np.vstack([radius, np.ones(len(radius))]).T, Kfield)
    slopeF2, _,_,_ = np.linalg.lstsq(radius[:,np.newaxis],Kfield)

    (c1,c2,c3,c4) = np.polyfit(radius,Khost,3)
    (d1,d2,d3,d4) = np.polyfit(radius2,Khost2,3)

    # try a logistics fit
    def logistic(x, a,b,k):
        return a/(1+ b* np.e**(-k * (x)))
    # end line fitting

    popt, pcov = curve_fit(logistic, radius2[1:], Khost2[1:])
    print logistic(1,popt[0],popt[1],popt[2]), 'logistic at 1. should be .0011' 


    #plt.plot(radius, interceptF+slopeF*radius,color='red',ls='--')
    #plt.plot(radius, slopeF2*radius,color='red',ls='--')
    
    norm = c4+c3+c2+c1
    print norm, 'should be close to K0 = 0.001133'
    #print (c4+c3+c2+c1)/norm, 'value at rvir is 1 rvir'
    #print c4,c3,c2,c1, 'c4,c3,c2,c1'
    print c4/norm, c3/norm, c2/norm, c1/norm, 'normed constants'


    norm = d4+d3+d2+d1
    #print (d4+d3+d2+d1)/norm, 'sum of d values. value at rvir is 1 rvir'
    #print d4,d3,d2,d1, 'd4,d3,d2,d1'
    print d4/norm, d3/norm, d2/norm, d1/norm, 'normed constants with no rmax cut'


    rr = np.arange(0.05,2,.01)
    plt.plot(rr, logistic(rr, popt[0],popt[1],popt[2]), color='magenta',ls='--',label='logistic fit')
    plt.plot(rr,d4+d3*rr+d2*rr**2+d1*rr**3,color='cyan',ls='--')
    plt.plot(rr,c4+c3*rr+c2*rr**2+c1*rr**3,color='green',ls='--')
    #plt.plot(radius,-.00015+.00071*radius+.00086*radius**2-.00031*radius**3,color='purple',ls='--')
    plt.plot(radius,Khost,color='blue',label='host')
    plt.plot(radius,Kfield,color='orange',label='field')
    plt.plot(radius2,Khost2,color='red',label='host no cut')
    plt.legend(loc='upper left')
    plt.xlabel('radius [kpc]')
    plt.xlim((0.1,2))
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig('radialfit_logscale')
    
    plt.close()


    rr = np.arange(0.0,2,.01)
    plt.plot(rr, logistic(rr, popt[0],popt[1],popt[2]), color='magenta',ls='--',label='logistic fit')
    plt.plot(rr,d4+d3*rr+d2*rr**2+d1*rr**3,color='cyan',ls='--')
    plt.plot(rr,c4+c3*rr+c2*rr**2+c1*rr**3,color='green',ls='--')
    print c4+c3*.1+c2*.1**2+c1*.1**3, 'value of green fit at 0.1'
    #plt.plot(radius,-.00015+.00071*radius+.00086*radius**2-.00031*radius**3,color='purple',ls='--')
    plt.plot(radius,Khost,color='blue',label='host')
    plt.plot(radius,Kfield,color='orange',label='field')
    plt.plot(radius2,Khost2,color='red',label='host no cut')
    plt.legend(loc='upper left')
    plt.xlabel('radius [kpc]')
    plt.xlim((0.0,2))
    plt.savefig('radialfit')
    plt.close()
    
    """
    plt.plot(volume,Khost,color='blue',label='host')
    plt.plot(volume,Kfield,color='orange',label='field')
    plt.xlabel('volume $[kpc^3]$')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper left')
    plt.savefig('radialfit_volume')
    plt.xlim((0,2))
    plt.close()
    """

def get_field_halos(cat,hpath,hostID,mlow=10,mhigh=11.5):
    # get halos beyond host halo virial radius, but less than 
    # contamination radius
    hosts = cat.get_hosts()
    dists = dm.distance(cat.ix[hostID][['posX','posY','posZ']], hosts[['posX','posY','posZ']])
    contam_dist = dm.get_contam_dist(hpath)
    mask = dists < contam_dist
    mass_mask = (10**mlow < np.array(hosts['mgrav']/cat.h0)) &(10**mhigh > np.array(hosts['mgrav']/cat.h0))
    return hosts[mask & mass_mask]


# values correspond to host halo values and the commented normed values k1=-.1366
def getK_not_normed(r):
    return -0.00015148503612 + 0.000708840718105*r + 0.000861092068406*r**2 - 0.00030949197861*r**3
    #return -.00015+.00071*r+.00086*r**2-.00031*r**3

# input the radial ratio r/rvir. normalized to 1 at r=1
def getK(r):
    return -0.0439735566925 + 0.39131605932*r + 0.996495770479*r**2 - 0.343838273107*r**3  # fit to rmaxcut = False. best so far

# need a piecewise function that fits the curve better at low r.
# maybe linear is ok, but need to define at what r to make the transition.


def getK_protected(R):
    discon = 0.15
    if type(R)==float:
        if r < discon:
            m = getK(discon)
            return m/discon * R
        else:
            return getK(R)
    else:
        # assume input is an array
        result = getK(R)
        mask = R < discon
        m = getK(discon)
        y = m/discon * R[mask]
        result[mask] = y
        return result


# normed values of k
#k1 = -0.136601512859; k2=0.639196563237; k3=0.77648909931; k4= -0.279084149688  # fit to rmaxcut=True above 2.0
k1 = -0.0439735566925; k2=0.39131605932; k3=0.996495770479; k4=-0.343838273107  # fit to rmaxcut = False. best so far
#k1=-0.0362134075824; k2=0.297178749054; k3=1.11253297655; k4=-0.37349831802  # fit to rmaxcut = True

def inversetan(x, const1, const2, const3):
    return const2 * np.arctan(x/const1 - const3)

def analytic_radial(R):
    d4 = 0.000212585034013
    d3 = -0.261479591837
    d2 = 6.88775510204
    d1 = -7.03463203463
    
    const1 = 0.529767613564
    intercept = 0.108858675195

    const1 = 0.529767613564 
    const2 = 0.966702936991
    const3 = 0.205483824242

    if type(R)==int or type(R)==float:
        if R < 0.2:
            return d4+d3*R+d2*R**2+d1*R**3
        else:
            return inversetan(R, const1, const2, const3)
    else:
        masklow = R<0.2
        maskhigh = R>=0.2
        firstpart = d4+d3*R[masklow]+d2*R[masklow]**2+d1*R[masklow]**3
        secondpart = inversetan(R[maskhigh], const1, const2, const3)
        return np.append(firstpart, secondpart)




def get_normed_factor(r):
    # constatns are a,b,c,d in the order written to correspond tomy math notes
    return k1 + k2*r + k3*r**2 + k4*r**3


def plot_normed_fit():
    radius = np.arange(0.2,2.1,.1)
    plt.plot(radius, get_normed_factor(radius))
    plt.savefig('normed_radial_fit')
    plt.close()


def get_normed_los(R,Z):
    g = np.arctan2(R,Z)
    term1 = Z*k2/2. * np.log(1+R**2/Z**2)
    term2 = k3*Z**2 * (np.sqrt(1+R**2/Z**2) - 1)
    term3 = k4*R**2*Z/2.

    term4 = k2*R*(np.pi/2 - g)
    term5 = k3*R**2 * np.log( Z/R *(1+np.sqrt(1+R**2/Z**2)))
    term6 = k4*R**2*Z

    return term1+term2+term3+term4+term5+term6 + k1  

def get_normed_los_protected(R,Z):
    if type(R)==float:
        if R < 0.1:
            m = get_normed_los(0.1, Z)
            return m/.1 * R
        else:
            return get_normed_los(R, Z)
    else:
        # assume input is an array
        result = get_normed_los(R,Z)
        mask = R < 0.1
        m = get_normed_los(0.1, Z)
        y = m/.1 * R[mask]
        result[mask] = y
        return result
        


# could add numerical integration to the plot to show how well it fits
# 
def get_numerical_normed_los(R,Z, tck):
    def dkdr(r):
        return interpolate.splev(r,tck,der=1) / 0.001133

    t0 = np.arctan2(R,Z)
    term1 = dblquad(lambda r, theta: dkdr(r)*np.sin(theta), 0, t0, lambda x: 0, lambda x: Z/np.cos(x))
    term2 = dblquad(lambda r, theta: dkdr(r)*np.sin(theta), t0, np.pi/2, lambda x: 0, lambda x: R/np.sin(x))
    print 'finished an iteration of integrating', R
    return term1[0] + term2[0]

def get_numerical_normed_los_better(R,Z):
    rrr = np.arange(0,2.001,.05)
    normed_n_within_lum = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/SHMF/normed_radial_survivors.npy')
    normed_n_within_lum[20:] = .4*rrr[20:]+.6
    tck =  interpolate.splrep(rrr,normed_n_within_lum, k=1)

    def dkdr(r):
        return interpolate.splev(r,tck,der=1)

    t0 = np.arctan2(R,Z)
    term1 = dblquad(lambda r, theta: dkdr(r)*np.sin(theta), 0, t0, lambda x: 0, lambda x: Z/np.cos(x))
    term2 = dblquad(lambda r, theta: dkdr(r)*np.sin(theta), t0, np.pi/2, lambda x: 0, lambda x: R/np.sin(x))
    print 'finished an iteration of integrating', R
    print term1[0] + term2[0], 'is K(R) for', R
    return term1[0] + term2[0]


def get_normed_los_better(R,Z):
    if Z == 1.5:
        ratio = np.linspace(0,1.5,31)
        K = np.array([0.0, 0.02749006,  0.09507474,  0.187685, 0.28999924, 0.3882824, 0.48242521, 0.56322701,0.63678851, 0.70413224, 0.76000307, 0.81019972, 0.85594824, 0.89798589, 0.9362271, 0.97067692, 1.00329909, 1.03375505, 1.06024407, 1.08245111, 1.10328478, 1.12271311, 1.14153286, 1.16037277, 1.17838655, 1.1962006, 1.2135096, 1.23041903, 1.24710864, 1.26325687, 1.27909509]) 
        tck =  interpolate.splrep(ratio, K, k=1)
        return interpolate.splev(R,tck)
    else:
        print 'not computed yet for this Z value'
        return 0

        
"""
ratio = np.linspace(0.05,1.5,30)
raw_values = np.array([get_numerical_normed_los_better(R,1.5) for R in ratio])
print raw_values, 'K(R) with Z=1.5. save these values for later'
print ratio
"""

# this is for khost2
#[ 0.01189194  0.03825338  0.1338373   0.40731802  0.69055662  0.94740314  1.16469735 
# 1.33790466  1.4688462   1.52461371] values 3 for [0.05,.1,.2,.4,.6,.8,1,1.2,1.4,1.5]

# this is for khost
#[ 0.01172879  0.03766099  0.12990579  0.39270223  0.67768313  0.93840293
#  1.15882717  1.34039798  1.48181768  1.53956582] values 3 for [0.05,.1,.2,.4,.6,.8,1,1.2,1.4,1.5]

# this is z=0.5 for khost
#[ 0.01007687  0.03110473  0.10446132  0.30095863  0.49895946  0.66964709
#  0.80700626  0.91751962  1.00130305  1.0353223 ] values 4 for [0.05,.1,.2,.4,.6,.8,1,1.2,1.4,1.5]

# this is z=0.5 for khost2
#[ 0.01027001  0.03181417  0.10881709  0.31673135  0.51349503  0.68128896
#  0.81923617  0.92638784  1.00304656  1.03662456] values 4 for [0.05,.1,.2,.4,.6,.8,1,1.2,1.4,1.5]


def plot_los_frac():
    radius = np.array([0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2])
    Khost2 = np.array([0, 0.000008, 0.000045, 0.000139, 0.000245, 0.000384, 0.000521,0.000672,0.000815,0.000992,0.001133,0.001282,0.001427,0.001545,0.001662,0.001784,0.001879,0.001962,0.002068,0.002148,0.002228])
    Khost = np.array([0.0, 0.000008, 0.000047,0.000129,0.00023,0.000362,0.000504,0.000655,0.000811,0.000971,0.001122,0.001269,0.001421,0.001549,0.001681,0.001802,0.001905,0.002003,0.002095,0.00218,0.002257])

    # numerical integration
    tck =  interpolate.splrep(radius, Khost2, k=1)

    Z1 = 1; Z2 = 1.5; Z3 = 0.5
    r = np.arange(0.1,1.91,.05)
    values = get_normed_los(r,Z1)
    values2 = get_normed_los(r,Z2)
    #values3 = get_normed_los(r,Z3)
    r3 = [0.05,.1,.2,.4,.6,.8,1,1.2,1.4,1.5]
    raw_values = [0.01189194, 0.03825338, 0.1338373, 0.40731802, 0.69055662, 0.94740314, 1.16469735, 1.33790466, 1.4688462, 1.52461371]
    raw_values2 = [0.01238067, 0.04019938, 0.1415055, 0.43624872, 0.74980424, 1.04072008, 1.29361954, 1.50023667, 1.66275806, 1.73320591] # this is for Z = 1.5
    #raw_values3 = [ 0.01027001, 0.03181417, 0.10881709, 0.31673135, 0.51349503, 0.68128896, 0.81923617, 0.92638784, 1.00304656, 1.03662456] # this is for z = 0.5

    #raw_values2 = np.array([get_numerical_normed_los(R,Z2,tck) for R in r3])
    #vvv = np.array([get_numerical_normed_los(R,0.5,tck) for R in r3])
    #print vvv, 'values for r3 Z=0.5'
    #print values4, 'values 4 for [0.05,.1,.2,.4,.6,.8,1,1.2,1.4,1.5]'
    # copy values 4 to not have to do again

    plt.plot(r3,raw_values2, lw=linewidth, label='$Z=%.1f$ (actual)' %Z2,color='orange',ls='--')
    plt.plot(r3,raw_values, lw=linewidth, label='$Z=%.1f$ (actual)' %Z1,color='red',ls='--')
    #plt.plot(r3,raw_values3, lw=linewidth, label='$Z=%.1f$ raw data' %Z3,color='saddlebrown',ls='--')


    plt.plot(r,values2, lw=linewidth, label='$Z=%.1f$ (analytic)' %Z2,color='green')
    plt.plot(r,values, lw=linewidth, label='$Z=%.1f$ (analytic)' %Z1,color='blue')
    #plt.plot(r,values3, lw=linewidth, label='$Z=%.1f$' %Z3,color='darkturquoise')


    #Z3=2
    #values2 = get_normed_los(r,Z3)
    #plt.plot(r,values2, lw=linewidth,label='$Z=%.1f$' %Z3,color='green')

    plt.ylabel('K(R)',fontsize=label_font)
    plt.xlabel('$R \equiv R_{\mathrm{fov}}/ R_{\mathrm{vir}}$', fontsize=label_font)
    plt.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.xlim((0,1.5))
    plt.ylim((0,1.8))
    plt.grid(alpha=0.4)
    plt.legend(frameon=False,loc='upper left',fontsize=legend_size)
    plt.savefig('FieldOfView_MultiplicativeFactor.pdf')
    plt.close()


def plot_los_frac_revised():
    plt.figure()
    f, ax = plt.subplots()
    Z1 = 1; Z2 = 1.5
    r = np.arange(0.1,1.91,.05)
    values = get_normed_los(r,Z1)
    values2 = get_normed_los(r,Z2)

    #r3 = [0.0, 0.05,.1,.2,.4,.6,.8,1,1.2,1.4,1.5]
    #raw_valuesZ1 = [0.0, 0.0291206531785, 0.102944601771, 0.300524769518, 0.639095896768, 0.845651561252, 0.973444262548, 1.05279195574, 1.11186828599,  1.16443114853, 1.18852089836]
    #raw_valuesZ2 = [0.0, 0.02928724, 0.10360976, 0.30315803, 0.64919431, 0.86710225, 1.00873837, 1.103422, 1.17847986, 1.24713824, 1.27914254]


    r3 = [0.0, 0.05,.1,.2, .3,.4,.5,.6,.7,.8,.9, 1, 1.1, 1.2, 1.3, 1.4,1.5]
    raw_valuesZ1 = [0.0, 0.0291206531785, 0.102944601771, 0.300524769518, 0.48968238, 0.639095896768, 0.75530574, 0.845651561252, 0.91827683, 0.973444262548, 1.01803111, 1.05279195574, 1.08319479, 1.11186828599, 1.13880562, 1.16443114853, 1.18852089836]
    raw_valuesZ2 = [0.0, 0.02928724, 0.10360976, 0.30315803, 0.495490586338, 0.64919431, 0.770664268383, 0.86710225, 0.946403823717, 1.00873837, 1.06091190553, 1.103422, 1.14177189992, 1.17847986, 1.21349950337, 1.24713824, 1.27914254]


    ax.plot(r3,raw_valuesZ2, lw=linewidth, label='$Z=%.1f$' %Z2,color='green',ls='-')
    ax.plot(r3,raw_valuesZ1, lw=linewidth, label='$Z=%.1f$' %Z1,color='blue',ls='-')
    print 'done with first plots'
    
    #r3 = [.3,.5,.7,.9,1.1,1.3]
    #raw_valuesZ1 = np.array([get_numerical_normed_los_better(R,Z1) for R in r3])
    #print raw_valuesZ1  # [ 0.48968238  0.75530574  0.91827683  1.01803111  1.08319479  1.13880562]
    #raw_valuesZ2 = np.array([get_numerical_normed_los_better(R,Z2) for R in r3])
    #print raw_valuesZ2 # [ 0.49549059  0.77066427  0.94640382  1.06091191  1.1417719   1.2134995 ]
    

    #plt.plot(r,values2, lw=linewidth, label='$Z=%.1f$ (analytic)' %Z2,color='green')
    #plt.plot(r,values, lw=linewidth, label='$Z=%.1f$ (analytic)' %Z1,color='blue')
    #print 'done with 2nd plots'

    ax.set_ylabel('$\mathrm{K_{los}(R)}$',fontsize=label_font)
    ax.set_xlabel('$R \equiv R_{\mathrm{fov}}/ R_{\mathrm{vir}}$', fontsize=label_font)
    
    yvals = ((0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1,1.1,1.2,1.3,1.4))
    ax.set_yticks(yvals)
    ax.set_yticklabels([ '0.0', '','0.2', '', '0.4','','0.6','','0.8','','1.0','','1.2','','1.4'] )

    xvals = ((0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1,1.1,1.2,1.3,1.4,1.5))
    ax.set_xticks(xvals)
    ax.set_xticklabels([ '0.0', '','0.2', '', '0.4','','0.6','','0.8','','1.0','','1.2','','1.4',''] )


    plt.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.xlim((0,1.5))
    plt.ylim((0,1.4))
    plt.grid(alpha=0.4)
    plt.legend(frameon=False,loc='upper left',fontsize=legend_size)
    plt.savefig('FieldOfView_MultiplicativeFactor_revised.pdf')
    plt.close()





def plot_radial_fits_careful():
    radius =  np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2])
    Khost2 = np.array([0, 0.000008, 0.000045, 0.000139, 0.000245, 0.000384, 0.000521,0.000672,0.000815,0.000992,0.001133,0.001282,0.001427,0.001545,0.001662,0.001784,0.001879,0.001962,0.002068,0.002148,0.002228])  # no rmax cut
    volume = 4/3.*np.pi*radius**3

    # fit a polynomial to the data
    (d1,d2,d3,d4) = np.polyfit(radius,Khost2,3)

    # try a logistics fit
    def logistic(x, a,b,k):
        return a/(1+ b* np.e**(-k * (x)))
    # end line fitting
    popt, pcov = curve_fit(logistic, radius[1:], Khost2[1:])
    logistic_norm = logistic(1,popt[0],popt[1],popt[2])
    norm = d4+d3+d2+d1
    (d1,d2,d3,d4) = (d1/norm, d2/norm, d3/norm, d4/norm)

    R = np.arange(0,2.001,.05)
    normed_n_within = np.load('normed_n_within.npy')
    normed_n_within_lum = np.load('normed_n_within_lum.npy')
    normed_n_within_lmc = np.load('normed_n_within_lmc.npy')
    plt.scatter(R, normed_n_within, label='all subhalos', marker='*')
    plt.scatter(R, normed_n_within_lum, label='subhalos $> 10^8$', marker='o')
    plt.scatter(R, normed_n_within_lmc, label='lmc subhalos $> 10^8$', marker='^')


    rr = np.arange(0.0,2,.01)
    plt.plot(rr, logistic(rr, popt[0],popt[1],popt[2]) / logistic_norm  , color='magenta',ls='--',label='logistic fit')
    plt.plot(rr,d4+d3*rr+d2*rr**2+d1*rr**3,color='cyan',ls='--', label='polynomial fit')
    plt.scatter(radius,Khost2 / 0.001133 ,color='red',label='host no cut', marker='x')
    plt.legend(loc='upper left')
    plt.xlabel('radius [kpc]')
    plt.xlim((0.0,2))
    plt.savefig('radialfit_morecareful')
    plt.xlim((0,1))
    plt.ylim((0,1))
    plt.grid(alpha=0.6)
    plt.savefig('radialfit_morecareful_1rvir')
    plt.close()



# beyond ratio = 1, values are extrapolated
def get_normed_radial_abundance(ratio, version=''):
    R = np.arange(0,2.001,.05)
    if version == '':
        path='/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/SHMF/normed_radial_survivors.npy'
    if version == 'classical':
        path='/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/SHMF/normed_radial_classical.npy'
    if version=='UFD':
        path='/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/SHMF/normed_radial_UFD_noreion.npy'
    if version == 'z0':
        path='/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/SHMF/normed_n_within_lum.npy'


    if version=='classical_noreion':
        path='/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/SHMF/normed_n_within_high.npy'
    if version=='UFD_noreion':
        path='/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/SHMF/normed_n_within_low.npy'
    if version=='LMC':
        path='/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/SHMF/normed_radial_survivors_LMC.npy'   
    if version=='dwarf':
        path='/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/SHMF/normed_radial_survivors_DWARF.npy'   


    normed_n_within_lum = np.load(path)
    normed_n_within_lum[20:] = .4*R[20:]+.6

    tck =  interpolate.splrep(R,normed_n_within_lum, k=1)
    return interpolate.splev(ratio,tck)

# no std provided for ratio > 1
def get_normed_radial_std(ratio, version=''):
    R = np.arange(0,2.001,.05)
    if version=='':
        path = '/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/SHMF/normed_radial_survivors_std.npy'
    if version=='LMC':
        path='/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/SHMF/normed_radial_survivors_std_LMC.npy'
    if version=='dwarf':
        path='/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/SHMF/normed_radial_survivors_std_DWARF.npy'

    n_within_lum_std = np.load(path)
    tck = interpolate.splrep(R,n_within_lum_std, k=1)
    return interpolate.splev(ratio,tck)


# must multiply by Ngreater(rvir)/(4pir**2) to get density
def get_normed_sat_derivative(ratio, rvir):
    R = np.arange(0,2.001,.05)
    normed_n_within_lum = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/SHMF/normed_radial_survivors.npy')
    normed_n_within_lum[20:] = .4*R[20:]+.6
    tck = interpolate.splrep(R*rvir,normed_n_within_lum)
    return interpolate.splev(ratio*rvir,tck, der=1)



# next do the same but for LMC sized halos. perhaps the number of galaxies near the center is higher.
def get_cum_radial_hist():
    hpaths = dm.get_hpaths(field=False, lx=14) 
    radial_ratio = []; radial_ratio_lum = []
    radial_lmcsized = []

    for hpath in hpaths:
        snap_z0 = htils.get_numsnaps(hpath)-1
        cat=htils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID = htils.load_zoomid(hpath)
        radius = cat.ix[hostID]['rvir']/cat.h0   #probably should divide by h0
        mvir = cat.ix[hostID]['mvir']/cat.h0
    
        dists = dm.distance(np.array(cat[['posX','posY','posZ']]), np.array(cat.ix[hostID][['posX','posY','posZ']]))*1000/cat.h0
        mask = (dists< 2*radius)*(dists>0) # exclude host 
        mask2 = np.logical_or(cat['hostID']==hostID, cat['hostID']==-1) # only take 1 level deep halos
        mass_mask = np.array(cat['mgrav']/cat.h0 > 10**8)
        mask3 = np.array(mask & mask2)
        radial_ratio = np.append(radial_ratio, dists[mask3]/radius)
        
        # NEXT WORK IS HERE
        radial_ratio_lum.append(dists[mask3 & mass_mask]/radius)
        #np.append(radial_ratio_lum, dists[mask3 & mass_mask]/radius)
        # instead of appending with numpy, append with regular. then do the n_within ratio for all, then find variance

        ## Now compute same thing but on LMCish sized hosts
        field_halos = fshmf.get_field_halos(cat,hpath,hostID,mlow=np.log10(5*10**10),mhigh=np.log10(5*10**11))
        for i in range(len(field_halos)):
            field_halo = field_halos[i:i+1]
            radius = float(field_halo['rvir']/cat.h0)
            dists = dm.distance(np.array(cat[['posX','posY','posZ']]), np.array(field_halo[['posX','posY','posZ']]))*1000/cat.h0
            mask = (dists< 2*radius)*(dists>0) # exclude host 
            mask2 = np.logical_or(cat['hostID']==int(field_halo['id']), cat['hostID']==-1)
            mask3 = np.array(mask & mask2)
            radial_lmcsized = np.append(radial_lmcsized, dists[mask3 & mass_mask]/radius)

        print 'done with cat', htils.hpath_catnum(hpath)

    # now make cum hist of it
    R = np.arange(0,2.001,.05)
    n_within = np.array([np.sum(radial_ratio < R[i]) for i in range(len(R))])
    norm = float(np.sum(radial_ratio < 1))
    normed_n_within = n_within / norm
    np.save('normed_n_within', normed_n_within)
    
    n_within_lum = np.array([[np.sum(ratio < R[i]) for i in range(len(R))] for ratio in radial_ratio_lum])
    n_within_lum_sum = np.sum(n_within_lum, axis=0)
    normed_n_within_lum = n_within_lum_sum / float(n_within_lum_sum[20])
    n_within_lum_norm = np.array([n_within_lum[i]/float(n_within_lum[i][20]) for i in range(len(hpaths))])
    n_within_lum_std = np.std(n_within_lum_norm, axis=0)
    np.save('normed_n_within_lum', normed_n_within_lum)
    np.save('n_within_lum_std', n_within_lum_std)

    """
    n_within_lum = np.array([np.sum(radial_ratio_lum < R[i]) for i in range(len(R))])
    norm = float(np.sum(radial_ratio_lum < 1))
    normed_n_within_lum = n_within_lum / norm
    np.save('normed_n_within_lum', normed_n_within_lum)
    """
    
    # replace this in a similar manner to get the standard deviations.
    # are the about the same?
    n_within_lmc = np.array([np.sum(radial_lmcsized < R[i]) for i in range(len(R))])
    norm = float(np.sum(radial_lmcsized < 1))
    normed_n_within_lmc = n_within_lmc / norm
    np.save('normed_n_within_lmc', normed_n_within_lmc)
    


def get_cum_radial_hist_diffmass():
    hpaths = dm.get_hpaths(field=False, lx=14)
    radial_ratio_lowmass = []
    radial_ratio_highmass = []
    radial_ratio_lum = []
    rr_67 = []; rr_78=[]; rr_89=[]; rr_9up=[]
    rr_inMT=[]; rr_notMT=[]

    
    import abundance_matching as am
    model = am.GK16_grow(reionization=True, hostSHMF=True)
    mlow = model.stellar_to_halo_mass(10**2.6)
    mhigh = model.mstar_to_mhalo(10**6)

    for hpath in hpaths:
        snap_z0 = htils.get_numsnaps(hpath)-1
        cat=htils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID = htils.load_zoomid(hpath)
        radius = cat.ix[hostID]['rvir']/cat.h0   #probably should divide by h0
        mvir = cat.ix[hostID]['mvir']/cat.h0

        data=dm.get_extant_data(hpath,field=False)
        rsids = np.array(data['rsid'])
        dists = dm.distance(np.array(cat.ix[rsids][['posX','posY','posZ']]), np.array(cat.ix[hostID][['posX','posY','posZ']]))*1000/cat.h0        

        mask = (dists< radius)*(dists>0) # exclude host 
        mask2 = np.logical_or(cat.ix[rsids]['hostID']==hostID, cat.ix[rsids]['hostID']==-1) # only take 1 level deep halos
        bigmass = np.array(data['peak_mvir']/cat.h0 > 2*10**9)
        lowmass = np.array(data['peak_mvir']/cat.h0 < 2*10**9)
        lum_mask = np.array(  (data['peak_mvir']/cat.h0 < mhigh)& (data['peak_mvir']/cat.h0 > mlow))

        submasses = np.array(cat.ix[rsids]['mgrav']/cat.h0)
        mask67 = np.array((submasses > 10**6)&(submasses < 10**7))
        mask78 = np.array((submasses > 10**7)&(submasses < 10**8))
        mask89 = np.array((submasses > 10**8)&(submasses < 10**9))
        mask9up = np.array((submasses > 10**9))
        MTmask = np.array((submasses > 10**6.1)&(submasses < 10**9.4))
        print np.log10(np.median(data['peak_mvir'][lum_mask]/cat.h0)), np.log10(np.median(data['peak_mvir'][MTmask]/cat.h0)), 'mean mass in lum mask and in MTmask'
        print len(lum_mask), len(MTmask), np.log10(mlow), 'lengthts of mask should be =. low end mass of lum mask'
        ## Include halos not in MT
        distsFULL = dm.distance(np.array(cat[['posX','posY','posZ']]), np.array(cat.ix[hostID][['posX','posY','posZ']]))*1000/cat.h0    
        maskFULL = (distsFULL< radius)*(distsFULL>0) # exclude host 
        mask2FULL = np.logical_or(cat['hostID']==hostID, cat['hostID']==-1) # only take 1 
        submassesFULL = np.array(cat['mgrav']/cat.h0)
        MTmaskFULL = np.array((submassesFULL > 10**6.1)&(submassesFULL < 10**9.4))
        ####
        

        mask_big = np.array(bigmass & mask & mask2)
        mask_low = np.array(lowmass & mask & mask2)
        mask_lum = np.array(lum_mask & mask & mask2)

        print np.sum(mask_lum), np.sum(np.array(MTmask&mask&mask2)), np.sum(np.array(MTmaskFULL&maskFULL&mask2FULL)), 'num samples in lum, z=0, z=0FULL'

        radial_ratio_lowmass = np.append(radial_ratio_lowmass, dists[mask_low]/radius)
        radial_ratio_highmass = np.append(radial_ratio_highmass, dists[mask_big]/radius) 
        radial_ratio_lum.append(dists[mask_lum]/radius)
        rr_67 = np.append(rr_67, dists[np.array(mask67&mask&mask2)]/radius)
        rr_78 = np.append(rr_78, dists[np.array(mask78&mask&mask2)]/radius)
        rr_89 = np.append(rr_89, dists[np.array(mask89&mask&mask2)]/radius)
        rr_9up = np.append(rr_9up, dists[np.array(mask9up&mask&mask2)]/radius)
        rr_inMT = np.append(rr_inMT, dists[np.array(MTmask&mask&mask2)]/radius)
        rr_notMT = np.append(rr_notMT, distsFULL[np.array(MTmaskFULL&maskFULL&mask2FULL)]/radius)


    R = np.arange(0,2.001,.05)
    n_within_low = np.array([np.sum(radial_ratio_lowmass < R[i]) for i in range(len(R))])
    n_within_high = np.array([np.sum(radial_ratio_highmass < R[i]) for i in range(len(R))])
    nw_67 = np.array([np.sum(rr_67 < R[i]) for i in range(len(R))])
    nw_78 = np.array([np.sum(rr_78 < R[i]) for i in range(len(R))])
    nw_89 = np.array([np.sum(rr_89 < R[i]) for i in range(len(R))])
    nw_9up = np.array([np.sum(rr_9up < R[i]) for i in range(len(R))])
    nw_inMT = np.array([np.sum(rr_inMT < R[i]) for i in range(len(R))])
    nw_notMT = np.array([np.sum(rr_notMT < R[i]) for i in range(len(R))])

    norm_low = float(np.sum(radial_ratio_lowmass < 1))
    norm_high = float(np.sum(radial_ratio_highmass < 1))
    norm_67 = float(np.sum(rr_67 < 1))
    norm_78 = float(np.sum(rr_78 < 1))
    norm_89 = float(np.sum(rr_89 < 1))
    norm_9up = float(np.sum(rr_9up < 1))
    norm_inMT = float(np.sum(rr_inMT < 1))
    norm_notMT = float(np.sum(rr_notMT < 1))

    normed_n_within_low = n_within_low / norm_low
    normed_n_within_high = n_within_high / norm_high

    np.save('normed_n_within_low', normed_n_within_low)
    np.save('normed_n_within_high', normed_n_within_high)
    np.save('nw_67', nw_67/norm_67)
    np.save('nw_78', nw_78/norm_78)
    np.save('nw_89', nw_89/norm_89)
    np.save('nw_9up', nw_9up/norm_9up)
    np.save('nw_inMT', nw_inMT/norm_inMT)

    #print nw_notMT
    #print norm_notMT
    #print nw_notMT/norm_notMT
    np.save('nw_notMT', nw_notMT/norm_notMT)


    n_within_lum = np.array([[np.sum(ratio < R[i]) for i in range(len(R))] for ratio in radial_ratio_lum])
    n_within_lum_sum = np.sum(n_within_lum, axis=0)
    normed_n_within_lum = n_within_lum_sum / float(n_within_lum_sum[20])
    n_within_lum_norm = np.array([n_within_lum[i]/float(n_within_lum[i][20]) for i in range(len(hpaths))])
    n_within_lum_std = np.std(n_within_lum_norm, axis=0)
    np.save('normed_n_within_lum_peakmass', normed_n_within_lum)
    np.save('n_within_lum_std_peakmass', n_within_lum_std)
    



def get_cum_radial_hist_reionization():
    hpaths = dm.get_hpaths(field=False, lx=14)
    radial_ratio_lum = []
    
    import abundance_matching as am
    model = am.GK16_grow(reionization=True, hostSHMF=True)
    mlow = model.stellar_to_halo_mass(10**2.6)
    mhigh = model.mstar_to_mhalo(10**6)

    for hpath in hpaths:
        snap_z0 = htils.get_numsnaps(hpath)-1
        cat=htils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID = htils.load_zoomid(hpath)
        radius = cat.ix[hostID]['rvir']/cat.h0   #probably should divide by h0
        mvir = cat.ix[hostID]['mvir']/cat.h0
        data=dm.get_extant_data(hpath,field=False)

        # now only get rsids of halos that survive reionization #        
        vmax_ach = 9.48535156 
        vmax_filt = 23.54248047
        mask_pre = np.array(data['vmax_12'] > vmax_ach)  # big enough at reionization
        mask_post = np.array(data['peak_vmax'] >= vmax_filt)  # or big enough later
        rsids = np.array(data['rsid'][mask_pre | mask_post])
        print np.sum(mask_pre | mask_post), '/', len(mask_pre), 'survive reionization'


        dists = dm.distance(np.array(cat.ix[rsids][['posX','posY','posZ']]), np.array(cat.ix[hostID][['posX','posY','posZ']]))*1000/cat.h0        
        mask = (dists< 2*radius)*(dists>0) # exclude host 
        #mask2 = np.logical_or(cat.ix[rsids]['hostID']==hostID, cat.ix[rsids]['hostID']==-1) # only take 1 level deep halos
        #lum_mask = np.array(  (data['peak_mvir']/cat.h0 < mhigh)& (data['peak_mvir']/cat.h0 > mlow))
        #mask_lum = np.array(lum_mask & mask & mask2)
        radial_ratio_lum.append(dists[mask]/radius)

    R = np.arange(0,2.001,.05)
    n_within_lum = np.array([[np.sum(ratio < R[i]) for i in range(len(R))] for ratio in radial_ratio_lum])
    n_within_lum_sum = np.sum(n_within_lum, axis=0)
    normed_n_within_lum = n_within_lum_sum / float(n_within_lum_sum[20])
    n_within_lum_norm = np.array([n_within_lum[i]/float(n_within_lum[i][20]) for i in range(len(hpaths))])
    n_within_lum_std = np.std(n_within_lum_norm, axis=0)
    np.save('normed_radial_survivors', normed_n_within_lum)
    np.save('normed_radial_survivors_std', n_within_lum_std)




## for LMC - choose 10^10.7 - 10^11.7
## for Isolatd Dwarfs in 1st paper - 10^6 - 10^8 Msun in stellar mass = 9.618 to 10.69

def get_cum_radial_hist_LMC():
    hpaths = dm.get_hpaths(field=False, lx=14)
    radial_ratio_lum = []
    
    import abundance_matching as am
    model = am.GK16_grow(reionization=True, hostSHMF=True)
    #mlow = model.stellar_to_halo_mass(10**2.6)
    #mhigh = model.mstar_to_mhalo(10**6)

    for hpath in hpaths:
        snap_z0 = htils.get_numsnaps(hpath)-1
        cat=htils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID = htils.load_zoomid(hpath)
        field_halos = get_field_halos(cat,hpath,hostID,mlow=10,mhigh=10.7)  #MANUALLY ADJUST FOR FIELD SIZES
        for i in range(len(field_halos)):
            field_halo=field_halos[i:i+1]
            fieldID = int(field_halo['id'])
            data = dm.get_extant_data(hpath, field=True)
            mask = data['field_rsid']== fieldID
            data = data[mask]
    
            radius = cat.ix[fieldID]['rvir']/cat.h0   #probably should divide by h0
            mvir = cat.ix[fieldID]['mvir']/cat.h0

            # now only get rsids of halos that survive reionization #        
            vmax_ach = 9.48535156 
            vmax_filt = 23.54248047 * 0.95
            mask_pre = np.array(data['vmax_12'] > vmax_ach)
            #mask_post = np.array(data['peak_vmax'] >= vmax_filt)
            mask_post = np.array(data['max_mass_vmax'] >= vmax_filt)
            rsids = np.array(data['rsid'][mask_pre | mask_post])
            print np.sum(mask_pre | mask_post), '/', len(mask_pre), 'survive reionization'

            dists = dm.distance(np.array(cat.ix[rsids][['posX','posY','posZ']]), np.array(cat.ix[fieldID][['posX','posY','posZ']]))*1000/cat.h0 
            mask = (dists< 2*radius)*(dists>0) # exclude host 
            #lum_mask = np.array(  (data['peak_mvir']/cat.h0 < mhigh)& (data['peak_mvir']/cat.h0 > mlow))
            #mask_lum = np.array(lum_mask & mask & mask2)
            radial_ratio_lum.append(dists[mask]/radius)

    R = np.arange(0,2.001,.05)
    n_within_lum = np.array([[np.sum(ratio < R[i]) for i in range(len(R))] for ratio in radial_ratio_lum])
    n_within_lum_sum = np.sum(n_within_lum, axis=0)
    normed_n_within_lum = n_within_lum_sum / float(n_within_lum_sum[20])
    n_within_lum_norm = np.array([n_within_lum[i]/float(n_within_lum[i][20]) for i in range(len(n_within_lum))])
    n_within_lum_std = np.std(n_within_lum_norm, axis=0)
    np.save('normed_radial_survivors_DWARF', normed_n_within_lum)
    np.save('normed_radial_survivors_std_DWARF', n_within_lum_std)




def plot_radial_fits_massdiff():
    R = np.arange(0,2.001,.05)
    normed_n_within_low = np.load('normed_n_within_low.npy')
    normed_n_within_high = np.load('normed_n_within_high.npy')
    normed_n_within_lum_peak = np.load('normed_n_within_lum_peakmass.npy')
    normed_n_within_lum = np.load('normed_n_within_lum.npy')

    nw_67 = np.load('nw_67.npy')
    nw_78 = np.load('nw_78.npy')
    nw_89 = np.load('nw_89.npy')
    nw_9up = np.load('nw_9up.npy')
    nw_inMT = np.load('nw_inMT.npy')
    nw_notMT = np.load('nw_notMT.npy')


    plt.scatter(R, normed_n_within_low, label='subhalos $< 10^{8.5}$', marker='*')
    plt.scatter(R, normed_n_within_high, label='subhalos $> 10^9$', marker='o')
    plt.scatter(R, normed_n_within_lum_peak, label='lum peak mass', marker='x')
    plt.scatter(R, normed_n_within_lum, label='$>10^8 at z=0$', marker='^')

    plt.plot(R, nw_67, label='6-7')
    plt.plot(R, nw_78, label='7-8')
    plt.plot(R, nw_89, label='8-9')
    plt.plot(R, nw_9up, label='9up')
    plt.plot(R, nw_inMT, label='in MT')
    plt.plot(R, nw_notMT, label='Full Halos')

    rr = np.arange(0.0,2,.01)
    plt.legend(loc='upper left')
    plt.xlabel('radius [kpc]')
    plt.xlim((0,1))
    plt.ylim((0,1))
    plt.grid(alpha=0.6)
    plt.savefig('radialfit_1rvir_massdiff')
    plt.close()


def plot_radial_fits_reionization():
    R = np.arange(0,2.001,.05)
    normed_n_within_lum_peak = np.load('normed_n_within_lum_peakmass.npy')
    normed_n_within_lum = np.load('normed_n_within_lum.npy')
    normed_reion = np.load('normed_radial_survivors.npy')
    
    # fill in values from 1 - 2 for normed reion
    normed_reion[20:] = .4*R[20:]+.6
    print normed_reion

    (c1,c2,c3,c4) = np.polyfit(R, normed_reion, 3)
    norm = c4+c3+c2+c1
    #print c4/norm, c3/norm, c2/norm, c1/norm, 'normed constants'
    plt.plot(R,c4+c3*R+c2*R**2+c1*R**3,color='green',ls='--')

    plt.scatter(R, normed_n_within_lum_peak, label='lum peak mass', marker='x')
    plt.scatter(R, normed_n_within_lum, label='$>10^8 at z=0$', marker='^')
    plt.scatter(R, normed_reion, label='reionization survivors', marker='o')

    # try a logistics fit
    #def logistic(x, a,b,k):
    #    return a/(1+ b* np.e**(-k * (x)))

    def logistic(x, b,k):
        a = (1+b)*(1-np.e**-k)/ (1-np.e**-k * (1/b + 1))
        return a/(1+ (a/b-1) * np.e**(-k * (x))) - b
    
    popt, pcov = curve_fit(logistic, R, normed_reion)
    plt.plot(R, logistic(R, popt[0],popt[1]), color='magenta',ls='--',label='logistic fit')

    tck =  interpolate.splrep(R,normed_reion, k=1)
    plt.plot(R, interpolate.splev(R,tck), label='interp')


    rr = np.arange(0.0,2,.01)
    plt.legend(loc='upper left')
    plt.xlabel('radius [kpc]')
    plt.xlim((0,1))
    plt.ylim((0,1.))
    plt.grid(alpha=0.6)
    plt.savefig('radialfit_1rvir_reionization')
    plt.close()





#get_cum_radial_hist()

#plot_radial_fits()


# FIGURE 7:
#plot_los_frac()


"""
#plot_normed_fit()
K = get_normed_los(1./np.sqrt(2),1./np.sqrt(2))
print K, 'should be less than 1 and greater than 0.26'
K = get_normed_los(1.,1.)
print K,'should be less than 1.5 and greater than 1'
K = get_normed_los(0.001,0.001)
print  K, 'should get close to 0'
"""

#K = get_normed_los(0.5,1)
#print K, 'value for half of the virial radius'



