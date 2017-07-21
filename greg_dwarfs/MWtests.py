import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import abundance_matching as am
from PlotParams import *
import DwarfMethods as dm
import RadialDependence as rd
from scipy.integrate import dblquad
from scipy.integrate import quad
import AMbase as ambase
from matplotlib import gridspec

LMC_stellar_mass = 2.6*10**9
SMC_stellar_mass = 7.1*10**8


vmax_ach = 9.48535156   # multiply by 0.6 for GK16 with z=13.  For z=9, Moster use 0.8, GK16 use 0.85
vmax_filt = 23.54248047  # has little effect
z = 11  # choices of 14,13,12,11,9

# plot number of observed galaxies in MW > 10^5

def MWngreater():    # ngreater_v_minlum_2panel
    #w=8; h=6
    #fig, (ax1, ax2) = plt.subplots(nrows=2,ncols=1,sharex=True, figsize=(w,h*1.6) )
    fig, ax1 = plt.subplots()
    min_lums = 10**np.linspace(5,8,20)


    # Plot Observed Galaxies 
    galx = dm.load_nearby_gx()
    conf_galx = ~np.array([galx['Name'][i][0]=='*' for i in range(len(galx))])
    #near_LMC[(galx['Name']=='LMC')|(galx['Name']=='SMC')] = False     # take out the LMC and SMC from the list
    smass = galx['mstar']
    #conf_smass = galx['mstar'][near_LMC & conf_galx]
    n_greater = np.array([np.sum(smass > min_lums[i]) for i in range(len(min_lums))])
    #n_greater_conf = np.array([np.sum(conf_smass > min_lums[i]) for i in range(len(min_lums))])
    ax1.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Known Satellites', color='black',lw=linewidth+1,linestyle='-')


    # FIRST DO TOP PANEL
    for model in [am.GarrisonKimmel(hostSHMF=True,reionization=True,  catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt ), am.GK16_grow(hostSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt),  am.Moster(hostSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt), am.Brook(hostSHMF=True,reionization=True,catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt), am.GK16_grow(hostSHMF=True,reionization=True, catreion=True,z=12, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma = -0.5), am.GK16_grow(hostSHMF=True,reionization=True, catreion=True,z=7, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma = -1.0), am.Behroozi(hostSHMF=True,reionization=True)]:
        print model.label, 'on this model'
        halo_mass = 1.4*10**12
        halo_mass = dm.convert_Mvir(halo_mass,model)
        N_visible, N_std = model.ngreater(halo_mass,min_lums)

        if isinstance(model, am.GK16_grow):
            if model.catreion and model.z == 11: 
                ax1.fill_between(min_lums, N_visible+N_std, N_visible-N_std, facecolor=model.color, alpha=0.15)
                ax1.plot(min_lums, N_visible,label=model.label,color=model.color,lw=linewidth,ls=model.ls)
            elif model.z ==12:
                ax1.plot(min_lums, N_visible,label='Slope = -2.43,z=12',color=model.color,lw=linewidth,ls='--')
                ax1.fill_between(min_lums, N_visible+N_std, N_visible-N_std, facecolor=model.color, alpha=0.15)
            else:
                ax1.plot(min_lums, N_visible,label='Slope = -3.31,z=7',color=model.color,lw=linewidth,ls='-.')
        else:
            ax1.plot(min_lums, N_visible,label=model.label,color=model.color,lw=linewidth,ls=model.ls)

    ax1.text(.2,.90,'MW Sized Host',transform=ax1.transAxes,fontsize=legend_size-2)
    ax1.set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=label_font)
    plt.xscale('log')
    #ax1.set_ylim(ymin=0)
    ax1.set_ylim((0,35))
    

    ax1.legend(loc='upper right',frameon=False,fontsize=legend_size-2)
    ax1.set_xlabel('$\mathrm{M_*} \, [\mathrm{M_\odot}]$', fontsize=label_font)
    plt.xlim((min_lums[0],min_lums[-1]))

    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.subplots_adjust(hspace = 0.1)
    plt.gcf().subplots_adjust(bottom=0.15)
    extra=''
    plt.savefig('LMCPlots/MW_Ngreater_vs_minlum'+extra+'.png')
    plt.close()




# get proper error bars on this
def get_std(nsubs, N=2000):
    from scipy import interpolate

    R = np.arange(0,1.0001,.02)
    factor =  rd.get_normed_radial_abundance(R)
    #given factor, find R
    tck = interpolate.splrep(factor,R, k=1)

    #nsubs = 11
    random_set = []
    for i in range(N):
        randnums = np.random.rand(nsubs)
        random_r = interpolate.splev(randnums,tck)
        random_set.append(random_r)

    n_within_lum = np.array([[np.sum(ratio < R[i]) for i in range(len(R))] for ratio in random_set])
    n_within_lum_std = np.std(n_within_lum, axis=0)
    n_within_mean = np.mean(n_within_lum, axis=0)
    
    #n_within_lum_sum = np.sum(n_within_lum, axis=0)
    #normed_n_within_lum = n_within_lum_sum / float(n_within_lum_sum[-1])
    #n_within_lum_norm = np.array([n_within_lum[i]/float(n_within_lum[i][-1]) for i in range(len(n_within_lum))])
    #n_within_lum_std = np.std(n_within_lum_norm, axis=0)
    return R, n_within_lum_std, n_within_mean



### NEED TO ADD BACK LMC AND SMC
def MW_radial(re=True, min_lum=2*10**5):
    f, ax = plt.subplots()
    radius = np.arange(0,301,4)
    rvir = 300.0

    # Plot Observed Galaxies First
    galx = dm.load_nearby_gx()
    conf_galx = ~np.array([galx['Name'][i][0]=='*' for i in range(len(galx))])
    mass_cut  = galx['mstar']> min_lum
    mcut2 = galx['mstar'] < 10**8
    dists = galx['dist_GC'][mass_cut & mcut2]
    #conf_dists = galx['dist_GC'][mass_cut & conf_galx]

    ngreater = float(len(dists))
    print ngreater, 'num satellites'
    print np.sum(mass_cut), 'with the LMC and SMC'
    
    n_within = np.array([np.sum(dists < radius[i]) for i in range(len(radius))])
    #n_within_conf = np.array([np.sum(conf_dists < radius[i]) for i in range(len(radius))])
    ax.plot(dm.step_plotX(radius)/rvir, dm.step_plotY(n_within)/ngreater, label=r'Known MW Satellites ($M_* > 2 \times 10^5 \, \mathrm{M_\odot}$)'+'\nExcluding LMC & SMC', color='black', lw=linewidth)

    # Now Plot Theoretical Model
    #dm_mass = 1.4*10**12   
    model = am.GK16_grow(reionization=True, hostSHMF=True)
    #rvir = model.mass_to_rvir(dm_mass)  # this does correctly convery to virial radius from M200 for Moster, etc.


    factor =  rd.get_normed_radial_abundance(radius / rvir)
    factor_bad = rd.get_normed_radial_abundance(radius/rvir, 'z0')
    #rd.getK_protected(radius / rvir)  #  'normed_n_within_lum' is > 10^8 at z=0
    #factor_classical_noreion = rd.get_normed_radial_abundance(radius / rvir,'classical_noreion')
    #factor_UFD_noreion = rd.get_normed_radial_abundance(radius / rvir,'UFD_noreion')
    factor_LMC = rd.get_normed_radial_abundance(radius / rvir,'LMC')
    factor_dwarf = rd.get_normed_radial_abundance(radius / rvir,'dwarf')





    R, std, mean = get_std( int(ngreater) )
    #ax.plot(R*rvir,mean,color='black',lw=1,label='rand method')
    ax.fill_between(R, (mean+std) / ngreater ,(mean-std) / ngreater, facecolor=model.color, alpha=0.2)


    ax.plot(radius/rvir, factor, color=model.color,linestyle='-',lw=linewidth+2, label='Luminous Satellites (predicted)') 
    plt.plot(R, analytic_radial(R),ls=':',label='Fitting Function', lw=3.5, color='blue')
    ax.plot(radius/rvir, factor_bad, color=model.color,linestyle='--',lw=linewidth+1, label='All DM Subhalos')
    #ax.plot(radius, factor_classical_noreion*ngreater, linestyle='--',lw=linewidth+1, label='classical no reionization')
    #ax.plot(radius, factor_UFD_noreion*ngreater, linestyle='--',lw=linewidth+1, label='UFD no reionization')
    #ax.plot(radius, factor_LMC*ngreater, linestyle='--',lw=linewidth+1, label='LMC Prediction')
    #ax.plot(radius, factor_dwarf*ngreater, linestyle='--',lw=linewidth+1, label='Dwarf Prediction')


    r_sigma = rd.get_normed_radial_std(radius / rvir)
    #r_sigmaLMC = rd.get_normed_radial_std(radius / rvir, 'LMC')
    #r_sigmaDWARF = rd.get_normed_radial_std(radius / rvir,'dwarf')
    
    #print r_sigma, 'radial uncertainty'
    sig2 = ngreater*r_sigma
    full_sigma = sig2
    #ax.fill_between(radius, (factor+r_sigma)*ngreater, (factor-r_sigma)*ngreater, facecolor=model.color, alpha=0.20)
    #ax.fill_between(radius, (factor_LMC+r_sigmaLMC)*ngreater, (factor_LMC-r_sigmaLMC)*ngreater, alpha=0.20)
    #ax.fill_between(radius, (factor_dwarf+r_sigmaDWARF)*ngreater, (factor_dwarf-r_sigmaDWARF)*ngreater,facecolor='green', alpha=0.20)

 
   
    ax.set_ylabel('$\mathrm{ N_{sats} / N_{tot}} > \mathrm{R}$', fontsize=label_font)
    #ax.set_xlabel('Radius [kpc]',fontsize=label_font)
    ax.set_xlabel('$R \equiv r/R_{\mathrm{vir}}$',fontsize=label_font)
    #ax.text(.04,.93,'$M_* > 2 \times 10^5 \, \mathrm{M_\odot}$ in MW-like host',transform=ax.transAxes,fontsize=legend_size-3)

    ax.tick_params(axis='both', which='major', labelsize=tick_size)
    ax.legend(loc=(0.02,.67),frameon=False,fontsize=legend_size-3.8)
    ax.set_ylim((0,1.05))
    #ax.set_xlim((0,300))
    ax.set_xlim((0,1))
    
    plt.gcf().subplots_adjust(bottom=0.14)
    plt.savefig('MWradial_test.pdf')
    plt.close()



def inversetan(x, const1, const2, const3):
    return const2 * np.arctan(x/const1 - const3)

def analytic_radial(R):
    d4 = 0.0 # 0.000212585034013
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


def plot_radial_fits_reionization():
    from scipy.optimize import curve_fit
    from scipy import interpolate

    R = np.arange(0,2.001,.05)
    factor =  rd.get_normed_radial_abundance(R)
    plt.plot(R, factor, label='Luminous Satellites', lw=5)

    def inversetan_func(x, const1, intercept):
        #const1 = 0.5; intercept = 0.12
        const3 = intercept / const1
        const2 = 1. / np.arctan(1/const1 - const3)
        return const2 * np.arctan(x/const1 - const3)

    popt, pcov = curve_fit(inversetan_func, R[4:26], factor[4:26], p0=(0.5,0.12))
    print popt[0], popt[1], 'const1, intercept values'
    print popt[0], popt[1]/popt[0], 1./ np.arctan(1/popt[0] - popt[1]/popt[0]), 'const1, const3, const2'
    plt.plot(R, inversetan_func(R, popt[0],popt[1]),ls='--',label='arctanh fit', lw=2.5)

    start = R[0:5]
    (d1,d2,d3,d4) = np.polyfit(start, factor[0:5], 3)
    print d4,d3,d2,d1, 'ds'
    #print c4/norm, c3/norm, c2/norm, c1/norm, 'normed constants'
    #plt.plot(start,d4+d3*start+d2*start**2+d1*start**3,color='red',ls='--',lw=2.5,label='polynomial fit')
    plt.plot(R,d4+d3*R+d2*R**2+d1*R**3,color='red',ls='--',lw=2.5,label='polynomial fit')
    plt.plot(R, analytic_radial(R),ls=':',label='analytic splice', lw=2.5, color='orange')

    rr = np.arange(0.0,2,.01)
    plt.legend(loc='upper left')
    plt.xlabel('radius [kpc]')
    plt.xlim((0,1.3))
    plt.ylim((0,1.25))
    plt.grid(alpha=0.6)
    plt.savefig('LMCPlots/RadialFit')
    plt.close()









"""
# To fit the function, need a spliced function
# start with 3rd order polynomial to start - 0- 0.2
# Then go to the arctan function
const1 = 0.5; intercept = 0.12
const3 = intercept / const1
const2 = 1. / np.arctan(1/const1 - const3)
ax.plot(R*rvir, ngreater * const2 * np.arctan(R / const1 - const3),color='green',lw=2)

k1 = 0.2 /(0.3**3)
ax.plot(R*rvir, ngreater * k1 * R**3, color='black')
"""

"""
const1 = 0.65; intercept = 0.1
const3 = intercept / const1
const2 = 1. / np.arctan(1/const1 - const3)
ax.plot(R*rvir, ngreater * const2 * np.arctan(R / const1 - const3),color='blue',lw=2)

const1 = 0.4; intercept = 0.15
const3 = intercept / const1
const2 = 1. / np.arctan(1/const1 - const3)
ax.plot(R*rvir, ngreater * const2 * np.arctan(R / const1 - const3),color='red',lw=2)
"""
