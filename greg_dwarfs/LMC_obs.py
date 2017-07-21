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


### NEED TO ADD BACK LMC AND SMC
def plotMW_radial(re=True, min_lum=10**3):
    f, (ax2, ax1) = plt.subplots(nrows=2, figsize=(8,9))
    #fig = plt.figure(figsize=(8,9))
    #gs = gridspec.GridSpec(2,1, height_ratios=[1,2])
    #ax1 = plt.subplot(gs[0])
    #ax2 = plt.subplot(gs[1])
    
    radius = np.arange(0,301,4)
    # Compute the stellar mass observable as a function of m = 19
    #m = 19
    #M = m - 5*(np.log10(radius*1000)-1)

    M_DES = (1.45 - np.log10(radius))/0.228
    M_DR5 = (1.1 - np.log10(radius))/0.228

    #print zip(radius,M)
    min_mstar_DES = 10**( (4.86-M_DES)/2.5)
    min_mstar_DR5 = 10**( (4.86-M_DR5)/2.5)
    #print zip(radius,min_mstar)
    min_mstar_DES[min_mstar_DES < 100] = 100  # impose a lower limit.
    min_mstar_DR5[min_mstar_DR5 < 100] = 100  # impose a lower limit.
    #min_mstar = np.array([10**3]*len(M))
    

    # Plot Observed Galaxies First
    galx = dm.load_nearby_gx()
    conf_galx = ~np.array([galx['Name'][i][0]=='*' for i in range(len(galx))])
    mass_cut  = galx['mstar']>10**3
    near_LMC = galx['dist_LMC'] < 50   # I need real tidal radius of LMC/SCM. i.e. what is sub-substructure according to rockstar?
    near_SMC = galx['dist_SMC'] < 35   # use my SIDM paper to find this tidal radius for LMC and SMC
    near_LMC[(galx['Name']=='LMC')|(galx['Name']=='SMC')] = False # add back in the SMC and LMC themselves
    near_SMC[(galx['Name']=='LMC')|(galx['Name']=='SMC')] = False # add back in the SMC and LMC themselves
    print np.sum(near_SMC | near_LMC), 'num in magellanic clouds'
    
    
    dists = galx['dist_GC'][mass_cut]
    conf_dists = galx['dist_GC'][mass_cut & conf_galx]
    conf_dists_noLMC = galx['dist_GC'][ (mass_cut & conf_galx) &  ~(near_SMC | near_LMC) ] # add back in the SMC and LMC!
    dists_noLMC = galx['dist_GC'][ mass_cut &  ~(near_SMC | near_LMC) ] 
    
    
    n_within = np.array([np.sum(dists < radius[i]) for i in range(len(radius))])
    n_within_conf = np.array([np.sum(conf_dists < radius[i]) for i in range(len(radius))])
    n_within_conf_noLMC = np.array([np.sum(conf_dists_noLMC < radius[i]) for i in range(len(radius))])
    n_within_noLMC = np.array([np.sum(dists_noLMC < radius[i]) for i in range(len(radius))])
    
    ax2.plot(dm.step_plotX(radius), dm.step_plotY(n_within_conf_noLMC), label='Confirmed Galaxies', color='black', lw=linewidth)
    ax2.plot(dm.step_plotX(radius), dm.step_plotY(n_within_noLMC), label='All Candidate Galaxies',color='black', linestyle='--', lw=linewidth)
    #ax1.plot(dm.step_plotX(radius), dm.step_plotY(n_within_conf_noLMC), label='Confirmed', color='black', lw=linewidth)
    #ax1.plot(dm.step_plotX(radius), dm.step_plotY(n_within_noLMC), label='Confirmed + Unconfirmed',color='black', linestyle='--', lw=linewidth)


    # Now Plot Theoretical Model
    vmax_ach = 9.48535156   # multiply by 0.6 for GK16 with z=13.  For z=9, Moster use 0.8, GK16 use 0.85
    vmax_filt = 23.54248047  # has little effect
    z = 11  # choices of 14,13,12,11,9

    # I should plot a few different choices of pg, reionization, tidal stripping that best fits the LMC data

    for model, ls, label in zip([am.GK16_grow(reionization=re, hostSHMF=True), am.GK16_grow(reionization=re, hostSHMF=True, catreion=True,z=11, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=-0.5), am.GK16_grow(reionization=re, hostSHMF=True, catreion=True,z=7, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=-1.0)], ['-','--','-.'], ['$\gamma = -0.2, z=13.3$', '$\gamma = -0.5, z=11.3$','$\gamma = -1.0, z=7.3$' ]): 

        #,am.Brook(reionization=re, lmcSHMF=True,catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt)]:, am.GarrisonKimmel(reionization=re, hostSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt),    am.Moster(reionization=re, hostSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt)
        dm_mass = 1.4*10**12   # model.mstar_to_mhalo(LMC_stellar_mass,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        rvir = model.mass_to_rvir(dm_mass)  # this does correctly convery to virial radius from M200 for Moster, etc.

        factor =  rd.get_normed_radial_abundance(radius / rvir)
        ngreater, N_std, p20, p80 = model.ngreater(dm_mass,min_lum=min_lum, percentile=True)


        #ax1.plot(radius, factor*ngreater, color=model.color,linestyle=model.ls,lw=linewidth, label=model.label)
        if isinstance(model, am.GK16_grow):
            ax2.plot(radius, factor*ngreater, color=model.color,linestyle=ls,lw=linewidth+1, label= label)

            #print N_std, model.superPoissonian, 'standard deviation of mean within rvir, and is superpoisson on?'
            r_sigma = rd.get_normed_radial_std(radius / rvir)
            #sig1 = factor*N_std
            sig2 = ngreater*r_sigma
            #full_sigma = np.sqrt( sig1**2 + sig2**2)
            #ax2.fill_between(radius_small, factor*(ngreater-N_std), factor*(ngreater+N_std), facecolor=model.color, alpha=0.25)
            #ax2.fill_between(radius, factor*p80, factor*p20, facecolor=model.color, alpha=0.25)

            # compute 1 sigma directly from N using the Boylan-Kolchin variance!
            N = factor*ngreater
            sig1 = np.sqrt(N + 0.14**2 * N**2)
            full_sigma = np.sqrt( sig1**2 + sig2**2)
            ax2.fill_between(radius, factor*ngreater+full_sigma, factor*ngreater-full_sigma, facecolor=model.color, alpha=0.20)
            #ax1.fill_between(radius, factor*ngreater+full_sigma, factor*ngreater-full_sigma, facecolor=model.color, alpha=0.25)


            # copmute with apparent mag = m limit
            ngreater_mag, N_std_mag = model.ngreater(dm_mass,min_lum=min_mstar_DES, percentile=False)
            diff_factor = factor - np.append(0,factor[:-1])
            diff_ngreater = diff_factor*ngreater_mag
            cum_ngreater = np.cumsum(diff_ngreater)
            ax1.plot(radius, cum_ngreater, color='blue',linestyle=ls,lw=linewidth+1, label='DES Sensitivity')

            ngreater_mag, N_std_mag = model.ngreater(dm_mass,min_lum=min_mstar_DR5, percentile=False)
            diff_factor = factor - np.append(0,factor[:-1])
            diff_ngreater = diff_factor*ngreater_mag
            cum_ngreater = np.cumsum(diff_ngreater)
            ax1.plot(radius, cum_ngreater, color=model.color,linestyle=ls,lw=linewidth+1, label='DR5 Sensitivity')

        else:
            ax2.plot(radius, factor*ngreater, color='red',linestyle=model.ls,lw=linewidth, label=model.label)
        print 'done with', model.label

     
    """
    ### INCLUDE ORIGINAL REIONIZATION MODEL 
    # Now Plot Theoretical Model
    for model in [am.GK16_grow(hostSHMF=True,reionization=True)]:
        dm_mass = 1.4*10**12   # model.mstar_to_mhalo(LMC_stellar_mass,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        rvir = model.mass_to_rvir(dm_mass)  # this is verified to be correct
        factor =  rd.get_normed_radial_abundance(radius / rvir)
        ngreater, N_std, p20, p80 = model.ngreater(dm_mass,min_lum=min_lum, percentile=True)
        ax2.plot(radius, factor*ngreater, color=model.color,linestyle='--',lw=linewidth, label='GK16 ($z_{\mathrm{reion}} = 13.3$)')
        #ax1.plot(radius, factor*ngreater, color=model.color,linestyle='--',lw=linewidth, label='GK16 + Late Reion')

        r_sigma = rd.get_normed_radial_std(radius / rvir)
        #sig1 = factor*N_std
        sig2 = ngreater*r_sigma
        N = factor*ngreater
        sig1 = np.sqrt(N + 0.14**2 * N**2)
        full_sigma = np.sqrt( sig1**2 + sig2**2)
        #ax2.fill_between(radius, factor*ngreater+full_sigma, factor*ngreater-full_sigma, facecolor=model.color, alpha=0.25)
        print 'done with', model.label
    """


    ax2.set_ylabel('$\mathrm{ N_{sats}} > 10^3 \, \mathrm{M_\odot}$', fontsize=label_font)
    ax1.set_ylabel('$\mathrm{ N_{sats}} > \mathrm{M_v}$', fontsize=label_font)
    ax2.set_xlabel('Radius [kpc]',fontsize=label_font)
    #ax1.set_ylabel('$\mathrm{ N_{sats}} > 10^3 \, \mathrm{M_\odot}$', fontsize=label_font)
    #ax1.set_xlabel('Radius [kpc]',fontsize=label_font)

    ax2.text(.04,.94,'$M_* > 10^3 \, \mathrm{M_\odot}$ in MW-like host',transform=ax2.transAxes,fontsize=legend_size-3)
    #ax1.text(.04,.94,'$M_* > 10^3 \, \mathrm{M_\odot}$ in MW-like host',transform=ax2.transAxes,fontsize=legend_size-3)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    #ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    ax2.legend(loc=(0.04,.35),frameon=False,fontsize=legend_size-4.3)
    ax1.legend(loc=(0.04,.35),frameon=False,fontsize=legend_size-4.3)
    ax2.set_ylim((0,120))
    ax2.set_xlim((0,300))
    
    ax1.set_ylim((0,120))
    ax1.set_xlim((0,300))

    plt.gcf().subplots_adjust(bottom=0.14)
    plt.savefig('LMCPlots/MWradial_figure.png')



def background_sats():
    #model = am.GK16_grow(reionization=True, hostSHMF=True) 
    #model=am.Moster(reionization=True, hostSHMF=True)
    #model = am.Brook(reionization=True, hostSHMF=True) 
    #model = am.GK16_grow(reionization=True, hostSHMF=True, plotgamma=-0.5) 
    #model = am.GK16_grow(reionization=True, hostSHMF=True, plotgamma=-0.2) 
    model = am.Behroozi(reionization=True,hostSHMF=True)

    dm_mass = 1.4*10**12  # for Milky Way
    R = 50 # distance from MW to LMC. should be 50

    #dm_mass = model.mstar_to_mhalo(SMC_stellar_mass,a=1.0) # for the SMC
    #R = 24.474  #distance from LMC to galactic center, or from LMC to SMC
    print dm_mass, 'mass SMC'

    rvir = model.mass_to_rvir(dm_mass)
    #factor =  rd.get_normed_radial_abundance(radius / rvir)
    ngreater, N_std = model.ngreater(dm_mass,min_lum=10**3, percentile=False)

    def f(r, theta):
        d = np.sqrt(R**2 + r**2 + 2*r*R*np.cos(theta))
        return ngreater/(4*np.pi*d**2) * rd.get_normed_sat_derivative(d/rvir, rvir)

    def f2(r,theta):
        d = np.sqrt(R**2 + r**2 + 2*r*R*np.cos(theta))
        dr=1.0
        up = rd.get_normed_radial_abundance( (d+dr/2) / rvir)
        down = rd.get_normed_radial_abundance( (d-dr/2) / rvir)
        return (up-down)*ngreater / (4*np.pi*d**2 * dr)
        
    def g(d):
        dr =1.0
        up = rd.get_normed_radial_abundance( (d+dr/2) / rvir)
        down = rd.get_normed_radial_abundance( (d-dr/2) / rvir)
        return up-down

    #term1 = dblquad(lambda r, theta: f(r,theta), 0, np.pi, lambda x: 0, lambda x: 50.0)  # within 50 kpc of the LMC
    #term1 = quad(lambda r: 4*np.pi*r**2*f(r,0), 0, rvir) this does not work
    #term2 = quad(lambda r: 4*np.pi*r**2*f2(r,0), 0, rvir)
    #term3 = dblquad(lambda r, theta: 2*np.pi*r**2 * np.sin(theta)*f(r,theta), 0, np.pi, lambda x: 0, lambda x: rvir) # full MW volume

    term4 = dblquad(lambda r, theta: 2*np.pi*r**2*np.sin(theta)*f(r,theta), 0, np.pi, lambda x: 0, lambda x: 50) # full MW volume when R=0

    #term4 = dblquad(lambda r, theta: 2*np.pi*r**2*np.sin(theta)*f2(r,theta), 0, np.pi, lambda x: 0, lambda x: 50)  # within 50 kpc of the LMC
    #term4 = dblquad(lambda r, theta: 2*np.pi*r**2*np.sin(theta)*f(r,theta), 0, np.pi, lambda x: 0, lambda x: 50)  # within 50 kpc of the LMC

    print term4
    print ngreater, 'ngreater'
    print term4[0]/ngreater
    return term4


def LMC_in_50kpc():  
    #ratio_MW = 0.022158  # value recored in notebook for MOSTER.
    #ratio_SMC = 0.228168 # value for MOSTER
    #ratio_MW = 0.0251748  # value for GK16
    #ratio_SMC = 0.1795   # value for GK16
    ratio_MW = 0.0720138  # value for GK16
    ratio_SMC = 0.38563   # value for GK16
    f_left = .05 # fraction of remaining stars after tidal stripping


    pg = -0.2  # for playing with plot gamma. -0.2 is default
    f, ax = plt.subplots()
    low=2; high=7.5
    min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)  

    # Plot Observed Galaxies 
    galx = dm.load_nearby_gx()
    conf_galx = ~np.array([galx['Name'][i][0]=='*' for i in range(len(galx))])
    near_LMC = galx['dist_LMC'] < 50   # I need real tidal radius of LMC/SCM. i.e. what is sub-substructure according to rockstar?
    near_LMC[(galx['Name']=='LMC')|(galx['Name']=='SMC')] = False     # take out the LMC and SMC from the list
    smass = galx['mstar'][near_LMC]
    conf_smass = galx['mstar'][near_LMC & conf_galx]
    n_greater = np.array([np.sum(smass > min_lums[i]) for i in range(len(min_lums))])
    n_greater_conf = np.array([np.sum(conf_smass > min_lums[i]) for i in range(len(min_lums))])

    ## OLD LABELS
    #ax.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater_conf), label='Confirmed Galaxies',color='black', linestyle='-',lw=linewidth+1)
    #ax.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Confirmed + To Be Confirmed', color='black',lw=linewidth+1,linestyle='--')    

    ax.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')
    


    # Then plot MW background
    model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg)
    halo_mass = 1.4*10**12
    #N_visible_MW, Nstd_MW,p20,p80 = model.ngreater(halo_mass,min_lums,percentile=True) # original method. wrong p20,p80
    samples_MW = model.generate_stellar_mass_samples(halo_mass, factor=ratio_MW)  
    N_visible_MW, Nstd_MW, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_MW)
    ax.plot(min_lums, N_visible_MW, label='MW background', color=BehrooziGreen,lw=linewidth,ls='-')

    # Next plot SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg)
    halo_mass = model.mstar_to_mhalo(SMC_stellar_mass,a=1.0)
    #N_visible_SMC, Nstd_SMC = model.ngreater(halo_mass,min_lums)  # original method. wrong p20, p80
    samples_SMC = model.generate_stellar_mass_samples(halo_mass, factor=ratio_SMC)  
    N_visible_SMC, Nstd_SMC, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_SMC)


    # Plot LMC Alone
    lmc = True  # want to use the lmc SHMF
    for model in [am.GK16_grow(lmcSHMF=lmc,reionization=True, plotgamma=pg)]: #,am.Moster(lmcSHMF=lmc,reionization=True)]:  #,am.Brook(lmcSHMF=lmc,reionization=True)]:   am.GarrisonKimmel(lmcSHMF=lmc,reionization=True),
        print model.label, 'on this model'
        halo_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1.0)
        lmc_rvir = model.mass_to_rvir(halo_mass)
        factor =  rd.get_normed_radial_abundance(50 / lmc_rvir)
        samples_LMC = model.generate_stellar_mass_samples(halo_mass, factor=factor)  
        N_visible, N_std, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_LMC)
        # N_visible,N_std = model.ngreater(halo_mass,min_lums) # former method. wrong N_std

        if isinstance(model, am.GK16_grow):
            ax.plot(min_lums, N_visible, color=BehrooziGreen,lw=linewidth,ls='-.',label='LMC alone')
            # combine all samples. keep the length N_iter the same.
            all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
            N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
            #ax.fill_between(min_lums, p80_full, p20_full, facecolor=model.color, alpha=0.25)
            ax.fill_between(min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.25)
            # Want to plot SMC after LMC for ordering
            ax.plot(min_lums, N_visible_SMC, label='SMC alone', color=BehrooziGreen,lw=linewidth,ls=':')
            #ax.fill_between(min_lums*f_left, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.25)

        ax.plot(min_lums, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls=model.ls,label='LMC+SMC+MW (Baseline)') #model.label)    
        ax.plot(min_lums * f_left, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='--',label='LMC+SMC+MW (Tidally Stripped)') #model.label)    



    #### Change Reionization Model
    vmax_ach = 9.48535156   # multiply by 0.6 for GK16 with z=13.  For z=9, Moster use 0.8, GK16 use 0.85
    vmax_filt = 23.54248047  # has little effect
    z = 11  # choices of 14,13,12,11,9

    # Then plot MW background
    model = am.GK16_grow(hostSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)
    halo_mass = 1.4*10**12
    samples_MW = model.generate_stellar_mass_samples(halo_mass, factor=ratio_MW)  
    N_visible_MW, Nstd_MW, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_MW)
    #ax.plot(min_lums, N_visible_MW, label='MW background', color=BehrooziGreen,lw=linewidth,ls='-')

    # SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)
    halo_mass = model.mstar_to_mhalo(SMC_stellar_mass,a=1.0)
    #N_visible_SMC, Nstd_SMC = model.ngreater(halo_mass,min_lums)  # original method. wrong p20, p80
    samples_SMC = model.generate_stellar_mass_samples(halo_mass, factor=ratio_SMC)  
    N_visible_SMC, Nstd_SMC, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_SMC)
    #ax.plot(min_lums, N_visible_SMC, label='SMC Weak Reion', color='grey',lw=linewidth,ls=':')


    lmc = True  # want to use the lmc SHMF
    for model in [am.GK16_grow(lmcSHMF=lmc,catreion=True,z=z,reionization=True, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)]:
        # am.Moster(lmcSHMF=lmc,reionization=True,catreion=True,z=9, vmax_ach=vmax_ach, vmax_filt=vmax_filt)
        print model.label, 'on this model'
        halo_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1.0)
        lmc_rvir = model.mass_to_rvir(halo_mass)
        factor =  rd.get_normed_radial_abundance(50 / lmc_rvir)
        samples_LMC = model.generate_stellar_mass_samples(halo_mass, factor=factor)  
        N_visible, N_std, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_LMC)
        #ax.plot(min_lums, N_visible,label=model.label,color=model.color,lw=linewidth,ls='--')
        #ax.plot(min_lums, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='--',label='LMC+SMC+MW (Late Reion)')
        #ax.plot(min_lums*0.5, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='--',label='LMC+SMC+MW (Late Reion)')
        ax.plot(min_lums*f_left, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='--',label='LMC+SMC+MW (Late Reion)')

        if isinstance(model, am.GK16_grow):
            # combine all samples. keep the length N_iter the same.
            all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
            N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
            #ax.fill_between(min_lums, p80_full, p20_full, facecolor=model.color, alpha=0.25)
            ax.fill_between(min_lums*f_left, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.25)

    ax.set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=label_font)
    ax.set_xlabel('$\mathrm{M_*} \, [\mathrm{M_\odot}]$', fontsize=label_font)
    plt.xscale('log')
    ax.legend(loc='upper right',frameon=False,fontsize=legend_size-4)

    plt.xlim((min_lums[0],10**5.5))
    plt.ylim((0,15))
    ax.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.subplots_adjust(hspace = 0.1)
    plt.gcf().subplots_adjust(bottom=0.13)
    extra=''
    plt.savefig('LMCPlots/LMC_in_50kpc'+extra+'.png')
    plt.close()


