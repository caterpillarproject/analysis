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
from matplotlib.font_manager import FontProperties


LMC_stellar_mass = 2.6*10**9
SMC_stellar_mass = 7.1*10**8


def Baseline(verbose=False):  
    ratio_MW = 0.0720138  # value for GK16
    ratio_SMC = 0.38563   # value for GK16
    f_left = 1 # fraction of remaining stars after tidal stripping

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
    ax.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')
    

    # Then plot MW background
    model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg)
    #model = am.Brook(hostSHMF=True,reionization=True)

    halo_mass = 1.4*10**12
    samples_MW = model.generate_stellar_mass_samples(halo_mass, factor=ratio_MW)  
    N_visible_MW, Nstd_MW, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_MW)
    ax.plot(min_lums, N_visible_MW, label='MW background', color=BehrooziGreen,lw=linewidth,ls='-')

    # Next plot SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg)
    #model = am.Brook(lmcSHMF=True,reionization=True)

    halo_mass = model.mstar_to_mhalo(SMC_stellar_mass,a=1.0)
    samples_SMC = model.generate_stellar_mass_samples(halo_mass, factor=ratio_SMC)  
    N_visible_SMC, Nstd_SMC, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_SMC)


    # Plot LMC Alone
    lmc = True  # want to use the lmc SHMF
    for model in [am.GK16_grow(lmcSHMF=lmc,reionization=True, plotgamma=pg)]:  #,am.Brook(lmcSHMF=lmc,reionization=True)]:   am.GarrisonKimmel(lmcSHMF=lmc,reionization=True),
    #for model in [am.Brook(lmcSHMF=lmc,reionization=True)]: 
        print model.label, 'on this model'
        halo_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1.0)
        lmc_rvir = model.mass_to_rvir(halo_mass)
        factor =  rd.get_normed_radial_abundance(50 / lmc_rvir)
        samples_LMC = model.generate_stellar_mass_samples(halo_mass, factor=factor)  
        N_visible, N_std, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_LMC)
        # N_visible,N_std = model.ngreater(halo_mass,min_lums) # former method. wrong N_std

        if isinstance(model, am.GK16_grow):
        #if isinstance(model, am.Brook):
            ax.plot(min_lums, N_visible, color=BehrooziGreen,lw=linewidth,ls='-.',label='LMC alone')
            # combine all samples. keep the length N_iter the same.
            all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
            N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
            ax.fill_between(min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.25)
            # Want to plot SMC after LMC for ordering
            ax.plot(min_lums, N_visible_SMC, label='SMC alone', color=BehrooziGreen,lw=linewidth,ls=':')
            #ax.fill_between(min_lums*f_left, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.25)

        ax.plot(min_lums, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='--',label='LMC+SMC+MW') #model.label)    
        #ax.plot(min_lums * f_left, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='--',label='LMC+SMC+MW (Tidally Stripped)') #model.label)    


    if verbose:
        print min_lums, 'min lums'
        print N_visible ,'LMC'
        print N_visible_SMC, 'SMC'
        print N_visible_MW, 'MW'
        print N_visible_full, 'total'
        print N_visible_full+N_std_full, 'up 1 sig'
        print N_visible_full-N_std_full, 'down 1 sig'
        plt.close()
        return min_lums, N_visible, N_visible_SMC, N_visible_MW, N_visible_full, N_visible_full+N_std_full, N_visible_full-N_std_full, all_samples


    font0 = FontProperties()
    font = font0.copy()
    font.set_weight('bold')
    ax.text(.6,.56,'Original Model',transform=ax.transAxes,fontsize=legend_size, fontproperties=font)

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
    plt.savefig('LMCPlots/LMC_in_50kpc1'+extra+'.png')
    plt.close()




def TidalStripped():  
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
    ax.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')
    

    # Then plot MW background
    model = am.GK16_grow(hostSHMF=True,reionization=True)
    halo_mass = 1.4*10**12
    samples_MW = model.generate_stellar_mass_samples(halo_mass, factor=ratio_MW)  
    N_visible_MW, Nstd_MW, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_MW)
    ax.plot(min_lums*f_left, N_visible_MW, label='MW background', color=BehrooziGreen,lw=linewidth,ls='-')

    # SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True)
    halo_mass = model.mstar_to_mhalo(SMC_stellar_mass,a=1.0)
    #N_visible_SMC, Nstd_SMC = model.ngreater(halo_mass,min_lums)  # original method. wrong p20, p80
    samples_SMC = model.generate_stellar_mass_samples(halo_mass, factor=ratio_SMC)  
    N_visible_SMC, Nstd_SMC, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_SMC)
    ax.plot(min_lums*f_left, N_visible_SMC, label='SMC alone', color='grey',lw=linewidth,ls=':')


    lmc = True  # want to use the lmc SHMF
    for model in [am.GK16_grow(lmcSHMF=lmc,reionization=True)]:
        # am.Moster(lmcSHMF=lmc,reionization=True,catreion=True,z=9, vmax_ach=vmax_ach, vmax_filt=vmax_filt)
        print model.label, 'on this model'
        halo_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1.0)
        lmc_rvir = model.mass_to_rvir(halo_mass)
        factor =  rd.get_normed_radial_abundance(50 / lmc_rvir)
        samples_LMC = model.generate_stellar_mass_samples(halo_mass, factor=factor)  
        N_visible, N_std, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_LMC)

        if isinstance(model, am.GK16_grow):
            ax.plot(min_lums*f_left, N_visible, color=BehrooziGreen,lw=linewidth,ls='-.',label='LMC alone')
            # combine all samples. keep the length N_iter the same.
            all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
            N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
            ax.fill_between(min_lums*f_left, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.25)

        ax.plot(min_lums*f_left, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='--',label='LMC+SMC+MW')


    #print N_visible ,'LMC'
    #print N_visible_SMC, 'SMC'
    #print N_visible_MW, 'MW'

    font0 = FontProperties()
    font = font0.copy()
    font.set_weight('bold')
    ax.text(.48,.56,'95% Tidally Stripped',transform=ax.transAxes,fontsize=legend_size, fontproperties=font)

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
    plt.savefig('LMCPlots/LMC_in_50kpc2'+extra+'.png')
    plt.close()



def TidalStripped_Reionization():  
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
    ax.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')
    

    #### Change Reionization Model
    vmax_ach = 9.48535156   # multiply by 0.6 for GK16 with z=13.  For z=9, Moster use 0.8, GK16 use 0.85
    vmax_filt = 23.54248047  # has little effect
    z = 11  # choices of 14,13,12,11,9

    # Then plot MW background
    #model = am.GK16_grow(hostSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)
    model = am.Brook(hostSHMF=True,reionization=False, catreion=True,z=9, vmax_ach=vmax_ach*0.2, vmax_filt=vmax_filt)
    halo_mass = 1.4*10**12
    samples_MW = model.generate_stellar_mass_samples(halo_mass, factor=ratio_MW)  
    N_visible_MW, Nstd_MW, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_MW)
    ax.plot(min_lums*f_left, N_visible_MW, label='MW background', color=BehrooziGreen,lw=linewidth,ls='-')

    # SMC Background
    #model = am.GK16_grow(lmcSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)
    model = am.Brook(lmcSHMF=True,reionization=False, catreion=True,z=9, vmax_ach=vmax_ach*.2, vmax_filt=vmax_filt)
    halo_mass = model.mstar_to_mhalo(SMC_stellar_mass,a=1.0)
    #N_visible_SMC, Nstd_SMC = model.ngreater(halo_mass,min_lums)  # original method. wrong p20, p80
    samples_SMC = model.generate_stellar_mass_samples(halo_mass, factor=ratio_SMC)  
    N_visible_SMC, Nstd_SMC, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_SMC)
    ax.plot(min_lums*f_left, N_visible_SMC, label='SMC alone', color='grey',lw=linewidth,ls=':')


    lmc = True  # want to use the lmc SHMF
    #for model in [am.GK16_grow(lmcSHMF=lmc,catreion=True,z=z,reionization=True, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)]:
    for model in [am.Brook(lmcSHMF=lmc,catreion=True,z=9,reionization=False, vmax_ach=vmax_ach*.2, vmax_filt=vmax_filt)]:
        # am.Moster(lmcSHMF=lmc,reionization=True,catreion=True,z=9, vmax_ach=vmax_ach, vmax_filt=vmax_filt)
        print model.label, 'on this model'
        halo_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1.0)
        lmc_rvir = model.mass_to_rvir(halo_mass)
        factor =  rd.get_normed_radial_abundance(50 / lmc_rvir)
        samples_LMC = model.generate_stellar_mass_samples(halo_mass, factor=factor)  
        N_visible, N_std, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_LMC)

        if isinstance(model, am.GK16_grow) or isinstance(model, am.Brook):
            ax.plot(min_lums*f_left, N_visible, color=BehrooziGreen,lw=linewidth,ls='-.',label='LMC alone')
            # combine all samples. keep the length N_iter the same.
            all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
            N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
            ax.fill_between(min_lums*f_left, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.25)

        ax.plot(min_lums*f_left, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='-',label='LMC+SMC+MW')


    #print N_visible ,'LMC'
    #print N_visible_SMC, 'SMC'
    #print N_visible_MW, 'MW'
    #print N_visible + N_visible_SMC + N_visible_MW



    font0 = FontProperties()
    font = font0.copy()
    font.set_weight('bold')
    ax.text(.48,.51,'95% Tidally Stripped \n + Later Reionization',transform=ax.transAxes,fontsize=legend_size, fontproperties=font)


    ax.set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=label_font)
    ax.set_xlabel('$\mathrm{M_*} \, [\mathrm{M_\odot}]$', fontsize=label_font)
    plt.xscale('log')
    ax.legend(loc='upper right',frameon=False,fontsize=legend_size-4)

    plt.xlim((min_lums[0],10**5.5))
    plt.ylim((0,15))
    ax.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.subplots_adjust(hspace = 0.1)
    plt.gcf().subplots_adjust(bottom=0.13)
    extra='Brook'
    plt.savefig('LMCPlots/LMC_in_50kpc3'+extra+'.png')
    plt.close()









def Steeper(verbose=False):  
    ratio_MW = 0.0720138  # value for GK16
    ratio_SMC = 0.38563   # value for GK16

    pg = -1.0  # for playing with plot gamma. -0.2 is default
    print pg, 'gamma for scatter'
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
    ax.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')
    

    # Then plot MW background
    model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg)
    halo_mass = 1.4*10**12
    samples_MW = model.generate_stellar_mass_samples(halo_mass, factor=ratio_MW)  
    N_visible_MW, Nstd_MW, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_MW)
    ax.plot(min_lums, N_visible_MW, label='MW background', color=BehrooziGreen,lw=linewidth,ls='-')

    # Next plot SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg)
    halo_mass = model.mstar_to_mhalo(SMC_stellar_mass,a=1.0)
    samples_SMC = model.generate_stellar_mass_samples(halo_mass, factor=ratio_SMC)  
    N_visible_SMC, Nstd_SMC, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_SMC)


    # Plot LMC Alone
    lmc = True  # want to use the lmc SHMF
    for model in [am.GK16_grow(lmcSHMF=lmc,reionization=True, plotgamma=pg)]:  #,am.Brook(lmcSHMF=lmc,reionization=True)]:   am.GarrisonKimmel(lmcSHMF=lmc,reionization=True),
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
            ax.fill_between(min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.25)
            # Want to plot SMC after LMC for ordering
            ax.plot(min_lums, N_visible_SMC, label='SMC alone', color=BehrooziGreen,lw=linewidth,ls=':')
            #ax.fill_between(min_lums*f_left, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.25)

        ax.plot(min_lums, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='--',label='LMC+SMC+MW') #model.label)    
        #ax.plot(min_lums * f_left, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='--',label='LMC+SMC+MW (Tidally Stripped)') #model.label)    


    if verbose:
        print min_lums, 'min lums'
        print N_visible ,'LMC'
        print N_visible_SMC, 'SMC'
        print N_visible_MW, 'MW'
        print N_visible_full, 'total'
        print N_visible_full+N_std_full, 'up 1 sig'
        print N_visible_full-N_std_full, 'down 1 sig'
        plt.close()
        return min_lums, N_visible, N_visible_SMC, N_visible_MW, N_visible_full, N_visible_full+N_std_full, N_visible_full-N_std_full, all_samples


    font0 = FontProperties()
    font = font0.copy()
    font.set_weight('bold')
    ax.text(.6,.56,'Original Model',transform=ax.transAxes,fontsize=legend_size, fontproperties=font)

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
    plt.savefig('LMCPlots/LMC_in_50kpc4'+extra+'.png')
    plt.close()







def SteeperTidal():  
    ratio_MW = 0.0720138  # value for GK16
    ratio_SMC = 0.38563   # value for GK16
    f_left = 0.3 # fraction of remaining stars after tidal stripping

    pg = -0.7  # for playing with plot gamma. -0.2 is default
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
    ax.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')
    

    # Then plot MW background
    model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg)
    halo_mass = 1.4*10**12
    samples_MW = model.generate_stellar_mass_samples(halo_mass, factor=ratio_MW)  
    N_visible_MW, Nstd_MW, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_MW)
    ax.plot(min_lums*f_left, N_visible_MW, label='MW background', color=BehrooziGreen,lw=linewidth,ls='-')

    # SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg)
    halo_mass = model.mstar_to_mhalo(SMC_stellar_mass,a=1.0)
    #N_visible_SMC, Nstd_SMC = model.ngreater(halo_mass,min_lums)  # original method. wrong p20, p80
    samples_SMC = model.generate_stellar_mass_samples(halo_mass, factor=ratio_SMC)  
    N_visible_SMC, Nstd_SMC, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_SMC)
    ax.plot(min_lums*f_left, N_visible_SMC, label='SMC alone', color='grey',lw=linewidth,ls=':')


    lmc = True  # want to use the lmc SHMF
    for model in [am.GK16_grow(lmcSHMF=lmc,reionization=True, plotgamma=pg)]:
        # am.Moster(lmcSHMF=lmc,reionization=True,catreion=True,z=9, vmax_ach=vmax_ach, vmax_filt=vmax_filt)
        print model.label, 'on this model'
        halo_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1.0)
        lmc_rvir = model.mass_to_rvir(halo_mass)
        factor =  rd.get_normed_radial_abundance(50 / lmc_rvir)
        samples_LMC = model.generate_stellar_mass_samples(halo_mass, factor=factor)  
        N_visible, N_std, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_LMC)

        if isinstance(model, am.GK16_grow):
            ax.plot(min_lums*f_left, N_visible, color=BehrooziGreen,lw=linewidth,ls='-.',label='LMC alone')
            # combine all samples. keep the length N_iter the same.
            all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
            N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
            ax.fill_between(min_lums*f_left, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.25)

        ax.plot(min_lums*f_left, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='--',label='LMC+SMC+MW')


    #print N_visible ,'LMC'
    #print N_visible_SMC, 'SMC'
    #print N_visible_MW, 'MW'

    font0 = FontProperties()
    font = font0.copy()
    font.set_weight('bold')
    ax.text(.48,.56,'70% Tidally Stripped',transform=ax.transAxes,fontsize=legend_size, fontproperties=font)

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
    plt.savefig('LMCPlots/LMC_in_50kpc5'+extra+'.png')
    plt.close()



def SteeperTidalReionization():  
    ratio_MW = 0.0720138  # value for GK16
    ratio_SMC = 0.38563   # value for GK16
    f_left = .30 # fraction of remaining stars after tidal stripping

    pg = -0.7  # for playing with plot gamma. -0.2 is default
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
    ax.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')
    

    #### Change Reionization Model
    vmax_ach = 9.48535156 * 0.85   # multiply by 0.6 for GK16 with z=13.  For z=9, Moster use 0.8, GK16 use 0.85
    vmax_filt = 23.54248047  # has little effect
    z = 9  # choices of 14,13,12,11,9

    # Then plot MW background
    model = am.GK16_grow(hostSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)
    #model = am.Brook(hostSHMF=True,reionization=False, catreion=True,z=9, vmax_ach=vmax_ach*0.2, vmax_filt=vmax_filt)
    halo_mass = 1.4*10**12
    samples_MW = model.generate_stellar_mass_samples(halo_mass, factor=ratio_MW)  
    N_visible_MW, Nstd_MW, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_MW)
    ax.plot(min_lums*f_left, N_visible_MW, label='MW background', color=BehrooziGreen,lw=linewidth,ls='-')

    # SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)
    #model = am.Brook(lmcSHMF=True,reionization=False, catreion=True,z=9, vmax_ach=vmax_ach*.2, vmax_filt=vmax_filt)
    halo_mass = model.mstar_to_mhalo(SMC_stellar_mass,a=1.0)
    #N_visible_SMC, Nstd_SMC = model.ngreater(halo_mass,min_lums)  # original method. wrong p20, p80
    samples_SMC = model.generate_stellar_mass_samples(halo_mass, factor=ratio_SMC)  
    N_visible_SMC, Nstd_SMC, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_SMC)
    ax.plot(min_lums*f_left, N_visible_SMC, label='SMC alone', color='grey',lw=linewidth,ls=':')


    lmc = True  # want to use the lmc SHMF
    for model in [am.GK16_grow(lmcSHMF=lmc,catreion=True,z=z,reionization=True, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)]:
    #for model in [am.Brook(lmcSHMF=lmc,catreion=True,z=9,reionization=False, vmax_ach=vmax_ach*.2, vmax_filt=vmax_filt)]:
        # am.Moster(lmcSHMF=lmc,reionization=True,catreion=True,z=9, vmax_ach=vmax_ach, vmax_filt=vmax_filt)
        print model.label, 'on this model'
        halo_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1.0)
        lmc_rvir = model.mass_to_rvir(halo_mass)
        factor =  rd.get_normed_radial_abundance(50 / lmc_rvir)
        samples_LMC = model.generate_stellar_mass_samples(halo_mass, factor=factor)  
        N_visible, N_std, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_LMC)

        if isinstance(model, am.GK16_grow) or isinstance(model, am.Brook):
            ax.plot(min_lums*f_left, N_visible, color=BehrooziGreen,lw=linewidth,ls='-.',label='LMC alone')
            # combine all samples. keep the length N_iter the same.
            all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
            N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
            ax.fill_between(min_lums*f_left, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.25)

        ax.plot(min_lums*f_left, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='-',label='LMC+SMC+MW')


    #print N_visible ,'LMC'
    #print N_visible_SMC, 'SMC'
    #print N_visible_MW, 'MW'
    #print N_visible + N_visible_SMC + N_visible_MW



    font0 = FontProperties()
    font = font0.copy()
    font.set_weight('bold')
    ax.text(.48,.51,'70% Tidally Stripped \n + Later Reionization',transform=ax.transAxes,fontsize=legend_size, fontproperties=font)


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
    plt.savefig('LMCPlots/LMC_in_50kpc6'+extra+'.png')
    plt.close()




################# AFTER TALKING TO ANNIKA  ##############################

# try to rule out the Behroozi and Brook models.
# Also identify that there is a problem!




def getMW(model, min_lums):
    ratio_MW = 0.0720138  # value for GK16
    # Get MW Component
    halo_mass = 1.4*10**12
    samples_MW = model.generate_stellar_mass_samples(halo_mass, factor=ratio_MW)  
    N_visible_MW, Nstd_MW, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_MW)
    return N_visible_MW, Nstd_MW, samples_MW

def getSMC(model,min_lums):
    ratio_SMC = 0.38563   # value for GK16
    halo_mass = model.mstar_to_mhalo(SMC_stellar_mass,a=1.0)
    samples_SMC = model.generate_stellar_mass_samples(halo_mass, factor=ratio_SMC)  
    N_visible_SMC, Nstd_SMC, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_SMC)
    return N_visible_SMC, Nstd_SMC, samples_SMC

def getLMC(model,min_lums):
        halo_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1.0)
        lmc_rvir = model.mass_to_rvir(halo_mass)
        factor =  rd.get_normed_radial_abundance(50 / lmc_rvir)
        samples_LMC = model.generate_stellar_mass_samples(halo_mass, factor=factor)  
        N_visible, N_std, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_LMC)
        return N_visible, N_std, samples_LMC
        


def Original(verbose=False):  
    f_left = 1 # fraction of remaining stars after tidal stripping
    pg = -0.2  # for playing with plot gamma. -0.2 is default
    f, ax = plt.subplots()
    low=2; high=7.5
    min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)  

    # Plot Observed Galaxies 
    galx = dm.load_nearby_gx()
    #conf_galx = ~np.array([galx['Name'][i][0]=='*' for i in range(len(galx))])
    near_LMC = galx['dist_LMC'] < 50   # I need real tidal radius of LMC/SCM. i.e. what is sub-substructure according to rockstar?
    near_LMC[(galx['Name']=='LMC')|(galx['Name']=='SMC')] = False     # take out the LMC and SMC from the list
    print galx['Name'][near_LMC]
    smass = galx['mstar'][near_LMC]
    #conf_smass = galx['mstar'][near_LMC & conf_galx]
    n_greater = np.array([np.sum(smass > min_lums[i]) for i in range(len(min_lums))])
    #n_greater_conf = np.array([np.sum(conf_smass > min_lums[i]) for i in range(len(min_lums))])
    ax.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')
    

    # Then plot MW background
    model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg)
    N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
    ax.plot(min_lums, N_visible_MW, label='MW background', color='grey',lw=linewidth,ls='-')

    model = am.Brook(hostSHMF=True,reionization=True)
    Brook_N_MW, Brook_std_MW, Brook_samples_MW = getMW(model,min_lums)
    model = am.Behroozi(hostSHMF=True,reionization=True)
    Beh_N_MW, Beh_std_MW, Beh_samples_MW = getMW(model,min_lums)
    #model = am.Brook(hostSHMF=True,reionization=False)
    #BrookNR_N_MW, BrookNR_std_MW, BrookNR_samples_MW = getMW(model,min_lums)


    # Next plot SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg)
    N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)

    model = am.Brook(lmcSHMF=True,reionization=True)
    Brook_N_SMC, Brook_std_SMC, Brook_samples_SMC = getSMC(model,min_lums)
    model = am.Behroozi(lmcSHMF=True,reionization=True)
    Beh_N_SMC, Beh_std_SMC, Beh_samples_SMC = getSMC(model,min_lums)
    #model = am.Brook(lmcSHMF=True,reionization=False)
    #BrookNR_N_SMC, BrookNR_std_SMC, BrookNR_samples_SMC = getSMC(model,min_lums)


    # Plot LMC Alone
    lmc = True  # want to use the lmc SHMF
    for model in [am.GK16_grow(lmcSHMF=lmc,reionization=True, plotgamma=pg), am.Behroozi(lmcSHMF=lmc,reionization=True), am.Brook(lmcSHMF=lmc,reionization=True)]:   #  am.Brook(lmcSHMF=lmc,reionization=False),
        print model.label, 'on this model'
        N_visible, N_std, samples_LMC = getLMC(model,min_lums)

        if isinstance(model, am.GK16_grow):
            ax.plot(min_lums, N_visible, color='grey',lw=linewidth,ls='-.',label='LMC alone')
            # combine all samples. keep the length N_iter the same.
            all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
            N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
            ax.fill_between(min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.25)
            # Want to plot SMC after LMC for ordering
            ax.plot(min_lums, N_visible_SMC, label='SMC alone', color='grey',lw=linewidth,ls=':')
            #ax.fill_between(min_lums*f_left, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.25)
            ax.plot(min_lums, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='-',label=model.label)   #'LMC+SMC+MW') 

        if isinstance(model, am.Brook):
            if model.reionization:
                ax.plot(min_lums, N_visible + Brook_N_SMC + Brook_N_MW, color=model.color,lw=linewidth+2,ls='-',label=model.label)
            else:
                ax.plot(min_lums, N_visible + BrookNR_N_SMC + BrookNR_N_MW, color=model.color,lw=linewidth+2,ls='--',label=model.label)
        if isinstance(model, am.Behroozi):
            ax.plot(min_lums, N_visible + Beh_N_SMC + Beh_N_MW, color=model.color,lw=linewidth+2,ls='-',label=model.label)


        #ax.plot(min_lums * f_left, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='--',label='LMC+SMC+MW (Tidally Stripped)') #model.label)    


    if verbose:
        print min_lums, 'min lums'
        print N_visible ,'LMC'
        print N_visible_SMC, 'SMC'
        print N_visible_MW, 'MW'
        print N_visible_full, 'total'
        print N_visible_full+N_std_full, 'up 1 sig'
        print N_visible_full-N_std_full, 'down 1 sig'
        plt.close()
        return min_lums, N_visible, N_visible_SMC, N_visible_MW, N_visible_full, N_visible_full+N_std_full, N_visible_full-N_std_full, all_samples


    font0 = FontProperties()
    font = font0.copy()
    font.set_weight('bold')
    #ax.text(.6,.56,'Original Model',transform=ax.transAxes,fontsize=legend_size, fontproperties=font)

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
    plt.savefig('LMCPlots/LMC_in_50kpcBehrooziBrook'+extra+'.png')
    plt.close()









def Grid(verbose=False):
    # 8 x 6 is default
    f, ax = plt.subplots(3, 3, sharex='col',sharey='row', figsize=(16,13))
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    font0 = FontProperties()
    font = font0.copy()
    #font.set_weight('bold')
    vmax_ach = 9.48535156  # multiply by 0.6 for GK16 with z=13.  For z=9, Moster use 0.8, GK16 use 0.85
    vmax_filt = 23.54248047  # has little effect
    redshift = [[14,12,7],
                [14,11,7],
                [13,10,6]]


    #f_left = 1 # fraction of remaining stars after tidal stripping
    #pg = -0.2  # for playing with plot gamma. -0.2 is default
    z = 13.3

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


    for row, f_left in zip(range(3), [1,.7,.5] ):
        for col, pg in zip(range(3), [-.2,-.5,-1.0] ):
            # convert pg to slope of GK line
            z = redshift[row][col]
 
            ax[row][col].plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')
            # Then plot MW background
            #model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg)
            model = am.GK16_grow(hostSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)
            N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
            ax[row][col].plot(min_lums*f_left, N_visible_MW, label='MW background', color='grey',lw=linewidth,ls='-')

    
            # Next plot SMC Background
            #model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg)
            model = am.GK16_grow(lmcSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)
            N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)

            # Plot LMC Alone
            lmc = True  # want to use the lmc SHMF
            for model in [am.GK16_grow(lmcSHMF=True,reionization=True,catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)]:
                #am.GK16_grow(lmcSHMF=lmc,reionization=True, plotgamma=pg)]:   #  am.Brook(lmcSHMF=lmc,reionization=False),
                print model.label, 'on this model'
                N_visible, N_std, samples_LMC = getLMC(model,min_lums)

                if isinstance(model, am.GK16_grow):
                    ax[row][col].plot(min_lums*f_left, N_visible, color='grey',lw=linewidth,ls='-.',label='LMC alone')
                    # combine all samples. keep the length N_iter the same.
                    all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
                    N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
                    ax[row][col].fill_between(min_lums*f_left, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.25)
                    ax[row][col].plot(min_lums*f_left, N_visible_SMC, label='SMC alone', color='grey',lw=linewidth,ls=':')
                    ax[row][col].plot(min_lums*f_left, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='-',label='LMC+SMC+MW') #model.label)

            # plotting things
            ax[row][col].set_xscale('log')
            ax[row][col].set_xlim((min_lums[0],10**5.5))
            ax[row][col].set_ylim((0,15))
            ax[row][col].tick_params(axis='both', which='major', labelsize=tick_size)
            ax[row][col].text(.5,.90,'AM slope = '+str(model.alphaAM0),transform=ax[row][col].transAxes,fontsize=legend_size-1, fontproperties=font)
            ax[row][col].text(.7,.8,'$\mathrm{f_{strip}} = $'+str(1-f_left),transform=ax[row][col].transAxes,fontsize=legend_size-1, fontproperties=font)
            ax[row][col].text(.65,.7,'$\mathrm{Z_{reion}} = $'+str(z)+'.3',transform=ax[row][col].transAxes,fontsize=legend_size-1, fontproperties=font)   
    
    for i in range(3):
        ax[i][0].set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=label_font+1)
        ax[2][i].set_xlabel('$\mathrm{M_*} \, [\mathrm{M_\odot}]$', fontsize=label_font+1)

    #ax[0][2].legend(loc='upper right',frameon=False,fontsize=legend_size-2)


    plt.subplots_adjust(hspace = 0.1)
    plt.gcf().subplots_adjust(bottom=0.13)
    extra=''
    plt.savefig('LMCPlots/LMC_in_50kpcGrid'+extra+'.png')
    plt.close()
