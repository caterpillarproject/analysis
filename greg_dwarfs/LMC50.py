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
vmax_ach = 9.48535156   # multiply by 0.6 for GK16 with z=13.  For z=9, Moster use 0.8, GK16 use 0.85
vmax_filt = 23.54248047  # has little effect


################# AFTER TALKING TO ANNIKA  ##############################

def get_ratio(model, galaxy):
    if galaxy == 'MW':
        if model.label == 'Moster + Reion':
            return 0.0698049317676
            print 'in moster'
        if model.label == 'Brook + Reion':
            return 0.0648588855887
        if model.label == 'GK16 + Reion':
            if model.plotgamma == -0.2:
                return 0.0766766197663
            if model.plotgamma == -0.5:
                return 0.0766766198533
            if model.plotgamma == -1.0:
                return 0.0766766198627
        if model.label == 'Behroozi + Reion':
            return 0.0766766197664 # Behroozi
        else:
            print 'in else - should only be here for bent plot'
            return  0.0766766198533
    if galaxy == 'SMC':
        if model.label == 'Moster + Reion':
            return 0.462038583656
            print 'in moster'
        if model.label == 'Brook + Reion':
            return 0.455996567841
        if model.label == 'GK16 + Reion':
            if model.plotgamma == -0.2:
                return 0.399382697282
            if model.plotgamma == -0.5:
                return 0.384700258627
            if model.plotgamma == -1.0:
                return 0.363079218041
        if model.label == 'Behroozi + Reion':
            return 0.429648136194 #Behroozi
        else:
            print 'in else - should only be here for bent plot'
            return 0.384700258627
        

# try to rule out the Behroozi and Brook models.
# Also identify that there is a problem!
def getMW(model, min_lums):
    ratio_MW = get_ratio(model, 'MW')   # 0.0720138 
    # Get MW Component
    halo_mass = 1.4*10**12
    samples_MW = model.generate_stellar_mass_samples(halo_mass, factor=ratio_MW)  
    N_visible_MW, Nstd_MW, p20, p80 = ambase.Nvisible_From_Ngreater(min_lums, samples_MW)
    return N_visible_MW, Nstd_MW, samples_MW

def getSMC(model,min_lums):
    ratio_SMC = get_ratio(model, 'SMC')  #  0.38563   # value for GK16
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
        
def get_observed(min_lums):
    galx = dm.load_nearby_gx()
    near_LMC = galx['dist_LMC'] < 50   # I need real tidal radius of LMC/SCM. i.e. what is sub-substructure according to rockstar?
    near_LMC[(galx['Name']=='LMC')|(galx['Name']=='SMC')] = False     # take out the LMC and SMC from the list
    #print galx['Name'][near_LMC]
    smass = galx['mstar'][near_LMC]
    n_greater = np.array([np.sum(smass > min_lums[i]) for i in range(len(min_lums))])
    return n_greater


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
            ax.fill_between(min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.2)
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
    plt.savefig('LMCPlots/LMC_in_50kpcBehrooziBrook'+extra+'.pdf')  #'.png')
    plt.close()





# Plot the GK16 original, and with 30% stripped and 70% stripped.
def TidallyStripped(verbose=False):  
    f_left = 0.02 # fraction of remaining stars after tidal stripping
    pg = -0.2  # for playing with plot gamma. -0.2 is default
    w=8; h=6
    f, (ax, ax2) = plt.subplots(nrows=2,ncols=1,sharex=True, figsize=(w,h*1.6))
    low=2; high=7.5
    min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)  

    # Plot Observed Galaxies
    n_greater = get_observed(min_lums)
    ax.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')
    ax2.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')

    # Then plot MW background
    model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg)
    N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
    #ax.plot(min_lums, N_visible_MW, label='MW background', color='grey',lw=linewidth,ls='-')

    # Next plot SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg)
    N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)

    # Plot LMC Alone
    lmc = True  # want to use the lmc SHMF
    for model in [am.GK16_grow(lmcSHMF=lmc,reionization=True, plotgamma=pg)]:
        print model.label, 'on this model'
        N_visible, N_std, samples_LMC = getLMC(model,min_lums)

        if isinstance(model, am.GK16_grow):
            #ax.plot(min_lums, N_visible, color='grey',lw=linewidth,ls='-.',label='LMC alone')
            # combine all samples. keep the length N_iter the same.
            all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
            N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
            ax.fill_between(f_left*min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.2)
            # Want to plot SMC after LMC for ordering
            #ax.plot(min_lums, N_visible_SMC, label='SMC alone', color='grey',lw=linewidth,ls=':')
            #ax.fill_between(min_lums*f_left, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.2)
            ax.plot(min_lums, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='-',label=model.label)
            ax.plot(0.7*min_lums, (N_visible + N_visible_SMC + N_visible_MW), color=model.color,lw=linewidth+2,ls='--',label='$30\%$ stripped')
            ax.plot(f_left*min_lums, (N_visible + N_visible_SMC + N_visible_MW), color=model.color,lw=linewidth+2,ls='-.',label='$98\%$ stripped')

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

    ##### Plot the 2nd panel
    # Then plot MW background
    z = 8
    model = am.GK16_grow(hostSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)
    N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
    ax2.plot(f_left*min_lums, N_visible_MW, label='MW background', color='grey',lw=linewidth,ls='-')
    
    # Next plot SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)
    N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)
    
    # Plot LMC Alone
    N_visible, N_std, samples_LMC = getLMC(model,min_lums)
    ax2.plot(f_left*min_lums, N_visible, color='grey',lw=linewidth,ls='-.',label='LMC alone')

    # combine all samples. keep the length N_iter the same.
    all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
    N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
    ax2.fill_between(f_left*min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.2)
    # Want to plot SMC after LMC for ordering
    ax2.plot(f_left*min_lums, N_visible_SMC, label='SMC alone', color='grey',lw=linewidth,ls=':')
    ax2.plot(f_left*min_lums, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='-',label='Adjusted\nGK16 + Reion')  #model.label)

    font0 = FontProperties()
    font = font0.copy()
    #font.set_weight('bold')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.2)
    ax2.text(.25,.8,'$98\%$ stripped \n $z_{\mathrm{reion}} = 8.3$',transform=ax2.transAxes,fontsize=legend_size-1, fontproperties=font, bbox=props)
    ##### End plot the 2nd panel


    font0 = FontProperties()
    font = font0.copy()
    font.set_weight('bold')
    #ax.text(.6,.56,'Original Model',transform=ax.transAxes,fontsize=legend_size, fontproperties=font)

    ax.set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=label_font)
    ax2.set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=label_font)
    ax2.set_xlabel('$\mathrm{M_*} \, [\mathrm{M_\odot}]$', fontsize=label_font)
    ax.set_xscale('log')
    ax2.set_xscale('log')
    ax.legend(loc='upper right',frameon=False,fontsize=legend_size-3)
    ax2.legend(loc='upper right',frameon=False,fontsize=legend_size-3)

    plt.xlim((min_lums[0],10**5.5))
    ax.set_ylim((0,15))
    ax2.set_ylim((0,15))
    ax.tick_params(axis='both', which='major', labelsize=tick_size)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.subplots_adjust(hspace = 0.1)
    plt.gcf().subplots_adjust(bottom=0.13)
    extra=''
    plt.savefig('LMCPlots/LMC_in_50kpcTidal'+extra+'.pdf') #'.png')
    plt.close()









# Plot the GK16 original, and with 30% stripped and 70% stripped.
def Steeper(verbose=False):  
    f_left = 1 # fraction of remaining stars after tidal stripping
    pg = -0.2  # for playing with plot gamma. -0.2 is default
    w=8; h=6
    f, (ax, ax2) = plt.subplots(nrows=2,ncols=1,sharex=True, figsize=(w,h*1.6))
    low=2; high=7.5
    min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)  
    
    # Plot Observed Galaxies
    n_greater = get_observed(min_lums)
    ax.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')
    ax2.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')
    
    
    ### Plot First Panel ###    
    # original
    pg = -.2  
    model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg)
    N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
    model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg)
    N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)
    N_visible, N_std, samples_LMC = getLMC(model,min_lums)
    all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
    N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
    #ax.fill_between(f_left*min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.2)
    ax.plot(min_lums, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='-',label='AM slope = 1.97')

    # steeper
    pg = -.5
    model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg)
    N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
    model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg)
    N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)
    N_visible, N_std, samples_LMC = getLMC(model,min_lums)
    all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
    N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
    #ax.fill_between(f_left*min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.2)
    ax.plot(min_lums, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='--',label='AM slope = 2.43')


    # steepest
    pg = -1.0
    model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg)
    N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
    model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg)
    N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)
    N_visible, N_std, samples_LMC = getLMC(model,min_lums)
    all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
    N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
    ax.fill_between(f_left*min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.2)
    ax.plot(min_lums, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='-.',label='AM slope = 3.31')
    ### End of first panel ###



    ##### Plot the 2nd panel
    # Then plot MW background
    z = 8
    pg = -1.0
    model = am.GK16_grow(hostSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)
    N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
    ax2.plot(f_left*min_lums, N_visible_MW, label='MW background', color='grey',lw=linewidth,ls='-')
    
    # Next plot SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)
    N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)
    
    # Plot LMC Alone
    N_visible, N_std, samples_LMC = getLMC(model,min_lums)
    ax2.plot(f_left*min_lums, N_visible, color='grey',lw=linewidth,ls='-.',label='LMC alone')

    # combine all samples. keep the length N_iter the same.
    all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
    N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
    ax2.fill_between(f_left*min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.2)
    # Want to plot SMC after LMC for ordering
    ax2.plot(f_left*min_lums, N_visible_SMC, label='SMC alone', color='grey',lw=linewidth,ls=':')
    ax2.plot(f_left*min_lums, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=linewidth+2,ls='-',label='Adjusted\nGK16 + Reion')   #model.label)

    font0 = FontProperties()
    font = font0.copy()
    #font.set_weight('bold')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.2)
    ax2.text(.23,.8,'AM slope = 3.31 \n $z_{\mathrm{reion}} = 8.3$',transform=ax2.transAxes,fontsize=legend_size-2, fontproperties=font, bbox=props)
    #ax2.text(.23,.89,'AM slope = 3.31',transform=ax2.transAxes,fontsize=legend_size-2, fontproperties=font, bbox=props)
    ##### End plot the 2nd panel


    ax.set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=label_font)
    ax2.set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=label_font)
    ax2.set_xlabel('$\mathrm{M_*} \, [\mathrm{M_\odot}]$', fontsize=label_font)
    ax.set_xscale('log')
    ax2.set_xscale('log')
    ax.legend(loc='upper right',frameon=False,fontsize=legend_size-3)
    ax2.legend(loc='upper right',frameon=False,fontsize=legend_size-3)

    plt.xlim((min_lums[0],10**5.5))
    ax.set_ylim((0,15))
    ax2.set_ylim((0,15))
    ax.tick_params(axis='both', which='major', labelsize=tick_size)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.subplots_adjust(hspace = 0.1)
    plt.gcf().subplots_adjust(bottom=0.13)
    extra=''
    plt.savefig('LMCPlots/LMC_in_50kpcSteeper'+extra+'.pdf')  #'.png')
    plt.close()






def Bent(verbose=False):  
    f_left = 1 # fraction of remaining stars after tidal stripping
    pg = -0.2  # for playing with plot gamma. -0.2 is default
    w=8; h=6
    f, (ax, ax2) = plt.subplots(nrows=2,ncols=1,sharex=False, figsize=(w,h*1.8))
    low=2; high=7.5
    min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)  
    
    # Plot Observed Galaxies
    n_greater = get_observed(min_lums)
    ax2.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')
    
    # variables
    sigdex = 0.4
    pg = -0.7
    # end variables

    ### Plot the 1st panel #####
    a = 1.0
    M = 10**np.arange(7.5,np.log10(5e10),.05)
    AMmodels= [am.Sawala(), am.GK16_grow(bent=True, plotgamma = pg, bentsigdex = sigdex)]
    AMlabels = ['Sawala AM model', 'Tuned Bent Model']
    colors = ['cyan','orchid']
    for model, label, color in zip(AMmodels, AMlabels,colors):
        mstar = model.getStellarMass(M,a)
        ax.plot(M, mstar, label=label,linewidth=5, color=color)
        if isinstance(model, am.GK16_grow):  #or isinstance(model, Sawala):
            mstar_up = model.getStellarMass_up1sig(M,a)
            mstar_down = model.getStellarMass_down1sig(M,a)
            ax.fill_between(M,mstar_up,mstar_down, facecolor=color, alpha=0.2)

    ax.legend(loc='upper left',frameon=False,fontsize=legend_size)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('$M_{\mathrm{peak}} \mathrm{[M_\odot]}$',fontsize=label_font)
    ax.set_ylabel('$M_* \mathrm{[M_\odot]}$',fontsize=label_font)
    ax.set_ylim((10**2,10**8))
    ax.set_xlim((10**7.5,4e10))
    ax.tick_params(axis='both', which='major', labelsize=tick_size)
    #plt.gca().tight_layout()	
    #plt.gcf().subplots_adjust(bottom=0.14)
    #plt.gcf().subplots_adjust(left=0.14)
    #plt.savefig('LMCPlots/BentByBaryons2') # compare to figure 4 in bent by baryons
    #plt.close()
    ############

    

    ##### Plot the 2nd panel
    # Then plot MW background
    model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg, bent=True,bentsigdex = sigdex)
    N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
    #ax2.plot(f_left*min_lums, N_visible_MW, label='MW background', color='grey',lw=linewidth,ls='-')
    
    # Next plot SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg, bent=True,bentsigdex = sigdex)
    N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)
    
    # Plot LMC Alone
    N_visible, N_std, samples_LMC = getLMC(model,min_lums)
    #ax2.plot(f_left*min_lums, N_visible, color='grey',lw=linewidth,ls='-.',label='LMC alone')

    # combine all samples. keep the length N_iter the same.
    all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
    N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
    ax2.fill_between(f_left*min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor='orchid', alpha=0.2)
    # Want to plot SMC after LMC for ordering
    #ax2.plot(f_left*min_lums, N_visible_SMC, label='SMC alone', color='grey',lw=linewidth,ls=':')
    ax2.plot(f_left*min_lums, N_visible + N_visible_SMC + N_visible_MW, color='orchid',lw=linewidth+2,ls='-',label='Tuned Bent Model\n+ Reion')

    ax2.set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=label_font)
    ax2.set_xlabel('$\mathrm{M_*} \, [\mathrm{M_\odot}]$', fontsize=label_font)
    ax2.set_xscale('log')
    ax2.legend(loc='upper right',frameon=False,fontsize=legend_size)
    ax2.set_xlim((min_lums[0],10**5.5))
    ax2.set_ylim((0,15))
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    ##### End plot the 2nd panel


    plt.subplots_adjust(hspace = 0.28)
    plt.gcf().subplots_adjust(bottom=0.1)
    extra=''
    plt.savefig('LMCPlots/LMC_in_50kpcBent'+extra+'.pdf')  #'.png')
    plt.close()


#Steeper()
#Bent()





def TalkPlot1(verbose=False): 
    f_left = 1 # fraction of remaining stars after tidal stripping
    pg = -0.2  # for playing with plot gamma. -0.2 is default
    f, ax = plt.subplots()
    low=2; high=7.5
    min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)  
    lw = 5

    # Plot Observed Galaxies 
    galx = dm.load_nearby_gx()
    near_LMC = galx['dist_LMC'] < 50   # I need real tidal radius of LMC/SCM. i.e. what is sub-substructure according to rockstar?
    near_LMC[(galx['Name']=='LMC')|(galx['Name']=='SMC')] = False     # take out the LMC and SMC from the list
    smass = galx['mstar'][near_LMC]
    n_greater = np.array([np.sum(smass > min_lums[i]) for i in range(len(min_lums))])
    ax.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=lw,linestyle='--')
    
    # Then plot MW background
    model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg)
    N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
    ax.plot(min_lums, N_visible_MW, color='grey',lw=lw,ls='-')

    # Next plot SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg)
    N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)

    # Plot LMC Alone
    lmc = True  # want to use the lmc SHMF
    for model in [am.GK16_grow(lmcSHMF=lmc,reionization=True, plotgamma=pg)]: 
        print model.label, 'on this model'
        N_visible, N_std, samples_LMC = getLMC(model,min_lums)

        if isinstance(model, am.GK16_grow):
            ax.plot(min_lums, N_visible, color='grey',lw=lw,ls='-.')
            # combine all samples. keep the length N_iter the same.
            all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
            N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
            ax.fill_between(min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.2)
            # Want to plot SMC after LMC for ordering
            ax.plot(min_lums, N_visible_SMC, color='grey',lw=lw,ls=':')
            ax.plot(min_lums, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=lw,ls='-',label=model.label)   #'LMC+SMC+MW') 

    ax.set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=25)
    ax.set_xlabel('$\mathrm{M_*} \, [\mathrm{M_\odot}]$', fontsize=25)
    plt.xscale('log')
    ax.legend(loc='upper right',frameon=False,fontsize=legend_size)

    plt.xlim((min_lums[0],10**5.5))
    plt.ylim((0,15))
    ax.tick_params(axis='both', which='major', labelsize=20)
    #plt.subplots_adjust(hspace = 0.1)
    plt.gcf().subplots_adjust(bottom=0.14,left=0.14,right=0.95,top=0.95)
    extra='1'
    plt.savefig('LMCPlots/LMC_in_50kpcTalk'+extra)  #'.png')
    plt.close()




def TalkPlot2(verbose=False): 
    f_left = 1 # fraction of remaining stars after tidal stripping
    pg = -0.2  # for playing with plot gamma. -0.2 is default
    f, ax = plt.subplots()
    low=2; high=7.5
    min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)  
    lw = 5

    # Plot Observed Galaxies 
    galx = dm.load_nearby_gx()
    near_LMC = galx['dist_LMC'] < 50   # I need real tidal radius of LMC/SCM. i.e. what is sub-substructure according to rockstar?
    near_LMC[(galx['Name']=='LMC')|(galx['Name']=='SMC')] = False     # take out the LMC and SMC from the list
    smass = galx['mstar'][near_LMC]
    n_greater = np.array([np.sum(smass > min_lums[i]) for i in range(len(min_lums))])
    ax.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=lw,linestyle='--')
    

    # Then plot MW background
    model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg)
    N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
    ax.plot(min_lums, N_visible_MW, color='grey',lw=lw,ls='-')

    model = am.Brook(hostSHMF=True,reionization=True)
    Brook_N_MW, Brook_std_MW, Brook_samples_MW = getMW(model,min_lums)
    model = am.Behroozi(hostSHMF=True,reionization=True)
    Beh_N_MW, Beh_std_MW, Beh_samples_MW = getMW(model,min_lums)


    # Next plot SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg)
    N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)

    model = am.Brook(lmcSHMF=True,reionization=True)
    Brook_N_SMC, Brook_std_SMC, Brook_samples_SMC = getSMC(model,min_lums)
    model = am.Behroozi(lmcSHMF=True,reionization=True)
    Beh_N_SMC, Beh_std_SMC, Beh_samples_SMC = getSMC(model,min_lums)


    # Plot LMC Alone
    lmc = True  # want to use the lmc SHMF
    for model in [am.GK16_grow(lmcSHMF=lmc,reionization=True, plotgamma=pg), am.Behroozi(lmcSHMF=lmc,reionization=True), am.Brook(lmcSHMF=lmc,reionization=True)]:   #  am.Brook(lmcSHMF=lmc,reionization=False),
        print model.label, 'on this model'
        N_visible, N_std, samples_LMC = getLMC(model,min_lums)

        if isinstance(model, am.GK16_grow):
            ax.plot(min_lums, N_visible, color='grey',lw=lw,ls='-.')
            # combine all samples. keep the length N_iter the same.
            all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
            N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
            ax.fill_between(min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor=model.color, alpha=0.2)
            # Want to plot SMC after LMC for ordering
            ax.plot(min_lums, N_visible_SMC, color='grey',lw=lw,ls=':')
            ax.plot(min_lums, N_visible + N_visible_SMC + N_visible_MW, color=model.color,lw=lw,ls='-',label=model.label)   #'LMC+SMC+MW') 

        if isinstance(model, am.Brook):
            if model.reionization:
                ax.plot(min_lums, N_visible + Brook_N_SMC + Brook_N_MW, color=model.color,lw=lw,ls='-',label=model.label)
            else:
                ax.plot(min_lums, N_visible + BrookNR_N_SMC + BrookNR_N_MW, color=model.color,lw=lw,ls='--',label=model.label)
        if isinstance(model, am.Behroozi):
            ax.plot(min_lums, N_visible + Beh_N_SMC + Beh_N_MW, color=model.color,lw=lw,ls='-',label=model.label)


    ax.set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=25)
    ax.set_xlabel('$\mathrm{M_*} \, [\mathrm{M_\odot}]$', fontsize=25)
    plt.xscale('log')
    ax.legend(loc='upper right',frameon=False,fontsize=legend_size)

    plt.xlim((min_lums[0],10**5.5))
    plt.ylim((0,15))
    ax.tick_params(axis='both', which='major', labelsize=20)
    plt.subplots_adjust(hspace = 0.1)
    plt.gcf().subplots_adjust(bottom=0.14,left=0.14,right=0.95,top=0.95)
    extra='2'
    plt.savefig('LMCPlots/LMC_in_50kpcTalk'+extra)  #'.png')
    plt.close()





def BentTalk(verbose=False):  
    f_left = 1 # fraction of remaining stars after tidal stripping
    pg = -0.2  # for playing with plot gamma. -0.2 is default
    w=8; h=6
    f, (ax, ax2) = plt.subplots(nrows=2,ncols=1,sharex=False, figsize=(w,h*1.8))
    low=2; high=7.5
    min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)  
    
    # Plot Observed Galaxies
    n_greater = get_observed(min_lums)
    ax2.plot(dm.step_plotX(min_lums), dm.step_plotY(n_greater), label='Galaxy Candidates', color='black',lw=linewidth+1,linestyle='--')
    
    # variables
    sigdex = 0.4
    pg = -0.7
    # end variables

    ### Plot the 1st panel #####
    a = 1.0
    M = 10**np.arange(7.5,np.log10(10**11.01),.05)
    AMmodels= [am.Sawala(), am.GK16_grow(), am.Moster(), am.Brook(),am.GK16_grow(bent=True, plotgamma = pg, bentsigdex = sigdex, color='orchid',label='Tuned Bent Model')]  # am.Behroozi(),  #  am.GarrisonKimmel(),
    for model in AMmodels:
        mstar = model.getStellarMass(M,a)
        ax.plot(M, mstar, label=model.label,linewidth=3.5, color=model.color)
        if isinstance(model, am.GK16_grow) and model.color == 'orchid':  #or isinstance(model, Sawala):
            mstar_up = model.getStellarMass_up1sig(M,a)
            mstar_down = model.getStellarMass_down1sig(M,a)
            ax.fill_between(M,mstar_up,mstar_down, facecolor=model.color, alpha=0.2)

    ax.legend(loc='upper left',frameon=False,fontsize=legend_size)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('$M_{\mathrm{peak}} \mathrm{[M_\odot]}$',fontsize=25 )
    ax.set_ylabel('$M_* \mathrm{[M_\odot]}$',fontsize=25)
    ax.set_ylim((10**2,10**9))
    ax.set_xlim((10**7.5,10**11))
    ax.tick_params(axis='both', which='major', labelsize=20 )
    #plt.gca().tight_layout()	
    #plt.gcf().subplots_adjust(bottom=0.14)
    #plt.gcf().subplots_adjust(left=0.14)
    #plt.savefig('LMCPlots/BentByBaryons2') # compare to figure 4 in bent by baryons
    #plt.close()
    ############

    

    ##### Plot the 2nd panel
    # Then plot MW background
    model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg, bent=True,bentsigdex = sigdex)
    N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
    #ax2.plot(f_left*min_lums, N_visible_MW, label='MW background', color='grey',lw=linewidth,ls='-')
    
    # Next plot SMC Background
    model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg, bent=True,bentsigdex = sigdex)
    N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)
    
    # Plot LMC Alone
    N_visible, N_std, samples_LMC = getLMC(model,min_lums)
    #ax2.plot(f_left*min_lums, N_visible, color='grey',lw=linewidth,ls='-.',label='LMC alone')

    # combine all samples. keep the length N_iter the same.
    all_samples = np.array([np.append( np.append(samples_MW[i], samples_SMC[i]), samples_LMC[i]) for i in range(len(samples_MW))])
    N_visible_full, N_std_full, p20_full, p80_full = ambase.Nvisible_From_Ngreater(min_lums, all_samples)
    ax2.fill_between(f_left*min_lums, N_visible_full+N_std_full, N_visible_full-N_std_full, facecolor='orchid', alpha=0.2)
    # Want to plot SMC after LMC for ordering
    #ax2.plot(f_left*min_lums, N_visible_SMC, label='SMC alone', color='grey',lw=linewidth,ls=':')
    ax2.plot(f_left*min_lums, N_visible + N_visible_SMC + N_visible_MW, color='orchid',lw=linewidth+2,ls='-',label='Tuned Bent Model\n+ Reion')

    ax2.set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=25)
    ax2.set_xlabel('$\mathrm{M_*} \, [\mathrm{M_\odot}]$', fontsize=25)
    ax2.set_xscale('log')
    ax2.legend(loc='upper right',frameon=False,fontsize=20)
    ax2.set_xlim((min_lums[0],10**5.5))
    ax2.set_ylim((0,15))
    ax2.tick_params(axis='both', which='major', labelsize=20)
    ##### End plot the 2nd panel


    plt.subplots_adjust(hspace = 0.28)
    plt.gcf().subplots_adjust(bottom=0.1,left=0.14,right=0.96,top=0.96)
    extra='TALK'
    plt.savefig('LMCPlots/LMC_in_50kpcBent'+extra)  #'.png')
    plt.close()
