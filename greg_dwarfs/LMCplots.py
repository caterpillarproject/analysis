import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import abundance_matching as am
from PlotParams import *
import DwarfMethods as dm
import RadialDependence as rd
import AMbase as ambase

vmax_ach = 9.48535156   # multiply by 0.6 for GK16 with z=13.  For z=9, Moster use 0.8, GK16 use 0.85
vmax_filt = 23.54248047  # has little effect
z = 11  # choices of 14,13,12,11,9

LMC_stellar_mass = 2.6*10**9
SMC_stellar_mass = 7.1*10**8


#luminosity = 10**np.linspace(4,9,16)
#halo_masses = 10**np.linspace(9,11.5,20)
ML_ratio = 1 #mass to luminosity ratio from Brook 2014 was 1.6. 


def convert_Mvir(Mvir, model):
    if isinstance(model, am.Moster) or isinstance(model, am.Sawala):
        return dm.convert_Mhalo_z0(Mvir, 200)
    if isinstance(model, am.GarrisonKimmel16) or isinstance(model, am.GarrisonKimmel) or isinstance(model, am.Behroozi) or isinstance(model, am.GK16_grow):
        return Mvir
    if isinstance(model, am.Brook):
        return dm.convert_Mhalo_z0(Mvir, 350)


# LMC Sized host and SMC sized host plot.
def plotFigure1(fixdmhalo=False):    # ngreater_v_minlum_2panel
    w=8; h=6
    fig, (ax1, ax2) = plt.subplots(nrows=2,ncols=1,sharex=True, figsize=(w,h*1.6) )
    min_lums = 10**np.linspace(3,7,25)
    lmc = True  # want to use the lmc SHMF

    # FIRST DO TOP PANEL
    #for model in [am.GarrisonKimmel(lmcSHMF=lmc,reionization=True,  catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt ), am.GK16_grow(lmcSHMF=lmc,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt),  am.Moster(lmcSHMF=lmc,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt), am.Brook(lmcSHMF=lmc,reionization=True,catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt),     am.GK16_grow(lmcSHMF=lmc,reionization=True),  am.GK16_grow(lmcSHMF=lmc,reionization=True, catreion=True,z=9, vmax_ach=vmax_ach, vmax_filt=vmax_filt)]:
    for model in [am.Behroozi(lmcSHMF=lmc,reionization=True), am.GarrisonKimmel(lmcSHMF=lmc,reionization=True), am.GK16_grow(lmcSHMF=lmc,reionization=True), am.Moster(lmcSHMF=lmc,reionization=True), am.Brook(lmcSHMF=lmc,reionization=True),   am.GK16_grow(lmcSHMF=lmc,reionization=True, catreion=True,z=11, vmax_ach=vmax_ach, vmax_filt=vmax_filt), am.GK16_grow(lmcSHMF=lmc,reionization=True, catreion=True,z=9, vmax_ach=vmax_ach, vmax_filt=vmax_filt)]:
        print model.label, 'on this model'
        if fixdmhalo:
            tmpmodel=am.Moster()
            halo_mass = tmpmodel.mstar_to_mhalo(LMC_stellar_mass, a=1.0)
            halo_mass = convert_Mvir(halo_mass,model)
        else:
            halo_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1.0)

        N_visible, N_std = model.ngreater(halo_mass,min_lums)

        if isinstance(model, am.GK16_grow):
            if not model.catreion:  # model.catreion and model.z == 11: 
                ax1.fill_between(min_lums, N_visible+N_std, N_visible-N_std, facecolor=model.color, alpha=0.2)
                ax1.plot(min_lums, N_visible,label=model.label,color=model.color,lw=linewidth,ls=model.ls)
            elif  model.z == 11:
                ax1.plot(min_lums, N_visible,label='GK16 ($z_{\mathrm{reion}} = 11.3$)',color=model.color,lw=linewidth,ls='--')
            else:
                ax1.plot(min_lums, N_visible,label='GK16 ($z_{\mathrm{reion}} = 9.3$)',color=model.color,lw=linewidth,ls='-.') # used to be 9
        else:
            ax1.plot(min_lums, N_visible,label=model.label,color=model.color,lw=linewidth,ls=model.ls)

    # THEN DO THE SECOND PANEL
    #for model in [am.GarrisonKimmel(lmcSHMF=lmc,reionization=True), am.GK16_grow(lmcSHMF=lmc,reionization=True),am.Moster(lmcSHMF=lmc,reionization=True),am.Brook(lmcSHMF=lmc,reionization=True)]:
    #for model in [am.GarrisonKimmel(lmcSHMF=lmc,reionization=True,  catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt ), am.GK16_grow(lmcSHMF=lmc,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt),  am.Moster(lmcSHMF=lmc,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt), am.Brook(lmcSHMF=lmc,reionization=True,catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt),     am.GK16_grow(lmcSHMF=lmc,reionization=True),  am.GK16_grow(lmcSHMF=lmc,reionization=True, catreion=True,z=9, vmax_ach=vmax_ach, vmax_filt=vmax_filt)]:
    for model in [am.Behroozi(lmcSHMF=lmc,reionization=True), am.GarrisonKimmel(lmcSHMF=lmc,reionization=True), am.GK16_grow(lmcSHMF=lmc,reionization=True),  am.Moster(lmcSHMF=lmc,reionization=True), am.Brook(lmcSHMF=lmc, reionization=True),  am.GK16_grow(lmcSHMF=lmc,reionization=True, catreion=True,z=11, vmax_ach=vmax_ach, vmax_filt=vmax_filt),  am.GK16_grow(lmcSHMF=lmc,reionization=True, catreion=True,z=9, vmax_ach=vmax_ach, vmax_filt=vmax_filt)]:
        print model.label, 'on this model'
        if fixdmhalo:
            tmpmodel=am.Moster()
            halo_mass = tmpmodel.mstar_to_mhalo(SMC_stellar_mass, a=1.0) 
            halo_mass = convert_Mvir(halo_mass,model)
        else:
            halo_mass = model.mstar_to_mhalo(SMC_stellar_mass, a=1.0)

        N_visible, N_std = model.ngreater(halo_mass,min_lums)

        if isinstance(model, am.GK16_grow):
            if not model.catreion:  #model.catreion and model.z == 11: 
                ax2.fill_between(min_lums, N_visible+N_std, N_visible-N_std, facecolor=model.color, alpha=0.2)
                ax2.plot(min_lums, N_visible,label=model.label,color=model.color,lw=linewidth,ls=model.ls)
            elif model.z == 11:   #not model.catreion:
                ax2.plot(min_lums, N_visible,label='GK16 ($z_{\mathrm{reion}} = 11.3$)',color=model.color,lw=linewidth,ls='--')
            else:
                ax2.plot(min_lums, N_visible,label='GK16 ($z_{\mathrm{reion}} = 9.3$)',color=model.color,lw=linewidth,ls='-.')
        else:
            ax2.plot(min_lums, N_visible,label=model.label,color=model.color,lw=linewidth,ls=model.ls)                


    ax1.text(.2,.90,'LMC Sized Host',transform=ax1.transAxes,fontsize=legend_size-2)
    ax2.text(.2,.90,'SMC Sized Host',transform=ax2.transAxes,fontsize=legend_size-2)
    ax1.set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=label_font)
    ax2.set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=label_font)
    plt.xscale('log')
    
    ax1.set_ylim(ymin=0)
    ax2.set_ylim(ymin=0)
    ax1.set_ylim(ymax=25)   # can be adjusted
    ax2.set_ylim(ymax=14)


    ax1.legend(loc='upper right',frameon=False,fontsize=legend_size-2)
    ax2.legend(loc='upper right',frameon=False,fontsize=legend_size-2)

    ax2.set_xlabel('$\mathrm{M_*} \, [\mathrm{M_\odot}]$', fontsize=label_font)
    plt.xlim((min_lums[0],min_lums[-1]))

    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.subplots_adjust(hspace = 0.1)
    plt.gcf().subplots_adjust(bottom=0.15)

    extra=''
    if fixdmhalo:
        extra+='fixed_dm_mass'
    plt.savefig('LMCPlots/LMC_Ngreater_vs_minlum'+extra+'.pdf')  #.png')
    plt.close()




# Plotting Ngreater 3 panel.
def plotFigure2(min_lum=[10**3, 10**4, 10**5]):
    w=8; h=6
    fig, (ax1, ax2,ax3) = plt.subplots(nrows=3,ncols=1,sharex=True, figsize=(w,h*1.95) )
    stellar_masses = 10**np.linspace(8.5,10,14)
    #halo_masses = 10**np.linspace(9,11.5,20)
    for model in [am.GarrisonKimmel(reionization=True, lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt), am.GK16_grow(reionization=True, lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt),am.Moster(reionization=True,lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt), am.Brook(reionization=True,lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt)]:
    #for model in [am.GarrisonKimmel(reionization=True)]:
        print model.label, 'on this model'
        halo_masses = model.mstar_to_mhalo(stellar_masses)    #stellar_to_halo_mass(stellar_masses)
        N_visible,_ = model.ngreater(halo_masses,min_lum[0])
        ax1.plot(stellar_masses, N_visible,label=model.label,color=model.color,linestyle=model.ls,lw=linewidth)
        
        N_visible,_ = model.ngreater(halo_masses,min_lum[1])
        ax2.plot(stellar_masses, N_visible,label=model.label,color=model.color,linestyle=model.ls,lw=linewidth)
        
        N_visible,_ = model.ngreater(halo_masses,min_lum[2])
        ax3.plot(stellar_masses, N_visible,label=model.label,color=model.color,linestyle=model.ls,lw=linewidth)

    ax1.set_ylabel('$\mathrm{ N_{sats}} > 10^%d \, \mathrm{M_\odot}$' %np.log10(min_lum[0]),fontsize=label_font)
    ax2.set_ylabel('$\mathrm{N_{sats}} > 10^%d \, \mathrm{M_\odot}$' %np.log10(min_lum[1]),fontsize=label_font)
    ax3.set_ylabel('$\mathrm{N_{sats}} > 10^%d \, \mathrm{M_\odot}$' %np.log10(min_lum[2]),fontsize=label_font)
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax3.set_xscale('log')

    ax1.text(.04,.37,'$M_* > 10^3 \, \mathrm{M_\odot}$',transform=ax1.transAxes,fontsize=legend_size)
    ax2.text(.04,.37,'$M_* > 10^4 \, \mathrm{M_\odot}$',transform=ax2.transAxes,fontsize=legend_size)
    ax3.text(.04,.37,'$M_* > 10^5 \, \mathrm{M_\odot}$',transform=ax3.transAxes,fontsize=legend_size)

    ax1.set_ylim((0,35)); ax2.set_ylim((0,23)); ax3.set_ylim((0,9))

    yvals = ((0,5,10,15,20,25,30,35))
    ax1.set_yticks(yvals)
    ax1.set_yticklabels([str(y) for y in yvals])
    yvals = ((0,4,8,12,16,20))
    ax2.set_yticks(yvals)
    ax2.set_yticklabels([str(y) for y in yvals])
    yvals = ((0,1,2,3,4,5,6,7,8))
    ax3.set_yticks(yvals)
    ax3.set_yticklabels([str(y) for y in yvals])

 
    ax1.legend(loc='upper left',frameon=False,fontsize=legend_size)
    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    ax3.tick_params(axis='both', which='major', labelsize=tick_size)
    
    ax1.vlines([LMC_stellar_mass, SMC_stellar_mass], ax1.get_ylim()[0], ax1.get_ylim()[1], linestyles='--',color='black')
    ax2.vlines([LMC_stellar_mass, SMC_stellar_mass], ax2.get_ylim()[0], ax2.get_ylim()[1], linestyles='--',color='black')
    ax3.vlines([LMC_stellar_mass, SMC_stellar_mass], ax3.get_ylim()[0], ax3.get_ylim()[1], linestyles='--',color='black')

    plt.subplots_adjust(hspace = 0.0)
    plt.gcf().subplots_adjust(bottom=0.15)

    ax3.set_xlabel('$M_* \ \mathrm{of \ Host} \ \mathrm{[M_\odot]}$',fontsize=label_font)
    
    ax1.set_xlim((stellar_masses[0],stellar_masses[-1]))
    ax2.set_xlim((stellar_masses[0],stellar_masses[-1]))
    ax3.set_xlim((stellar_masses[0],stellar_masses[-1]))

    plt.savefig('LMCPlots/LMC_Ngreater3panel.png')
    plt.close()







# Plotting Ngreater 2 panel.
def plotNgreater(min_lum=10**5):
    w=8; h=6
    fig, (ax1, ax2) = plt.subplots(nrows=2,ncols=1,sharex=True, figsize=(w,h*1.6) )
    stellar_masses = 10**np.linspace(8.5,10,14)
    #halo_masses = 10**np.linspace(9,11.5,20)
    for model in [am.GarrisonKimmel(reionization=True, lmcSHMF=True), am.GK16_grow(reionization=True, lmcSHMF=True),am.Moster(reionization=True,lmcSHMF=True), am.Brook(reionization=True,lmcSHMF=True)]:
        print model.label, 'on this model'
        halo_masses = model.mstar_to_mhalo(stellar_masses)    #stellar_to_halo_mass(stellar_masses)
        N_visible, N_std = model.ngreater(halo_masses,min_lum)
        ax1.plot(stellar_masses, N_visible,label=model.label,color=model.color,linestyle=model.ls,lw=linewidth)
        
        #N_visible,_ = model.ngreater(halo_masses,min_lum)
        # scale to 100 kpc
        rvir = model.mass_to_rvir(halo_masses)
        factor =  rd.get_normed_radial_abundance(100 / rvir)
        ax2.plot(stellar_masses, factor*N_visible,label=model.label,color=model.color,linestyle=model.ls,lw=linewidth)

        if isinstance(model, am.GK16_grow):
            ax1.fill_between(stellar_masses, N_visible+N_std, N_visible-N_std, facecolor=model.color, alpha=0.2)
            ax2.fill_between(stellar_masses, factor*(N_visible+N_std), factor*(N_visible-N_std), facecolor=model.color, alpha=0.20)
            

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.2)
    ax1.text(.70,.87,'Within $R_{\mathrm{vir}}$',transform=ax1.transAxes,fontsize=legend_size-1, bbox=props)
    ax2.text(.65,.87,'Within 100 kpc',transform=ax2.transAxes,fontsize=legend_size-1, bbox=props)
    ax1.text(.70,.75,'$M_* > 10^5 \, \mathrm{M_\odot}$',transform=ax1.transAxes,fontsize=legend_size)
    ax2.text(.65,.75,'$M_* > 10^5 \, \mathrm{M_\odot}$',transform=ax2.transAxes,fontsize=legend_size)


    ax1.set_ylabel('$\mathrm{ N_{sats}} > 10^%d \, \mathrm{M_\odot}$' %np.log10(min_lum),fontsize=label_font)
    ax2.set_ylabel('$\mathrm{N_{sats}} > 10^%d \, \mathrm{M_\odot}$' %np.log10(min_lum),fontsize=label_font)

    ax1.set_xscale('log')
    ax2.set_xscale('log')


    ax1.set_ylim((0,9));   ax2.set_ylim((0,9))

    yvals = ((0,1,2,3,4,5,6,7,8))
    ax1.set_yticks(yvals)
    ax1.set_yticklabels([str(y) for y in yvals])
    yvals = ((0,1,2,3,4,5,6,7,8))
    ax2.set_yticks(yvals)
    ax2.set_yticklabels([str(y) for y in yvals])

    ax1.legend(loc='upper left',frameon=False,fontsize=legend_size)
    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    
    ax1.vlines([LMC_stellar_mass, SMC_stellar_mass], ax1.get_ylim()[0], ax1.get_ylim()[1], linestyles='--',color='black')
    ax2.vlines([LMC_stellar_mass, SMC_stellar_mass], ax2.get_ylim()[0], ax2.get_ylim()[1], linestyles='--',color='black')

    plt.subplots_adjust(hspace = 0.0)
    plt.gcf().subplots_adjust(bottom=0.15)

    ax2.set_xlabel('$M_* \ \mathrm{of \ Host} \ \mathrm{[M_\odot]}$',fontsize=label_font)
    
    ax1.set_xlim((stellar_masses[0],stellar_masses[-1]))
    ax2.set_xlim((stellar_masses[0],stellar_masses[-1]))

    plt.savefig('LMCPlots/LMC_Ngreater2panel.pdf')  #.png
    plt.close()










def add_sum_distr_ax(ax,min_lum, re=True):
    StellarMassesSubset = LMC_stellar_mass
    import scipy.misc
    for model in [am.GarrisonKimmel(reionization=re, lmcSHMF=True), am.GK16_grow(reionization=re, lmcSHMF=True),am.Moster(reionization=re, lmcSHMF=True),am.Brook(reionization=re, lmcSHMF=True)]: 
        dm_masses = model.mstar_to_mhalo(StellarMassesSubset,a=1)
        print model.label
        samples = model.get_field_distr(dm_masses,min_lum)
        distr,nocc = np.histogram(samples,bins=np.arange(min(samples),max(samples)+2))
        prob = distr/float(len(samples))
        tmp=np.repeat(nocc[0:-1],2)  # last value of bins is not inclusive
        ax.plot(np.append(tmp[1:],tmp[-1]+1), np.repeat(prob,2),linewidth=linewidth,color=model.color,label=model.label)
    return


# Probability Distribution
# probability_distr.png
def plotFigure3():
    f=plt.figure()
    h=f.get_figheight()
    w=f.get_figwidth()
    plt.close()

    fig = plt.figure(figsize=(w,h*1.5))
    ax1 =fig.add_subplot(3,1,1)
    ax2 =fig.add_subplot(3,1,2)
    ax3 =fig.add_subplot(3,1,3)
    #f,(ax1,ax2,ax3) = plt.subplots(nrows=3)

    add_sum_distr_ax(ax1,min_lum=10**3)
    add_sum_distr_ax(ax2,min_lum=10**4)
    add_sum_distr_ax(ax3,min_lum=10**5)
    
    ax1.set_ylabel('Probability',fontsize=label_font)
    ax2.set_ylabel('Probability',fontsize=label_font)
    ax3.set_ylabel('Probability',fontsize=label_font)
    ax3.set_xlabel('N satellites total',fontsize=label_font)

    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    ax3.tick_params(axis='both', which='major', labelsize=tick_size)
    #plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().subplots_adjust(bottom=0.10)
    #plt.gcf().subplots_adjust(top=0.05)

    ax1.text(.35,.82,'$M_* > 10^3 \, \mathrm{M_\odot}$',transform=ax1.transAxes,fontsize=legend_size)
    ax2.text(.35,.8,'$M_* > 10^4 \, \mathrm{M_\odot}$',transform=ax2.transAxes,fontsize=legend_size)
    ax3.text(.35,.8,'$M_* > 10^5 \, \mathrm{M_\odot}$',transform=ax3.transAxes,fontsize=legend_size)
    
    #ax1.set_xlim((0,45))
    #ax2.set_xlim((0,35))
    #ax3.set_xlim((0,20))
    ax1.legend(loc='upper right',frameon=False,fontsize=legend_size-3)
    plt.savefig('LMCPlots/probability_distr.png')



def test_variables(re=True, min_lum=10**3):
    radius_small = np.arange(0,51)  # for the LMC
    model = am.GK16(reionization=re, lmcSHMF=True)
    dm_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
    lmc_rvir = model.mass_to_rvir(dm_mass)
    factor =  rd.get_normed_radial_abundance(radius_small / lmc_rvir)
    print min_lum, 'min_lum'
    ngreater, N_std = model.ngreater(dm_mass,min_lum=min_lum, percentile=False)
    print model.label
    print model.superPoissonian, 'super poisson?'
    print factor[50]*ngreater
    print factor[50]*(ngreater+N_std)

#test_variables(re=True)



##  3 CHANGE THIS PLOT! 
# Add step function of real data from the galaxies.
# Make the MW background uniform density - i.e. straight line from 0 to 1.75
# Make SMC radial distribution??? Or don't include this radial distribution because it is
# too complicated!!

def plotLMC_radial(re=True, min_lum=10**3):
    f, ax = plt.subplots()
    radius_small = np.arange(0,51)  # for the LMC
    for model in [am.GarrisonKimmel(reionization=re, lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt), am.GK16_grow(reionization=re, lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt),am.Moster(reionization=re, lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt),am.Brook(reionization=re, lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt)]:
        dm_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        lmc_rvir = model.mass_to_rvir(dm_mass)
        dm_mass_smc = model.mstar_to_mhalo(SMC_stellar_mass,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        smc_rvir = model.mass_to_rvir(dm_mass_smc)
        #print lmc_rvir, 'virial radius [kpc]'
        
        factor = rd.get_normed_radial_abundance(radius_small / lmc_rvir)
        ngreater, N_std = model.ngreater(dm_mass,min_lum=min_lum, percentile=False)
        factor_smc = rd.get_normed_radial_abundance(radius_small / smc_rvir)
        ngreater_smc, N_std_smc = model.ngreater(dm_mass_smc,min_lum=min_lum, percentile=False)
        ax.plot(radius_small, factor*ngreater, color=model.color,linestyle='-',lw=linewidth, label=model.label)
        ax.plot(radius_small, factor_smc*ngreater_smc, color=model.color,linestyle='--',lw=linewidth)
        if isinstance(model, am.Moster):
            print N_std, model.superPoissonian, 'standard deviation of mean within rvir, and is superpoisson on?'
            r_sigma = rd.get_normed_radial_std(radius_small / lmc_rvir)
            sig1 = factor*N_std
            sig2 = ngreater*r_sigma
            full_sigma = np.sqrt( sig1**2 + sig2**2)
            
            #ax.fill_between(radius_small, factor*(ngreater-N_std), factor*(ngreater+N_std), facecolor=model.color, alpha=0.25)
            ax.fill_between(radius_small, factor*ngreater+full_sigma, factor*ngreater-full_sigma, facecolor=model.color, alpha=0.25)
        print 'done with', model.label

    ax.set_ylabel('$\mathrm{ N_{sats}} > 10^3 \, \mathrm{M_\odot}$', fontsize=label_font)
    ax.set_xlabel('Radius [kpc]',fontsize=label_font)

    ax.text(.04,.86,'$M_* > 10^3 \, \mathrm{M_\odot}$ in LMC-like host',transform=ax.transAxes,fontsize=legend_size)
    ax.tick_params(axis='both', which='major', labelsize=tick_size)
    ax.vlines(35, ax.get_ylim()[0], ax.get_ylim()[1], linestyles='--',color='black')
    ax.legend(loc=(0.04,.45),frameon=False,fontsize=legend_size-3)
    plt.ylim((0,2.2))
    plt.xlim((0,50))
    plt.gcf().subplots_adjust(bottom=0.10)
    plt.savefig('LMCPlots/LMCradial_figure.png')


def plotLMC_los(re=True, min_lum=10**5):
    f, ax = plt.subplots()
    radius = np.arange(0,151)  # for LMC analogs
    #for model in [am.GarrisonKimmel(reionization=re, lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt), am.GK16_grow(reionization=re, lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt),am.Moster(reionization=re, lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt),am.Brook(reionization=re, lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt)]:
    for model in [am.GarrisonKimmel(reionization=re, lmcSHMF=True), am.GK16_grow(reionization=re, lmcSHMF=True),am.Moster(reionization=re, lmcSHMF=True),am.Brook(reionization=re, lmcSHMF=True)]:
        print model.label
        dm_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        #lmc_rvir = model.mass_to_rvir(dm_mass) # this was computing R200 and R350. I want the full Rvir
        dummy, lmc_rvir = dm.convert_Mhalo_d_d(dm_mass, model.delta_crit, 103.86)
        
        # for M_* > 10^5, distant LMC analogs
        ngreater, nstd = model.ngreater(dm_mass,min_lum=min_lum, percentile=False)
        #factor2 = rd.get_normed_los_protected(radius/lmc_rvir ,Z=1.5)    
        factor2 = rd.get_normed_los_better(radius/lmc_rvir, Z=1.5)

        ax.plot(radius, factor2*ngreater, color=model.color,linestyle=model.ls,lw=linewidth, label=model.label)
        if isinstance(model, am.GK16_grow):
            r_sig = np.append([0,5], np.arange(15,151,15))
            factors = rd.get_normed_los_better(r_sig/lmc_rvir, Z=1.5)
            samples_LMC = [model.generate_stellar_mass_samples(dm_mass, factor=factor) for factor in factors]
            N_std = [ambase.Nvisible_From_Ngreater(min_lum, samples)[1] for samples in samples_LMC]
            ax.fill_between(r_sig,  factor2[r_sig]*ngreater + N_std, factor2[r_sig]*ngreater - N_std, facecolor=model.color, alpha=0.2)
        print 'done with', model.label

    ax.set_ylabel('$\mathrm{ N_{sats}} > 10^5 \, \mathrm{M_\odot}$', fontsize=label_font)
    ax.set_xlabel('Field of view radius [kpc]',fontsize=label_font)

    ax.text(.04,.90,'$M_* > 10^5 \, \mathrm{M_\odot}$ in LMC-like host',transform=ax.transAxes,fontsize=legend_size)
    ax.tick_params(axis='both', which='major', labelsize=tick_size)
    #ax.vlines(35, ax.get_ylim()[0], ax.get_ylim()[1], linestyles='--',color='black')
    ax.legend(loc=(0.04,.57),frameon=False,fontsize=legend_size-3)
    plt.ylim((0,7))
    plt.xlim((0,150))
    plt.gcf().subplots_adjust(bottom=0.13)
    plt.savefig('LMCPlots/LMC_los_radial_figure.pdf')  #.png




# make the top panel just the LMC itself. Radius refers to spherical radius. and use > 10^3
# bottom panel for distance galaxies. Radius refers to the field of view radius with Z=1.5, and use >10^5
def plotFigure4(re=True, min_lum=10**5):
    f=plt.figure()
    h=f.get_figheight()
    w=f.get_figwidth()
    plt.close()

    fig = plt.figure(figsize=(w,h*1.5))
    ax1 =fig.add_subplot(2,1,1)
    ax2 =fig.add_subplot(2,1,2)
    radius = np.arange(0,150)  # for LMC analogs
    radius_small = np.arange(0,50)  # for the LMC

    for model in [am.GarrisonKimmel(reionization=re, lmcSHMF=True), am.GK16_grow(reionization=re, lmcSHMF=True),am.Moster(reionization=re, lmcSHMF=True),am.Brook(reionization=re, lmcSHMF=True)]:
        dm_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        lmc_rvir = model.mass_to_rvir(dm_mass)
        print lmc_rvir, 'virial radius [kpc]'

        # for M_* > 10^5, distant LMC analogs
        ngreater = model.ngreater(dm_mass,min_lum=10**5, percentile=False)[0]
        factor2 = rd.get_normed_los_protected(radius/lmc_rvir ,Z=1.5)
        ax2.plot(radius, factor2*ngreater, color=model.color,linestyle=model.ls,lw=linewidth )
        #print ngreater, lmc_rvir, 'ngreater and rvir'
        #factor = rd.getK_protected(radius_small / lmc_rvir)
        #ax1.plot(radius_small, factor*ngreater,label=model.label, color=model.color,linestyle=model.ls,lw=linewidth )

        # for M_* > 10^3, the actual LMC
        factor =  rd.get_normed_radial_abundance(radius_small / lmc_rvir) #rd.getK_protected(radius_small / lmc_rvir)
        ngreater = model.ngreater(dm_mass,min_lum=10**3, percentile=False)[0]
        ax1.plot(radius_small, factor*ngreater, color=model.color,linestyle='-',lw=linewidth, label=model.label )
        #factor2 = rd.get_normed_los_protected(radius/lmc_rvir ,Z=1.5)
        #ax2.plot(radius, factor2*ngreater, color=model.color,linestyle='--',lw=linewidth )
        print 'done with', model.label

    #ax1.set_ylim(0,3.5)
    ax1.set_ylabel('$\mathrm{ N_{sats}} > 10^3 \, \mathrm{M_\odot}$', fontsize=label_font)
    ax2.set_ylabel('$\mathrm{ N_{sats}} > 10^5 \, \mathrm{M_\odot}$', fontsize=label_font)
    ax2.set_xlabel('Radius [kpc]',fontsize=label_font)
    ax2.set_xlabel('Line of sight radius [kpc]',fontsize=label_font)

    ax1.text(.04,.86,'$M_* > 10^3 \, \mathrm{M_\odot}$ in LMC-like host',transform=ax1.transAxes,fontsize=legend_size)
    ax2.text(.04,.86,'$M_* > 10^5 \, \mathrm{M_\odot}$ in LMC-like host',transform=ax2.transAxes,fontsize=legend_size)

    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)

    ax1.vlines(35, ax1.get_ylim()[0], ax1.get_ylim()[1], linestyles='--',color='black')
    #ax2.vlines(35, ax2.get_ylim()[0], ax2.get_ylim()[1], linestyles='--',color='black')
    ax1.legend(loc=(0.04,.45),frameon=False,fontsize=legend_size-3)
    plt.gcf().subplots_adjust(bottom=0.10)
    plt.savefig('LMCPlots/radial_figure.png')




def plot_P_at_least_one_3panel(min_lum=[10**3, 10**4, 10**5]):
    w=8; h=6
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3,ncols=1,sharex=True, figsize=(w,h*1.95) )

    stellar_masses = 10**np.linspace(8,9.6,10)
    #halo_masses = 10**np.linspace(9,11.5,20)
    for model in [am.GarrisonKimmel(reionization=True),am.GK16_grow(reionization=True),am.Moster(reionization=True),am.Brook(reionization=True)]:
    #for model in [am.GarrisonKimmel(reionization=True)]:
        print model.label, 'on this model'
        halo_masses = model.mstar_to_mhalo(stellar_masses)   #stellar_to_halo_mass(stellar_masses)
        Pgt1 = model.P_at_least_one(halo_masses,min_lum[0])
        ax1.plot(stellar_masses, Pgt1,label=model.label,color=model.color,linestyle=model.ls,lw=linewidth)
        
        Pgt1 = model.P_at_least_one(halo_masses,min_lum[1])
        ax2.plot(stellar_masses, Pgt1,label=model.label,color=model.color,linestyle=model.ls,lw=linewidth)

        Pgt1 = model.P_at_least_one(halo_masses,min_lum[2])
        ax3.plot(stellar_masses, Pgt1,label=model.label,color=model.color,linestyle=model.ls,lw=linewidth)


    ax1.text(.08,.09,'$M_* > 10^3 \, \mathrm{M_\odot}$',transform=ax1.transAxes,fontsize=legend_size)
    ax2.text(.08,.09,'$M_* > 10^4 \, \mathrm{M_\odot}$',transform=ax2.transAxes,fontsize=legend_size)
    ax3.text(.08,.09,'$M_* > 10^5 \, \mathrm{M_\odot}$',transform=ax3.transAxes,fontsize=legend_size)

    ax1.set_ylabel('$\mathrm{P}(\geq 1 \ \mathrm{satellite})$',fontsize=label_font)
    ax2.set_ylabel('$\mathrm{P}(\geq 1 \ \mathrm{satellite})$',fontsize=label_font)
    ax3.set_ylabel('$\mathrm{P}(\geq 1 \ \mathrm{satellite})$',fontsize=label_font)
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax3.set_xscale('log')

    ax1.set_ylim((0.6,1));     ax2.set_ylim((0.4,1));     ax3.set_ylim((0.2,1))
    ax1.legend(loc='lower right',frameon=False,fontsize=legend_size-1)
    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)    
    ax3.tick_params(axis='both', which='major', labelsize=tick_size)    

    ax1.vlines([LMC_stellar_mass, SMC_stellar_mass], ax1.get_ylim()[0], ax1.get_ylim()[1], linestyles='--',color='black')
    ax2.vlines([LMC_stellar_mass, SMC_stellar_mass], ax2.get_ylim()[0], ax2.get_ylim()[1], linestyles='--',color='black')
    ax3.vlines([LMC_stellar_mass, SMC_stellar_mass], ax3.get_ylim()[0], ax3.get_ylim()[1], linestyles='--',color='black')

    plt.subplots_adjust(hspace = 0.0)
    plt.gcf().subplots_adjust(bottom=0.12)


    ax1.set_yticks((.6,.7,.8,.9,1))
    ax1.set_yticklabels(['0.6','0.7','0.8','0.9','1.0'])
    ax2.set_yticks((.4,.5,.6,.7,.8,.9))
    ax2.set_yticklabels(['0.4','0.5','0.6','0.7','0.8','0.9'])
    ax3.set_yticks((.2,.4,.6,.8))
    ax3.set_yticklabels(['0.2','0.4','0.6','0.8'])

    ax3.set_xlabel('$M_* \ \mathrm{of \ Host} \ [\mathrm{M_\odot}]$',fontsize=label_font)    
    ax1.set_xlim((stellar_masses[0],stellar_masses[-1]))
    ax2.set_xlim((stellar_masses[0],stellar_masses[-1]))
    ax3.set_xlim((stellar_masses[0],stellar_masses[-1]))
    
    plt.savefig('LMCPlots/LMC_ProbabilityOfOne3panel.png')
    plt.close()






def ngreater_v_minlum_1panel(fixdmhalo=False):
    w=8; h=6
    fig, (ax1) = plt.subplots(nrows=1,ncols=1,sharex=True, figsize=(w,h*.8) )
    min_lums = 10**np.linspace(3,7,25)

    # FIRST DO TOP PANEL
    mw = False  # want to use the field halo SHMF
    #for model in [am.GarrisonKimmel(hostSHMF=mw,reionization=True), am.GK16_grow(hostSHMF=mw,reionization=True)]:
    for model in [am.GarrisonKimmel(hostSHMF=mw,reionization=True), am.GK16_grow(hostSHMF=mw,reionization=True),am.Moster(hostSHMF=mw,reionization=True),am.Brook(hostSHMF=mw,reionization=True)]:
        print model.label, 'on this model'
        if fixdmhalo:
            tmpmodel=am.Moster()
            halo_mass = tmpmodel.mstar_to_mhalo(LMC_stellar_mass, a=1.0) 
            halo_mass = convert_Mvir(halo_mass,model)
        else:
            halo_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1.0)

        N_visible = np.array([model.ngreater(halo_mass,minlum)[0] for minlum in min_lums])
        ax1.plot(min_lums, N_visible,label=model.label,color=model.color,lw=linewidth,ls=model.ls)
        ax1.text(.2,.90,'LMC Sized Host',transform=ax1.transAxes,fontsize=legend_size-2)

   
    ax1.set_ylabel('$\mathrm{N_{sats}} > \mathrm{M_*}$',fontsize=label_font)
    plt.xscale('log')
    
    ax1.legend(loc='upper right',frameon=False,fontsize=legend_size-2)
    ax1.set_xlabel('$\mathrm{M_*} \, [\mathrm{M_\odot}]$', fontsize=label_font)
    plt.xlim((min_lums[0],min_lums[-1]))
    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.subplots_adjust(hspace = 0.1)
    plt.gcf().subplots_adjust(bottom=0.15)

    extra=''
    if fixdmhalo:
        extra+='fixed_dm_mass'
    plt.savefig('LMCPlots/LMC_Ngreater_vs_minlum_1panel.pdf')
    plt.close()





# FIGURE 1:  LMC_Ngreater_vs_minlum.png
plotFigure1(fixdmhalo=False)
#ngreater_v_minlum_1panel(fixdmhalo=False)

# FIGURE 2:
#plotFigure2()

# FIGURE 3: probability_distr.png
#plotFigure3()

# FIGURE 4: radial_figure.png
#plotLMC_radial(re=True, min_lum=10**3)
#plotFigure4()


# FIGURE 7: LMC_los_radial_figure.pdf
#plotLMC_los()

#plot_P_at_least_one_3panel()
