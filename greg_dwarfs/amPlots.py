import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PlotParams import *
from abundance_matching import *
#from AMbase import *



def plotMstar_v_Mhalo_talk():
    a = 1.0
    M = 10**np.arange(7.0,12.5,.1)
    #AMmodels=[GK16_grow()]
    AMmodels=[Behroozi(), GarrisonKimmel(), GK16_grow(), Moster(), Brook()]
    for model in AMmodels:
        if isinstance(model, Brook):
            # 10^7 to 10^8 Mstar
            Mlow = np.log10(model.stellar_to_halo_mass(10**7))
            Mhigh = np.log10(model.stellar_to_halo_mass(10**8))
            print Mlow, Mhigh, 'mlow and high'
            Mtmp = 10**np.arange(Mlow,Mhigh,.1)
            mstar = model.getStellarMass(Mtmp,a)
            plt.plot(Mtmp, mstar, label=model.label,lw=linewidth, color=model.color, ls='-')
            
            Mtmp = 10**np.arange(7.0,Mlow,.1)
            mstar = model.getStellarMass(Mtmp,a)
            plt.plot(Mtmp, mstar, lw=linewidth, color=model.color, ls='--')
            
            Mtmp = 10**np.arange(Mhigh,12.5,.1)
            mstar = model.getStellarMass(Mtmp,a)
            plt.plot(Mtmp, mstar, lw=linewidth, color=model.color, ls='--')

        else:

            if isinstance(model, GK16_grow):  #or isinstance(model, Sawala):
                mstar_up = model.getStellarMass_up1sig(M,a)
                mstar_down = model.getStellarMass_down1sig(M,a)
                plt.fill_between(M,mstar_up,mstar_down, facecolor=model.color, alpha=0.25)
                #mstar_rand = model.getStellarMass_random(M,a)
                #plt.scatter(M,mstar_rand,color=model.color)
                Mlow = np.log10(model.stellar_to_halo_mass(4.5*10**5))
                print Mlow, 'Mlow'

            if isinstance(model, GarrisonKimmel):
                Mlow = np.log10(model.stellar_to_halo_mass(10**8))
            if isinstance(model, Behroozi):
                Mlow = np.log10(8*10**10)
            if isinstance(model, Moster):
                Mlow = np.log10(model.stellar_to_halo_mass(10**7.4))

            # solid line
            Mtmp = 10**np.arange(Mlow,12.5,.1)
            mstar = model.getStellarMass(Mtmp,a)
            plt.plot(Mtmp, mstar, label=model.label,lw=linewidth, color=model.color)
        
            #extrapolation
            Mtmp = 10**np.arange(7.0,Mlow,.1)
            mstar = model.getStellarMass(Mtmp,a)
            plt.plot(Mtmp, mstar,lw=linewidth, color=model.color, ls='--')

        #else:
        #    mstar = model.getStellarMass(M,a)
        #    plt.plot(M, mstar, label=model.label,lw=linewidth, color=model.color)



    plt.legend(loc='lower right',fontsize=legend_size,frameon=False)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('$\mathrm{M_{halo}} \, \mathrm{[M_\odot]}$',fontsize=25)
    plt.ylabel('$\mathrm{M_*} \, \mathrm{[M_\odot]}$',fontsize=25)
    plt.ylim((10**2,10**11))
    plt.xlim((10**7, 2*10**12))
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.gcf().subplots_adjust(bottom=0.15,left=0.15,right=0.96,top=0.95)
    plt.savefig('smhm_relation_TALK') # compare to figure 5 in paper
    plt.close()


