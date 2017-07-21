import numpy as np
import abundance_matching as am
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PlotParams import *
import RadialDependence as rd


#And  XVIII  - distance 1355,   stellar mass 0.63
#And XXVIII -  distance 661    stellar mass 0.21



#Sag DIG  for Sagittarius dIrr
DwarfNames = ['Leo T', 'And XXVIII', 'KKR 3', 'Tucana', 'Leo P', 'And XVIII', 'Phoenix', 'KKH 86','Antlia','KKR 25','Aquarius','DDO 113','Cetus','ESO 294-G010','Sagittarius dIrr','ESO 410-G005','KKH 98','Leo A','GR 8','Pegasus dIrr','UGC 9128','UGC 4879', 'KK 258',  'UGCA 86','DDO 99','UKS 2323-326','UGC 8508', 'KKs 3', 'NGC 4163','WLM','Sextans A','DDO 125','DDO 190','Sextans B','IC 3104','NGC 3109','NGC 6822','IC 1613','IC 4662','IC 5152']

StellarMasses = 10**6*np.array([0.14,0.21,  0.54,0.56, 0.56, 0.63, 0.77,0.82,1.3,1.4,1.6,2.1,2.6,2.7,3.5,3.5,4.5,6.0,6.4,6.6,7.8,8.3,14, 16,16,17,19,23,37,43,44,47,51,52,62,76,100,100,190,270])

Distances = np.array([ 417, 661, 2188, 887, 1620, 1355, 415, 2582, 1349, 1905, 1072, 2951, 755, 2032, 1067, 1923, 2523, 798, 2178, 920, 2291, 1361, 2230, 2965, 2594, 2208, 2582, 2120,  2858, 933, 1432, 2582, 2793, 1426, 2270, 1300, 459, 755, 2443, 1950 ])  # in kpc

N_iter = 20000  # make 20000 for the real thing. 5000 min acceptable

def get_ngreater(mstar,min_lum,model):
        dm_masses = model.mstar_to_mhalo(mstar,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        return model.ngreater(dm_masses,min_lum, percentile=True)

def get_P_at_least_one(mstar,min_lum,model):
        dm_masses = model.mstar_to_mhalo(mstar,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        return model.P_at_least_one(dm_masses,min_lum)

def get_rvir(StellarMasses,model):
        dm_masses = model.mstar_to_mhalo(StellarMasses,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        #dummy, rvir = dm.convert_Mhalo_d_d(dm_masses, model.delta_crit, 103.86)
        #return rvir
        return model.mass_to_rvir(dm_masses)  # this does the correct conversion

# just want the multiplicative fraction
def get_in_fov(radii, dists, camera_radius = 0.56):
        Z = 1.5
        fov = camera_radius * np.pi/180 # convert to radians
        fov_radius = fov * dists # convert to kpc at dist of target
        R = fov_radius/radii
        frac = rd.get_normed_los_better(R,Z) # rd.get_normed_los_protected(R,Z) 
        return frac
        

def add_sum_distr_ax(ax,min_lum,subset=True,fixedDM=False, re=False):
    if subset:
        StellarMassesSubset = StellarMasses[-5:]
    else:
        StellarMassesSubset = StellarMasses

    import scipy.misc

    for model in [am.GarrisonKimmel(reionization=re), am.GK16_grow(reionization=re),am.Moster(reionization=re),am.Brook(reionization=re)]:  #am.Sawala(),am.GarrisonKimmel16(reionization=re)
        if fixedDM:  # fix the dm mass to just moster
            tmpmodel = am.Moster()
            dm_masses = tmpmodel.mstar_to_mhalo(StellarMassesSubset,a=1)  #stellar_to_halo_mass(StellarMassesSubset,a=1.0)
            if subset:
                print np.log10(dm_masses), 'moster masses'
                print model.mstar_to_mhalo(StellarMassesSubset,a=1.0), model.label, 'masses'
        else:
            dm_masses = model.mstar_to_mhalo(StellarMassesSubset,a=1) #model.mstar_to_mhalo(StellarMassesSubset,a=1)  #stellar_to_halo_mass(StellarMassesSubset,a=1.0)
        print model.label
        if model.isPoisson:
            mean = model.get_field_total(dm_masses,min_lum)
            nocc = np.arange(0,mean*3+1)
            prob = (mean**nocc * np.e**-mean /scipy.misc.factorial(nocc))
            tmp=np.repeat(nocc,2)
        else:
        # for monte carlo methods, need to do 10,000 instances, looping over each halo and collecting all 10,000 numbers
        # then loop over all halos, and keep adding to the running total list.
        # then make histogram of the final 10,000 long distribution, normalized to 1
            samples = model.get_field_distr(dm_masses,min_lum,N=N_iter)
            distr,nocc = np.histogram(samples,bins=np.arange(min(samples),max(samples)+2))
            prob = distr/float(len(samples))
            tmp=np.repeat(nocc[0:-1],2)  # last value of bins is not inclusive
        ax.plot(np.append(tmp[1:],tmp[-1]+1), np.repeat(prob,2),linewidth=5,color=model.color,label=model.label)

        #if isinstance(model, am.Moster):
        #    samples = model.get_field_distr(dm_massesFull,min_lum,N=N_iter)
        #    distr,nocc = np.histogram(samples,bins=np.arange(min(samples),max(samples)+2))
        #    prob = distr/float(len(samples))
        #    tmp=np.repeat(nocc[0:-1],2)  # last value of bins is not inclusive
        #    dall, = ax.plot(np.append(tmp[1:],tmp[-1]+1), np.repeat(prob,2),linewidth=linewidth,color='grey',linestyle='--')
    return



def plot_sum_distr_1panel(subset=True,fixedDM=False,re=True):
    #w,h = plt.figaspect(1.5)
    f=plt.figure()
    h=f.get_figheight()
    w=f.get_figwidth()
    plt.close()
    
    minlum =4

    fig = plt.figure(figsize=(w*1.2,h))
    ax1 =fig.add_subplot(1,1,1)
    add_sum_distr_ax(ax1,min_lum=10**minlum,subset=subset,fixedDM=fixedDM,re=re)
    
    ax1.set_ylabel('Probability',fontsize=27)
    ax1.set_xlabel('N satellites total',fontsize=27)

    ax1.tick_params(axis='both', which='major', labelsize=22)
    plt.gcf().subplots_adjust(bottom=0.14, left=0.14,right=0.96,top=0.95)
    ax1.text(.27,.85,'$M_* > 10^4 \, \mathrm{M_\odot}$',transform=ax1.transAxes,fontsize=legend_size+7)

    if minlum==5:
        ax1.set_xlim((0,22))    #((0,45))
    if minlum==3:
        ax1.set_xlim((0,45))
    if minlum==4:
        ax1.set_xlim((0,30))
        

    ax1.legend(loc='upper right',frameon=False,fontsize=legend_size+3)
    extra='TALK'+str(minlum)
    if re:
        extra+='reionization'
    if subset:
        plt.savefig('largest5_field_sats_1panel'+extra)
    else:
        plt.savefig('total_field_sats_1panel'+extra)




# Talk Figure
plot_sum_distr_1panel()
#plot_sum_distr_1panel(subset=False,re=True)

