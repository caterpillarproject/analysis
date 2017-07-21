import numpy as np
import abundance_matching as am
from PlotParams import *
import DwarfMethods as dm
import RadialDependence as rd
from scipy.integrate import dblquad
from scipy.integrate import quad
import AMbase as ambase

LMC_stellar_mass = 2.6*10**9
SMC_stellar_mass = 7.1*10**8
vmax_ach = 9.48535156   # multiply by 0.6 for GK16 with z=13.  For z=9, Moster use 0.8, GK16 use 0.85
vmax_filt = 23.54248047  # has little effect


def get_ratio(model, galaxy):
    if galaxy == 'MW':
        if model.label == 'Moster + Reion':
            print 'in moster'
            return 0.0698049317676
        if model.label == 'Brook + Reion':
            print 'in brook'
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
            print 'in moster'
            return 0.462038583656
        if model.label == 'Brook + Reion':
            print 'in brook'
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


        


# Find ratios in LMC_obs.py background_sats() 
model = am.GK16_grow(lmcSHMF=True,reionization=True)
halo_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1.0)
lmc_rvir = model.mass_to_rvir(halo_mass)
factor =  rd.get_normed_radial_abundance(50 / lmc_rvir)

MW at 1.4e12 Msun host mass: 0.0720138 
MW at 0.7e12 Msun host mass: 0.110016653198

LMC at normal mass: 0.3660271951182939
LMC at 1/2 mass: 0.48956870866880375

SMC at normal mass: 0.38563
SMC at 1/2 mass: 0.519757096168


# ratios for original GK16
SMC: 0.399382697282
MW: 0.0766766197663

# Ratios for Brook
SMC: 0.455996567841
LMC: 
MW: 0.0648588855887


# Ratios for Moster
SMC: 0.462038583656
LMC: 
MW: 0.0698049317676

# Behroozi
SMC: 429648136194
MW: 0.0766766197664

# Ratios for GK16 pg = -0.5
SMC: 0.384700258627
MW:  0.0766766198533

# Ratios for GK16 pg = -1.0
SMC: 0.363079218041
MW: 0.0766766198627




low=2; high=7.5
pg = -0.2
min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)
model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg)
N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)

model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg)
N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)
N_visible_LMC, N_std_LMC, samples_LMC = getLMC(model,min_lums)
min_lums[20]


# normal masses > 10^4:
lmc = N_visible_LMC[20]
smc = N_visible_SMC[20]
mw = N_visible_MW[20]
lmc / (lmc+smc+mw)
smc / (lmc+smc+mw)
mw / (lmc+smc+mw)
# SMC = 1.645250 - 21% or 21,21
# LMC = 2.761750 - 33% or 34, 33
# MW = 3.5354999 - 46% or 45, 46
# total = 7.94249

# normal masses > 10^3:
lmc = N_visible_LMC[10]
smc = N_visible_SMC[10]
mw = N_visible_MW[10]
lmc / (lmc+smc+mw)
smc / (lmc+smc+mw)
mw / (lmc+smc+mw)
# SMC = 1.645250 - 21% or 22
# LMC = 2.761750 - 34%
# MW = 3.5354999 - 45% or 44


# normal masses > 10^2:
lmc = N_visible_LMC[0]
smc = N_visible_SMC[0]
mw = N_visible_MW[0]
lmc / (lmc+smc+mw)
smc / (lmc+smc+mw)
mw / (lmc+smc+mw)
# SMC = 1.645250 - 22% or 22
# LMC = 2.761750 - 35%
# MW = 3.5354999 - 44% or 44


# half masses
# SMC = 1.1074999
# LMC = 1.83375
# MW = 2.6425
# total = 5.58375
# SMC = 1.645250 - 20% 
# LMC = 2.761750 - 33%
# MW = 3.5354999 - 47%




# Brook model ratios
# the ratios need to be different for the Brook model
low=2; high=7.5
pg = -0.2
min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)
model = am.Brook(hostSHMF=True,reionization=True)
N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
model = am.Brook(lmcSHMF=True,reionization=True)
N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)
N_visible_LMC, N_std_LMC, samples_LMC = getLMC(model,min_lums)
min_lums[20]
# normal masses > 10^4:
lmc = N_visible_LMC[20]
smc = N_visible_SMC[20]
mw = N_visible_MW[20]
lmc / (lmc+smc+mw)
smc / (lmc+smc+mw)
mw / (lmc+smc+mw)
# SMC =  - 18%  19
# LMC =  - 32%  32
# MW =  - 50% 49



# GK16 with steeper slope of 3.31
low=2; high=7.5
pg = -1.0
min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)
model = am.GK16_grow(hostSHMF=True,reionization=True, plotgamma=pg)
N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
model = am.GK16_grow(lmcSHMF=True,reionization=True, plotgamma=pg)
N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)
N_visible_LMC, N_std_LMC, samples_LMC = getLMC(model,min_lums)
min_lums[20]
# normal masses > 10^4:
lmc = N_visible_LMC[20]
smc = N_visible_SMC[20]
mw = N_visible_MW[20]
lmc / (lmc+smc+mw)
smc / (lmc+smc+mw)
mw / (lmc+smc+mw)
# SMC =  - 23%  
# LMC =  - 32%  
# MW =  - 45%  





# GK16 with reionization at 8.3
low=2; high=7.5
pg = -0.2
z=8
min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)
model = am.GK16_grow(hostSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)

N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
model = am.GK16_grow(lmcSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt, plotgamma=pg)
N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)
N_visible_LMC, N_std_LMC, samples_LMC = getLMC(model,min_lums)
min_lums[20]
# normal masses > 10^4:
lmc = N_visible_LMC[20]
smc = N_visible_SMC[20]
mw = N_visible_MW[20]
lmc / (lmc+smc+mw)
smc / (lmc+smc+mw)
mw / (lmc+smc+mw)
# SMC =   21 %  
# LMC =   34 %  
# MW =    45 %  



# GK16 with bend introduced
low=2; high=7.5
pg = -0.7
min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)
model = am.GK16_grow(hostSHMF=True,reionization=True, bent=True, plotgamma=pg)

N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
model = am.GK16_grow(lmcSHMF=True,reionization=True, bent=True, plotgamma=pg)
N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)
N_visible_LMC, N_std_LMC, samples_LMC = getLMC(model,min_lums)
min_lums[20]
# normal masses > 10^4:
lmc = N_visible_LMC[20]
smc = N_visible_SMC[20]
mw = N_visible_MW[20]
lmc / (lmc+smc+mw)
smc / (lmc+smc+mw)
mw / (lmc+smc+mw)
# SMC = 24  % 
# LMC = 31  %  
# MW =  45  % 





# Moster
low=2; high=7.5
min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)
model = am.Moster(hostSHMF=True,reionization=True)

N_visible_MW, Nstd_MW, samples_MW = getMW(model,min_lums)
model = am.Moster(lmcSHMF=True,reionization=True)
N_visible_SMC, Nstd_SMC, samples_SMC = getSMC(model,min_lums)
N_visible_LMC, N_std_LMC, samples_LMC = getLMC(model,min_lums)
min_lums[20]
# normal masses > 10^4:
lmc = N_visible_LMC[20]
smc = N_visible_SMC[20]
mw = N_visible_MW[20]
lmc / (lmc+smc+mw)
smc / (lmc+smc+mw)
mw / (lmc+smc+mw)
# SMC =   19%  
# LMC =   32%  
# MW =    49% 





# For GK16 Model
SMC = 1.645250 
LMC = 2.761750 
MW = 3.5354999 

MW_r= 0.0720138 
LMC_r =0.3660271951182939
SMC_r = 0.38563
(LMC/LMC_r + SMC/SMC_r) / (MW/MW_r + MW+LMC+SMC)  = 20%
(LMC/LMC_r + SMC/SMC_r) / (MW/MW_r)
# 24% not including the LMC vicinity as extras.


# For Moster Model
LMC = 2.6360000
MW = 3.97825000000
SMC = 1.5325

MW_r= 0.0698049317676
LMC_r = 0.4209565996022283
SMC_r = 0.462038583656
(LMC/LMC_r + SMC/SMC_r) / (MW/MW_r + MW+LMC+SMC) = 15%
(LMC/LMC_r + SMC/SMC_r) / (MW/MW_r) = 17%
