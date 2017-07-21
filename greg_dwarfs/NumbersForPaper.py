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



## 1)  Radius of the SMC in Brook model?

model = am.Brook(reionization=True)
Mvir = model.mstar_to_mhalo(SMC_stellar_mass)

G = 4.157e-39  # kpc^3/Msun/s^2
H = dm.Hubble(a=1) # in 1/s
rho_crit = (3*H**2)/(8*np.pi*G) # in Msun/kpc^3
rvir = (Mvir*3/(4*np.pi*103.86*rho_crit))**(1/3.) # in kpc
# solution is 108 kpc



# Find Mvir in all models
model = am.GarrisonKimmel(reionization=True)
mvir = model.mstar_to_mhalo(LMC_stellar_mass)
rvir = (mvir*3/(4*np.pi*103.86*rho_crit))**(1/3.) # in kpc

model = am.GK16_grow(reionization=True)
mvir = model.mstar_to_mhalo(LMC_stellar_mass)
rvir = (mvir*3/(4*np.pi*103.86*rho_crit))**(1/3.) # in kpc

model = am.Moster(reionization=True)
mvir = model.mstar_to_mhalo(LMC_stellar_mass)
mvir, rvir = dm.convert_Mhalo_d_d(mvir, 200, 103.86)
#rvir = (mvir*3/(4*np.pi*103.86*rho_crit))**(1/3.) # in kpc

model = am.Brook(reionization=True)
mvir = model.mstar_to_mhalo(LMC_stellar_mass)
mvir, rvir = dm.convert_Mhalo_d_d(mvir, 350, 103.86)
#rvir = (mvir*3/(4*np.pi*103.86*rho_crit))**(1/3.) # in kpc

model = am.Behroozi(reionization=True)
mvir = model.mstar_to_mhalo(LMC_stellar_mass)
rvir = (mvir*3/(4*np.pi*103.86*rho_crit))**(1/3.) # in kpc


# SMC values
model = am.GarrisonKimmel(reionization=True)
mvir = model.mstar_to_mhalo(SMC_stellar_mass)
rvir = (mvir*3/(4*np.pi*103.86*rho_crit))**(1/3.) # in kpc

model = am.GK16_grow(reionization=True)
mvir = model.mstar_to_mhalo(SMC_stellar_mass)
rvir = (mvir*3/(4*np.pi*103.86*rho_crit))**(1/3.) # in kpc

model = am.Moster(reionization=True)
mvir = model.mstar_to_mhalo(SMC_stellar_mass)
mvir, rvir = dm.convert_Mhalo_d_d(mvir, 200, 103.86)
#rvir = (mvir*3/(4*np.pi*103.86*rho_crit))**(1/3.) # in kpc

model = am.Brook(reionization=True)
mvir = model.mstar_to_mhalo(SMC_stellar_mass)
mvir, rvir = dm.convert_Mhalo_d_d(mvir, 350, 103.86)
#rvir = (mvir*3/(4*np.pi*103.86*rho_crit))**(1/3.) # in kpc

model = am.Behroozi(reionization=True)
mvir = model.mstar_to_mhalo(SMC_stellar_mass)
rvir = (mvir*3/(4*np.pi*103.86*rho_crit))**(1/3.) # in kpc



## 2) convert 2*10^11 and 2.5*10^11 from Mvir to M200
print dm.convert_Mhalo_d_d(2*10**11, 200, 103.86)
print dm.convert_Mhalo_d_d(2.5*10**11, 200, 103.86)



## number of satellites above 10^5 in the LMC and SMC.


##Also get the probability




### PROBLEM!! FOR MW MASS, use the same mass within 300 kpc.
# Convert Mvir to M350, M200, etc. for that plot.
### MORE PROBLEMS!! WHEN USING RADIUS in my plots, must convert from
# mass to Mvir and then to Rvir.  MIGHT BE A PROBLEM IN OLD PAPER!!
## yes, fixed this now in the Table generating functions for old paper and current paper.
# Both need to be run still.
# Fixed in current figure 3


# get what fraction have 0 > 10^5 to report as a percentage.
#model = am.GK16_grow(lmcSHMF=True,reionization=True)
#model = am.GK16_grow(reionization=True, lmcSHMF=True, catreion=True,z=9, vmax_ach=vmax_ach, vmax_filt=vmax_filt)

model = am.Brook(lmcSHMF=True,reionization=True)
halo_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1.0)
#N_visible, N_std = model.ngreater(halo_mass,10**5)
samples_LMC = model.generate_stellar_mass_samples(halo_mass, factor=1)
N_visible, N_std, p20, p80 = ambase.Nvisible_From_Ngreater(10**5, samples_LMC)
print N_visible+N_std, N_visible, N_visible-N_std  # 5.66 and 1.69  with 6,000 examples. mean is 3.68
 # 3.1332956827 1.764 0.394704317298 for Brook

Ngreater = np.array([np.sum(lum>10**5) for lum in samples_LMC])
numlow = np.sum(Ngreater==0)
print numlow, len(Ngreater)
print float(numlow)/len(samples_LMC)*100  # 2.0 % chance



halo_mass = model.mstar_to_mhalo(SMC_stellar_mass,a=1.0)
#N_visible, N_std = model.ngreater(halo_mass,10**5)
samples_SMC = model.generate_stellar_mass_samples(halo_mass, factor=1)
N_visible, N_std, p20, p80 = ambase.Nvisible_From_Ngreater(10**5, samples_SMC)
print N_visible+N_std, N_visible, N_visible-N_std  # 3.65 and 0.67 with 6,000 examples. mean is 2.16
Ngreater = np.array([np.sum(lum>10**5) for lum in samples_SMC])
numlow = np.sum(Ngreater==0)
print numlow, len(Ngreater)
print float(numlow)/len(samples_SMC)*100  # 10.0 % chance




## Get slope of Mhalo and Stellar mass
#model = am.GarrisonKimmel(reionization=True, lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt)
model = am.GK16_grow(reionization=True, lmcSHMF=True)
#model = am.Moster(reionization=True, lmcSHMF=True)

loglmc = np.log10(SMC_stellar_mass)
stellar_masses = 10**np.array([loglmc-0.1, loglmc, loglmc+0.1])  # 10**np.linspace(8.5,10,14)
halo_mass = np.log10(model.mstar_to_mhalo(stellar_masses))
print (halo_mass[2]-halo_mass[1])/.1
print (halo_mass[1]-halo_mass[0])/.1
print (halo_mass[2]-halo_mass[0])/.2
# slope is 0.40 for both the LMC and the SMC
# 2**0.4 should equal the values computed just below, 1.3, which it is.

## doubling stellar mass of the LMC results in what increase of ngreater?

loglmc = np.log10(LMC_stellar_mass)

stellar_masses = np.array([LMC_stellar_mass*0.5, LMC_stellar_mass, 2*LMC_stellar_mass])
halo_masses = model.mstar_to_mhalo(stellar_masses) 
N_visible, N_std = model.ngreater(halo_masses[0],np.array([10**3, 10**4, 10**5]))
N_visible2, N_std2 = model.ngreater(halo_masses[1],np.array([10**3, 10**4, 10**5]))
N_visible3, N_std3 = model.ngreater(halo_masses[2],np.array([10**3, 10**4, 10**5]))

print N_visible3/N_visible2
print N_visible2/N_visible
# Answer seems to be about 1.3

# now within 100 kpc
rvir = model.mass_to_rvir(halo_masses)
factor =  rd.get_normed_radial_abundance(50 / rvir)
print (N_visible3*factor[2])/(N_visible2*factor[1])
print (N_visible2*factor[1])/(N_visible*factor[0])


#slope = (N_visible3-N_visible)/0.2
#slope/N_visible[1]









### Line of Sight radial distribution values
model = am.Brook(reionization=True, lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt)
dm_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1)
dummy, lmc_rvir = dm.convert_Mhalo_d_d(dm_mass, model.delta_crit, 103.86)
ngreater = model.ngreater(dm_mass,min_lum=10**5, percentile=False)[0]
print rd.get_normed_los_better(50./lmc_rvir, Z=1.5)*ngreater  # trying to guess distance to get this to 1
# 50 is just about perfect.


print rd.get_normed_los_better(0.5, Z=1.5)
print rd.get_normed_los_better(0.5, Z=1.5) / rd.get_normed_los_better(1.0, Z=1.5)
# 69 % contained within half of the virial radius. Round to 70%


########### SECTION 3.2 CALCULATIONS #############################
model = am.GarrisonKimmel(reionization=True, lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt)
dm_mass = model.mstar_to_mhalo(LMC_stellar_mass)
ngreater = model.ngreater(dm_mass,min_lum=10**5, percentile=False)[0]
dummy, lmc_rvir = dm.convert_Mhalo_d_d(dm_mass, model.delta_crit, 103.86)
factor1=rd.get_normed_radial_abundance(50 / lmc_rvir) 
in50kpc = ngreater*factor1
print in50kpc, 'orig stellar mass'

dm_mass = model.mstar_to_mhalo(0.7*LMC_stellar_mass)
ngreater = model.ngreater(dm_mass,min_lum=10**5, percentile=False)[0]
dummy, lmc_rvir = dm.convert_Mhalo_d_d(dm_mass, model.delta_crit, 103.86)
factor2 = rd.get_normed_radial_abundance(50 / lmc_rvir) 
in50kpc2 = ngreater *factor2
print in50kpc2, 'with lower stellar mass'
print 1 - in50kpc2/in50kpc
print factor2/factor1


# 30% decrease in stellar mass of host increases satellite abundances by how much?
print (1 - 0.7**0.4)*100 # percent decrease
print (1 - ( (0.7**0.4) * factor2/factor1 ))*100 # percent decrease when limiting to a fixed volume



### How many satellites come from LMC, SMC, and MW?
import LMC50kpc as lmc50
min_lums, N_visible, N_visible_SMC, N_visible_MW, N_visible_full, up, down,all_samples = lmc50.Baseline(True)
min_lums[20] # should be 10**4
N_visible_full[20]
up[20] - N_visible_full[20]
N_visible_full[20] - down[20]

Ngreater = np.array([np.sum(lum>10**4) for lum in all_samples])
numlow = np.sum(Ngreater==0)
num1 = np.sum(Ngreater==1)
print (float(numlow + num1))/len(Ngreater)*100  # this is the % chance of getting 0 or 1 > 10^4 Msun. answer is 0.3%
print (float(numlow))/len(Ngreater)*100  # this is the % chance of getting 0 > 10^4 Msun. 
# GK16 chance of 0:  .0575% with 40,000 iterations, then .045%, then .03%  average the 3 get .044
# Brook chance of 0: 1.4975% with 40,000 iterations. then 1.58     7.5% of 0 or 1.

## percent chance after the tidal stripping shift
Ngreater = np.array([np.sum(lum> (10**4 / .05) ) for lum in all_samples])
numlow = np.sum(Ngreater==0)
num1 = np.sum(Ngreater==1)
print (float(numlow + num1))/len(Ngreater)*100  # this is the % chance of getting 0 or 1 > 10^4 Msun
# answer after tidal stripping of 95% is now 20%

# Fraction in LMC, in SMC, and in MW
# make sure lmc50.Baseline() does what I want it to.

# above 10^2? above 10^3? above 10^4? Does it change?

# Now make reionization later. Does it change? (it shouldn't)

# Make slope steeper - does it change?  (it shouldn't)






# looking at 10^2 numbers:  4.2705, 2.7355, 5.484 , LMC, SMC, MW

total = 4.2705 + 2.7355 + 5.484
print 4.2705/total,  2.7355/total, 5.484/total
low=2; high=7.5
min_lums = 10**np.linspace(low,high, (high-low)/.10 +1)  


### Mean distance of satellite subsets to the MW, LMC and SMC
## Want to compare subhalos like the largest 11 in the MW
## to the subhalos within 100 kpc of MW, 50 kpc of LMC, and ~50 kpc of SMC

galx = dm.load_nearby_gx()
mass_cut  = galx['mstar']>2.2*10**5
argsort = np.argsort(galx[mass_cut]['mstar'])[::-1]
# now also remove the LMC and SMC
classical_dwarfs = galx[mass_cut][argsort][3:]
mean_MW_dist = np.mean(classical_dwarfs['dist_GC'])    # 122 to the MW
mean_LMC_dist = np.mean(classical_dwarfs['dist_LMC'])  # 129 kpc to the LMC
mean_SMC_dist = np.mean(classical_dwarfs['dist_SMC'])  # 138 kpc to the SMC


near_LMCmask = galx['dist_LMC'] < 50
masssort = np.argsort(galx[near_LMCmask]['mstar'])
namesort = np.argsort(galx[near_LMCmask]['Name'])
near_LMC = galx[near_LMCmask][masssort][0:-2]


mean_MW_dist = np.mean(near_LMC['dist_GC'])    # 53 to the MW
mean_LMC_dist = np.mean(near_LMC['dist_LMC'])  # 36 kpc to the LMC
mean_SMC_dist = np.mean(near_LMC['dist_SMC'])  # 36 kpc to the SMC



## what is mpeak if mstar is 10^6
model = am.GarrisonKimmel(reionization=True, lmcSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt)
dm_mass = model.mstar_to_mhalo(10**6)
dm_mass = model.mstar_to_mhalo(10**4)



### MW Radial Figure Calculations:
model = am.GK16_grow(reionization=True, hostSHMF=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt)
dm_mass = 1.4*10**12 
rvir = model.mass_to_rvir(dm_mass)
radius = 300 # focus on full 300

factor =  rd.get_normed_radial_abundance(radius / rvir)
ngreater, N_std, p20, p80 = model.ngreater(dm_mass,min_lum=10**3, percentile=True)
r_sigma = rd.get_normed_radial_std(radius / rvir)
sig2 = ngreater*r_sigma
N = factor*ngreater
sig1 = np.sqrt(N + 0.14**2 * N**2)
full_sigma = np.sqrt( sig1**2 + sig2**2)
print factor*ngreater+full_sigma, factor*ngreater, factor*ngreater-full_sigma
# 121.768578586 103.993323198 86.2180678108


### How many LMC and SMC satellites M* > 10^3
model = am.GK16_grow(lmcSHMF=True,reionization=True, catreion=True,z=z, vmax_ach=vmax_ach, vmax_filt=vmax_filt)

halo_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1.0)
N_visible, N_std = model.ngreater(halo_mass,10**3)
print N_visible+N_std, N_visible, N_visible-N_std  # 6.12 and 1.95 with 100,000 examples. mean is 4.0
# FOR LMC:  20.9, 16.4, 11.8


halo_mass = model.mstar_to_mhalo(SMC_stellar_mass,a=1.0)
N_visible, N_std = model.ngreater(halo_mass,10**3)
print N_visible+N_std, N_visible, N_visible-N_std  # 3.9 and 0.78 with 100,000 examples. mean is 2.35
# FOR SMC: 12.9, 9.6, 6.2



#### Ratio of min pericenter and current distance of dwarfs surviving reionization ###




# uncertainty in radial distribution for LMC
model = am.GK16_grow(lmcSHMF=True,reionization=True)
halo_mass = model.mstar_to_mhalo(LMC_stellar_mass,a=1.0)
lmc_rvir = model.mass_to_rvir(halo_mass)
r_sigma = rd.get_normed_radial_std(50 / lmc_rvir)
factor =  rd.get_normed_radial_abundance(50 / lmc_rvir)
print r_sigma / factor

model = am.GK16_grow(lmcSHMF=True,reionization=True)
halo_mass = model.mstar_to_mhalo(SMC_stellar_mass,a=1.0)
smc_rvir = model.mass_to_rvir(halo_mass)
r_sigma = rd.get_normed_radial_std( 26 / smc_rvir)
factor =  rd.get_normed_radial_abundance(50 / smc_rvir)
print r_sigma / factor

halo_mass = 1.4e12
mw_rvir =  model.mass_to_rvir(halo_mass)
r_sigma = rd.get_normed_radial_std(50 / mw_rvir)
factor =  rd.get_normed_radial_abundance(50 / mw_rvir)
print r_sigma / factor



#### How many 10^4 satellits are within 100 kpc, and how many shoould be in lmc subvolume?
galx = dm.load_nearby_gx()
mass_cut  = galx['mstar']>10**4
MWmask = galx['dist_GC'] < 100
N = np.sum(MWmask & mass_cut) - 2
halo_mass = 1.4e12
model = am.GK16_grow(lmcSHMF=True,reionization=True)
rvir = model.mass_to_rvir(halo_mass)
factor =  rd.get_normed_radial_abundance(100. / rvir)
ratio_MW = 0.0720138 
print ratio_MW / factor
print N
print N * ratio_MW / factor
