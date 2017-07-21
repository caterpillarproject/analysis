import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import abundance_matching as am
from PlotParams import *


MW_mass = 1.6e12 # from GK14 paper
min_lum = 4.5e5 # from GK16 paper 
h=0.67
G = 4.325e-6 # kpc * km^2 /s^2 /Msun # in 
rho_crit = 3*(100*h/1000.)**2/(8*np.pi*G)
del_crit = 103.86
MW_rvir = (3*MW_mass / (4*np.pi*rho_crit*del_crit))**(1/3.)
rvir=300 # form garrison kimmel 2016
SHMF_K = am.radial_dependence(rvir/MW_rvir)
SHMF_alpha = 1.90





#for model in [am.Moster(hostSHMF=True), am.GarrisonKimmel(hostSHMF=True), am.Brook(hostSHMF=True), am.GarrisonKimmel16(hostSHMF=True), am.Sawala(hostSHMF=True)]:
model = am.GarrisonKimmel(hostSHMF=True, hostAlpha=SHMF_alpha, hostK=SHMF_K)
N_visible,_ = model.ngreater(MW_mass, min_lum)
print N_visible, 'at 300 kpc'

model = am.GarrisonKimmel(hostSHMF=False)
N_visible,_ = model.ngreater(MW_mass, min_lum)
print N_visible, 'at rvir, default SHMF'




model = am.GarrisonKimmel(hostSHMF=True, hostAlpha=1.895, hostK=.00102)
N_visible,_ = model.ngreater(MW_mass, min_lum)
print N_visible, 'host halo'
model = am.GarrisonKimmel(hostSHMF=True, hostAlpha=1.90, hostK=.001122)
N_visible,_ = model.ngreater(MW_mass, min_lum)
print N_visible, 'host halo radial'
model = am.GarrisonKimmel(hostSHMF=True, hostAlpha=1.867, hostK=.000703)
N_visible,_ = model.ngreater(MW_mass, min_lum)
print N_visible, 'field halo'
model = am.GarrisonKimmel(hostSHMF=True, hostAlpha=1.90, hostK=.001314)
N_visible,_ = model.ngreater(MW_mass, min_lum)
print N_visible, 'field halo radial'




# N_visible should be 9 to 11.
