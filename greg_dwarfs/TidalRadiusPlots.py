import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import abundance_matching as am
from TidalRadius import *

LMC_stellar_mass = 2.6*10**9
SMC_stellar_mass = 7.1*10**8
model = am.GK16_grow(reionization=True)
lmc_mvir = model.mstar_to_mhalo(LMC_stellar_mass)
MW_mvir = 1.4e12

#MW_LMC = halo(lmc_mvir, MW_mvir) 
#MW_LMC.tidal_radius_disk(50)
#MW_LMC.tidal_radius(50)

dwarf_mvir = model.mstar_to_mhalo(10**4)
radius = np.arange(0,70,2)

MW_dwarf = halo(dwarf_mvir, MW_mvir) 
r_tidal_MW = [MW_dwarf.tidal_radius_disk(r) for r in radius]
r_tidal_MW_nodisk = [MW_dwarf.tidal_radius(r) for r in radius]

lmc_dwarf = halo(dwarf_mvir, lmc_mvir) 
r_tidal_lmc = [lmc_dwarf.tidal_radius_disk(r) for r in radius]


plt.plot(radius, r_tidal_MW, lw=3, label='MW Host')
plt.plot(radius, r_tidal_MW_nodisk, lw=3, label='MW Host No Disk')
plt.plot(radius, r_tidal_lmc, lw=3,label='LMC Host')
plt.legend(loc='upper left', frameon=False)
plt.xlabel('Distance to Host [kpc]')
plt.ylabel('Tidal Radius [kpc]')
plt.savefig('LMCPlots/TidalRadiusMWandLMC')
plt.close()

