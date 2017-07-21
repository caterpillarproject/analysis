import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import haloutils as htils
from scipy import stats
import methods
import MTanalysis_field as mtaf
import DwarfMethods as dm


hpaths = dm.get_hpaths(field=True)
h0=0.67
for hpath in hpaths:
    FieldData = mtaf.FieldHaloSubstructureFirstPass()
    data = FieldData.read(hpath)

    FieldDataAE = mtaf.AllExtantFieldData()
    dataAE = FieldDataAE.read(hpath)

    hostmass = data['host_mgravz0']/h0
    max_masses = data['max_mass']/h0
    infall_m200 = dataAE['infall_mass200']/h0
    maxmassNFW350 = dataAE['max_mass350NFW']/h0

    if len(hostmass) > 0:
        #plt.scatter(maxmassNFW350, infall_m200/maxmassNFW350)
        #plt.scatter(max_masses, infall_m200/max_masses)
        plt.scatter(max_masses, max_masses/infall_m200)



plt.ylabel('m200_infall / other mass')
plt.xlabel('other mass')
plt.xscale('log')
plt.ylim((0,100))
#plt.yscale('log')
#plt.legend()
plt.savefig('convert_mass_defs')
plt.close()




"""


I need to modify my code.
1st, extend the semi-analytic data to all mass definitions. This is too hard!
2nd, modify my code to get the peak vmax and peak Tvir of the satellite before the redshift.



Currently I have it at the redshift, but if it was accreted by z=10, it may be shrinking already.
Use the thresholds for vmax, for 10^6 Msun M200 and the Tvir proposed by other authors to get the lum_frac plots


data = dm.get_extant_data(hpath,False)

E = mtadd.ExtantDataReionization()
data_reion = E.read(hpath)
maskR = data_reion['depth_reion'] == 0 
data_reion = data_reion[maskR]
return pandas.concat((dataE,data_reion),axis=1)


import MTaddition as mtadd
Must add these things to MTaddition



"""
















