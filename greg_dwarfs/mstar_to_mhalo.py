
"""
Need to get a density metric for every data point.
metric includes dn/dm(Mhalo) from seth-tormen.
For a fixed Mstar, choose the Mhalo
then also need gauss(mstar, mean=mstar(mhalo), var = sigma(Mhalo))
"""
import abundance_matching as am
import random
import AnalyticMF as MF
import numpy as np
from scipy import interpolate

# me: its on basketball now, essentially the NBA all star team vs random guys from other countries.
# sara: I'd be so nervous to take utensils from there


def gauss(masses, sigma, mean, model):
    mstar_dist = np.log10(mean) - np.log10(model.getStellarMass(masses))
    return 1./np.sqrt(2*np.pi * sigma**2)  * np.e**(-( mstar_dist**2 )/(2*sigma**2))

def mstar_to_mhalo(model, Mstar):
    Mhalo = model.stellar_to_halo_mass(Mstar, a=1)
    halo_masses = np.logspace(np.log10(Mhalo)-.8,  np.log10(Mhalo)+0.3, 30)
    sigmas = model.sigma(halo_masses)
    #print sigmas
    weight1 = gauss(halo_masses, sigmas, Mstar, model)

    h0 = 0.6711; Om = 0.31; sig8 = 0.83; ns = 0.96

    Theoretical = MF.AnalyticMF(h=h0, Omega_m=Om,sigma8=sig8, n=ns)
    #STdndlogM = Theoretical.STdNdLogM(Mhalo, z=1.)
    weight2 = Theoretical.STdNdLogM(halo_masses,z=1)

    lum_tck, _, _ = reion.get_lum_fraction_mpeak()
    # shift lum frac by 0.1 dex
    halo_masses = 10**(np.log10(halo_masses)-0.1)
    lum_chances = reion.frac_luminous(halo_masses, self.lum_tck) # double check this is between 0 and 1
    weight3 = lum_chances
    #print weight3

    values = weight1*weight2*weight3
    values = values / np.sum(values)
    cum = np.cumsum(values)
    #print np.log10(cum)


    tck = interpolate.splrep(np.log10(cum)[3:-3] , np.log10(halo_masses)[3:-3], k=1)
    mass = 10**interpolate.splev(np.log10(0.5), tck)  # THIS IS THE FINAL ANSWER
    print mass/Mhalo, 'ratio of new mass to 0 scatter mass'
    if mass/Mhalo < 0.5 or mass/Mhalo > 1.0:
        raise Exception('mstar to mhalo failure!!!!')
    return mass



#model = am.GK16_grow()
#Mstar = 10**8
#print mstar_to_mhalo(model, Mstar)

# now want to find the median of values. or the max??
# plot the distribution. Is it wide enough to go to 0?
# If yes, I can normalize it

import numpy as np
import abundance_matching as am
from PlotParams import *
import DwarfMethods as dm
import RadialDependence as rd

stellar_masses = np.logspace(5.1,8.45,25)
#stellar_masses = np.array([10**5, 10**6])
#model = am.GK16_grow(reionization=True)
#model = am.GarrisonKimmel(reionization=True)
model = am.Moster(reionization=True)
#model = am.Brook(reionization=True)
halo_masses = model.stellar_to_halo_mass(stellar_masses)
hmass = model.mstar_to_mhalo(stellar_masses)



"""
# testing out some weirdness
model = am.Moster(reionization=True)
mstar = 270*10**6
halo_mass = model.stellar_to_halo_mass(mstar ,a=1.0)
ng = model.ngreater(halo_mass,10**3)[0]
print ng, 'should be 3.1'

halo_mass2 = model.mstar_to_mhalo(mstar,a=1)
print halo_mass2/halo_mass, 'ratio of halo masses'
ng2 = model.ngreater(halo_mass2,10**3)[0]
print ng2, 'should be', ng*halo_mass2/halo_mass



mstar = 190*10**6
halo_mass = model.stellar_to_halo_mass(mstar ,a=1.0)
ng = model.ngreater(halo_mass,10**4)[0]
print ng, 'should be 2.6'

halo_mass2 = model.mstar_to_mhalo(mstar,a=1)
print halo_mass2/halo_mass, 'ratio of halo masses'
ng2 = model.ngreater(halo_mass2,10**4)[0]
print ng2, 'should be', ng*halo_mass2/halo_mass
"""








