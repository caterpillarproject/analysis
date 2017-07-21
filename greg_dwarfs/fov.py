import RadialDependence as rd
import numpy as np
import haloutils as htils
import DwarfMethods as dm
import abundance_matching as am


"""
# IC 5152
dist = 1700 # distance in km
mstar = 270 * 10**6
fov = 1.1*np.pi/180  #radius of field of view is 1.1 degrees. convert to radians
fov_radius = fov * dist
Z = 1

#model = am.Moster(reionization=True)
model = am.GK16_grow(reionization=True)
#model = am.GarrisonKimmel(reionization=True)
Mvir =  model.mstar_to_mhalo(mstar)
rvir = model.mass_to_rvir(Mvir)
R = fov_radius/rvir
frac = rd.get_normed_los(R,Z)
print frac   # this is 0.33

"""

# 3.72 for GK16
# 2.87 for moster
# 4.23 for GK14


dist = 1700 # distance in km
mstar = 270 * 10**6
fov = 1.1*np.pi/180  #radius of field of view is 1.1 degrees. convert to radians
fov_radius = fov * dist
Z = 1.5
for model, N in zip([am.Moster(reionization=True), am.GK16_grow(reionization=True), am.GarrisonKimmel(reionization=True), am.Brook(reionization=True)], [2.84,3.73,4.25,1.31]  ):
    Mvir =  model.mstar_to_mhalo(mstar)
    rvir = model.mass_to_rvir(Mvir)
    fov = 1.1*np.pi/180  #radius of field of view is 1.1 degrees. convert to radians
    fov_radius = fov * dist
    R = fov_radius/rvir
    #print R, 'R', model.label
    frac_old = rd.get_normed_los(R,Z)
    print frac_old, 'original frac'
    frac = rd.get_numerical_normed_los_better(R,Z)
    print frac, 'new frac', model.label
    print frac/frac_old, 'ratio to multiply by'

    #rd.get_normed_los()
    print frac   # this is 0.33
    print N*frac, model.label, 'answer'




#G = 4.157e-39  # kpc^3/Msun/s^2
#H = dm.Hubble(a=1) # in 1/s
#rho_crit = (3*H**2)/(8*np.pi*G) # in Msun/kpc^3
#rvir = (Mvir*3/(4*np.pi*103.86*rho_crit))**(1/3.) # in kpc
