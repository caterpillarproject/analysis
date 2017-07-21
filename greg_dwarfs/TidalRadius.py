import numpy as np
from scipy import interpolate


G = 4.49225e-12 # kpc^3/Msun/myr^2                                                                   
rho_o = 3588168 # Msun/kpc^3                                                                         
r_s = 29.767
kpc_myr_2_km_s = 978.462
cm2_g_2_kpc2_msun = 2.088384e-10
rho_crit = 147.9 # Msun/kpc^3 

class halo(object):
    def __init__(self,mvir, mvir_host=1.4e12):
        delta = 104
        c=80.08*mvir**(-1/11.)
        self.mvir = mvir
        self.M = mvir*(1+c)**2/c**2
        self.rvir = (3*mvir/(4*np.pi*delta*rho_crit))**(1/3.)
        self.a = self.rvir/c
        #self.M = float(M)  # in solar masses                                                        
        #self.a = 1.39*M**(1./9) #float(a)  # in kpc                                                 
        #mvir_host = 1.4e12  # 2.165156e12
        self.mvir_host = mvir_host
        c_host = 80.08*mvir_host**(-1/11.)
        self.M_host = mvir_host*(1+c_host)**2/c_host**2
        self.rvir_host = (3*mvir_host/(4*np.pi*delta*rho_crit))**(1/3.)
        self.a_host = self.rvir_host/c_host

        maxR = np.log(self.rvir_host*1.5)
        lnR = np.linspace(0,maxR,400)
        lnM = np.log(self.Menc_host(np.e**lnR))
        self.tck_dlnMdlnR = interpolate.splrep(lnR,lnM)


    # r in kpc, return in Msun                                                                       
    def Menc(self,r):
        return self.M*r**2/(r+self.a)**2

    def Menc_host(self,r):
        return self.M_host*r**2/(r+self.a_host)**2
    # r in kpc                                                                                       
    # return: Msun/kpc^3                                                                             
    def rho(self,r):
        return self.M*self.a/(2*np.pi*r)/(r+self.a)**3

    def tidal_radius(self,r_peri):
        #from scipy.optimize import fsolve                                                           
        #import time               
        #st = time.time()      
        #def f(x):    
        #    return self.Menc(x)/(4/3.*np.pi*x**3) - 3*self.Menc_host(r_peri)/ (4/3.*np.pi*r_peri**3\)            
        #sol1 = fsolve(f,0.4)[0]                                                                     
        #et1 = time.time()-st                                                                        
        ## analytic solution                                                                         
        #st2=time.time()                                                                             
        b=self.M*r_peri*(r_peri+self.a_host)**2/(3*self.M_host)
        a=self.a
        const = (3*np.sqrt(3)*np.sqrt(4*a**3*b+27*b**2)+2*a**3+27*b)**(1/3.)
        sol2 = 1/3.*(const/2**(1/3.) + 2**(1/3.)*a**2/const) - 2*a/3.
        return sol2

    def approx_tidal_radius(self,r_peri):
        return .7*(self.M/self.M_host)**(1/3.) * r_peri

    def tidal_radius_core(self,r_peri,core_radius):
        from scipy.optimize import fsolve
        tck = self.remake_cored_profile(core_radius)
        def f(x):
            return interpolate.splev(x,tck)/(4/3.*np.pi*x**3) - 3*self.Menc_host(r_peri)/ (4/3.*np.pi*r_peri**3)
        sol1 = fsolve(f,0.4)[0]
        return sol1


    # add in disk potential. Disk extends to 9.7 kpc, then density scales as 1/r^3                   
    # Make the host halo NFW density 90%, to make total mass still 1e12 Msun.            
    def tidal_radius_disk(self,r_peri):
        a=6.5;b=0.25;Md=1.e11 # disk is 1/10 host virial mass                                        
        Md = 0.1 * self.mvir_host
        if r_peri < 9.7:
            rho_disk = 3*a*Md/(r_peri*8*np.pi*(a+b)**3)
        else:
            rho_disk = 3*Md/(4*np.pi*r_peri**3)
        from scipy.optimize import fsolve
        def f(x):
            return self.Menc(x)/(4/3.*np.pi*x**3) - 3* (0.9*self.Menc_host(r_peri)/ (4/3.*np.pi*r_peri**3) + rho_disk)
        sol1 = fsolve(f,0.2)[0]
        return sol1



    def tidal_radius_core_disk(self,r_peri,core_radius):
        a=6.5;b=0.25;Md=1.e11 # disk is 1/10 host virial mass                                        
        Md= 0.1 * self.mvir_host
        if r_peri < 9.7:
            rho_disk = 3*a*Md/(r_peri*8*np.pi*(a+b)**3)
        else:
            rho_disk = 3*Md/(4*np.pi*r_peri**3)

        from scipy.optimize import fsolve
        tck = self.remake_cored_profile(core_radius)
        def f(x):
            return interpolate.splev(x,tck)/(4/3.*np.pi*x**3) - 3*(0.9*self.Menc_host(r_peri)/ (4/3.*np.pi*r_peri**3) + rho_disk)
        #print (3*self.Menc_host(r_peri)/ (4/3.*np.pi*r_peri**3)), 'host dens'                       
        #print 3*a*Md/(r_peri*8*np.pi*(a+b)**3), 'disk dens'                                         
        sol1 = fsolve(f,0.02)[0]
        return sol1

    def dlnMdlnR(self,R):
        return interpolate.splev(np.log(R),self.tck_dlnMdlnR,der=1)

    def tidal_radius_BT(self,R):
        # R is orbit radius                                                                          
	# use binney and tremaine version from 8.107 and 8.108                                       
	# find dln M/dlnR                                                                            
        from scipy.optimize import fsolve
        f = (1-1/3.*self.dlnMdlnR(R))**(-1/3.)
	const = (1/(3*self.Menc_host(R)))**(1/3.)*R
        def f(x):
            return x - const*self.Menc(x)**(1/3.)
        sol1 = fsolve(f,3)[0]
        return sol1




## this is to use raw particle data for the host profile
def tidal_radius(v,rsid):
    rarr,mltr,tck_mltr=get_mltr(v,rsid)
    from scipy.optimize import fsolve
    import time
    st = time.time()
    h = insc.halo(mvir=10**9)
    r_peri=30
    def f(x):
        return interpolate.splev(x,tck_mltr)/(4/3.*np.pi*x**3) - 3*h.Menc_host(r_peri)/ (4/3.*np.pi*r_peri**3)
    sol1 = fsolve(f,0.4)[0]
    et1 = time.time()-st
    return sol1

"""

# make a plot of tidal radius vs distance to the MW host.

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import abundance_matching as am


LMC_stellar_mass = 2.6*10**9
SMC_stellar_mass = 7.1*10**8
model = am.GK16_grow(reionization=True)
lmc_mvir = model.mstar_to_mhalo(LMC_stellar_mass)
MW_mvir = 1.4e12

MW_LMC = halo(lmc_mvir, MW_mvir) 
MW_LMC.tidal_radius_disk(50)
MW_LMC.tidal_radius(50)


dwarf_mvir = model.mstar_to_mhalo(10**4)
radius = np.arange(0,70,2)

MW_dwarf = halo(dwarf_mvir, MW_mvir) 
r_tidal_MW = [MW_dwarf.tidal_radius_disk(r) for r in radius]

lmc_dwarf = halo(dwarf_mvir, lmc_mvir) 
r_tidal_lmc = [lmc_dwarf.tidal_radius_disk(r) for r in radius]


plt.plot(radius, r_tidal_MW, lw=3, label='LMC Host')
plt.plot(radius, r_tidal_lmc, lw=3,label='MW Host')

"""
