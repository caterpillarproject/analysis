import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy import interpolate

class AnalyticMF:
    
    def __init__(self, h=0.7, Omega_m=.27, sigma8=.809, n=0.95, delta_c=1.686, omega_lambda = 0.73, omega_b=0.045):
        """
        @param h: hubble parameter.
        @param omega_m: current matter fraction of universe.
        @param sigma8: Value of variance of density field with spherical smoothing at 8 Mpc/h
        @param n: spectral scalar index for primordial power spectrum
        @param delta_c: amplitude of perturbation at collapse using linear theory.
        @param rho_c: critical density of universe in Msun/Mpc in an h = 1.0 univserse
        """
        rho_c = 2.77494623782e11*h**2
        self.h = h
        self.Omega_m = Omega_m
        self.n = n
        self.delta_c = delta_c
        self.rho_m = rho_c*Omega_m/self.h**2 #put in units of Msun/h /(Mpc/h)^3
        self.sigma8 = sigma8
        self.omega_lambda = omega_lambda
        self.omega_b = omega_b
        
    # Used for Eisenstein & Hu transfer function.
    # Translated from Wayne Hu's C-code.
    # Sets many global variables for use in TFmdm_onek_mpc()
    def TFmdm_set_cosm(self, redshift, omega_hdm =0.0, degen_hdm=1):
        """This routine takes cosmological parameters and a redshift and sets up
        all the internal scalar quantities needed to compute the transfer function.
        INPUT: Omega_matter -- Density of CDM, baryons, and massive neutrinos,
        in units of the critical density.
        omega_baryon -- Density of baryons, in units of critical.
        omega_hdm    -- Density of massive neutrinos, in units of critical
        degen_hdm    -- (Int) Number of degenerate massive neutrino species
        omega_lambda -- Cosmological constant 
        hubble       -- Hubble constant, in units of 100 km/s/Mpc
        redshift     -- The redshift at which to evaluate
        OUTPUT: Returns 0 if all is well, 1 if a warning was issued.  Otherwise,
        sets many global variables for use in TFmdm_onek_mpc()"""
        hubble = self.h
        self.omega_hdm = 0.0
        qwarn = 0
        self.theta_cmb = 2.728/2.7;	# Assuming T_cmb = 2.728 

        # Look for strange input
        if (self.omega_b<0.0):
            print "TFmdm_set_cosm(): Negative omega_baryon set to trace amount."
            qwarn = 1
        if (omega_hdm<0.0):
            print "TFmdm_set_cosm(): Negative omega_hdm set to trace amount."
            qwarn = 1
        if (hubble<=0.0):
            print "TFmdm_set_cosm(): Negative Hubble constant illegal."
            sys.exit()
        elif (hubble>2.0):
            print "TFmdm_set_cosm(): Hubble constant should be in units of 100 km/s/Mpc."
            qwarn = 1;
        if (redshift<=-1.0):
            print "TFmdm_set_cosm(): Redshift < -1 is illegal."
            sys.exit()
        elif (redshift>99.0):
            print "TFmdm_set_cosm(): Large redshift entered.  TF may be inaccurate."
            qwarn = 1
        if (degen_hdm<1):
            degen_hdm=1
        self.num_degen_hdm = float(degen_hdm)	
	# Have to save this for TFmdm_onek_mpc()
        # This routine would crash if baryons or neutrinos were zero, so don't allow that
        if (self.omega_b<=0):
            self.omega_b=1e-5
        if (self.omega_hdm<=0):
            self.omega_hdm=1e-5

        self.omega_curv = 1.0-self.Omega_m-self.omega_lambda
        self.omhh = self.Omega_m*hubble**2
        obhh = self.omega_b*hubble**2
        self.onhh = self.omega_hdm*hubble**2
        self.f_baryon = self.omega_b/self.Omega_m
        self.f_hdm = self.omega_hdm/self.Omega_m
        f_cdm = 1.0-self.f_baryon-self.f_hdm
        self.f_cb = f_cdm+self.f_baryon
        f_bnu = self.f_baryon+self.f_hdm

        # Compute the equality scale.
        z_equality = 25000.0*self.omhh/(self.theta_cmb**4) # Actually 1+z_eq 
        self.k_equality = 0.0746*self.omhh/(self.theta_cmb**2)

        # Compute the drag epoch and sound horizon
        z_drag_b1 = 0.313*(self.omhh**-0.419)*(1+0.607*(self.omhh**0.674))
        z_drag_b2 = 0.238*(self.omhh**0.223)
        z_drag = 1291*(self.omhh**0.251)/(1.0+0.659*(self.omhh**0.828))*(1.0+z_drag_b1*(obhh**z_drag_b2))
        y_drag = z_equality/(1.0+z_drag)

        self.sound_horizon_fit = 44.5*np.log(9.83/self.omhh)/np.sqrt(1.0+10.0*(obhh**0.75))

        # Set up for the free-streaming & infall growth function
        p_c = 0.25*(5.0-np.sqrt(1+24.0*f_cdm))
        self.p_cb = 0.25*(5.0-np.sqrt(1+24.0*self.f_cb))

        omega_denom = self.omega_lambda+((1.0+redshift)**2)*(self.omega_curv+self.Omega_m*(1.0+redshift))
        omega_lambda_z = self.omega_lambda/omega_denom
        omega_matter_z = self.Omega_m*((1.0+redshift)**2)*(1.0+redshift)/omega_denom;
        self.growth_k0 = z_equality/(1.0+redshift)*2.5*omega_matter_z/((omega_matter_z**(4.0/7.0))-omega_lambda_z+(1.0+omega_matter_z/2.0)*(1.0+omega_lambda_z/70.0))
    
        growth_to_z0 = z_equality*2.5*self.Omega_m/((self.Omega_m**(4.0/7.0))-self.omega_lambda + (1.0+self.Omega_m/2.0)*(1.0+self.omega_lambda/70.0))
        growth_to_z0 = self.growth_k0/growth_to_z0	
    
        # Compute small-scale suppression
        alpha_nu = f_cdm/self.f_cb*(5.0-2.*(p_c+self.p_cb))/(5.-4.*self.p_cb)*pow(1+y_drag,self.p_cb-p_c)*(1+f_bnu*(-0.553+0.126*f_bnu*f_bnu))/(1-0.193*np.sqrt(self.f_hdm*self.num_degen_hdm)+0.169*self.f_hdm*(self.num_degen_hdm**0.2))*(1+(p_c-self.p_cb)/2*(1+1/(3.-4.*p_c)/(7.-4.*self.p_cb))/(1+y_drag))
        self.alpha_gamma = np.sqrt(alpha_nu)
        self.beta_c = 1/(1-0.949*f_bnu)
        # Done setting scalar variables
        return qwarn

    # Transfer Function from Eisenstein & Hu
    def T_EH(self, kk, z):
        """Given a wavenumber in Mpc^-1, return the transfer function for the
        cosmology held in the global variables.
        Input: kk -- Wavenumber in Mpc^-1
        Output: The following are set as global variables:
	growth_cb -- the transfer function for density-weighted
        CDM + Baryon perturbations. 
 	growth_cbnu -- the transfer function for density-weighted
        CDM + Baryon + Massive Neutrino perturbations.
        The function returns growth_cb"""

        # Compute lots of constants first.
        error = self.TFmdm_set_cosm(redshift=z, omega_hdm =0.0, degen_hdm=1)

        qq = kk/self.omhh*(self.theta_cmb**2)
        # Compute the scale-dependent growth functions
        y_freestream = 17.2*self.f_hdm*(1+0.488*(self.f_hdm**(-7.0/6.0)))*((self.num_degen_hdm*qq/self.f_hdm)**2)
        temp1 = (self.growth_k0**(1.0-self.p_cb))
        temp2 = (self.growth_k0/(1+y_freestream))**0.7

        growth_cb = ((1.0+temp2)**(self.p_cb/0.7))*temp1
        growth_cbnu = (((self.f_cb**(0.7/self.p_cb))+temp2)**(self.p_cb/0.7))*temp1
        
        # Compute the master function
        gamma_eff =self.omhh*(self.alpha_gamma+(1-self.alpha_gamma)/(1+(kk*self.sound_horizon_fit*0.43)**4))
        qq_eff = qq*self.omhh/gamma_eff

        tf_sup_L = np.log(2.71828+1.84*self.beta_c*self.alpha_gamma*qq_eff)
        tf_sup_C = 14.4+325./(1+60.5*(qq_eff**1.11))
        tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*(qq_eff**2))

        qq_nu = 3.92*qq*np.sqrt(self.num_degen_hdm/self.f_hdm)
        max_fs_correction = 1+1.2*(self.f_hdm**0.64)*(self.num_degen_hdm**(0.3+0.6*self.f_hdm))/((qq_nu**-1.6)+(qq_nu**0.8))
        tf_master = tf_sup*max_fs_correction

        # Now compute the CDM+HDM+baryon transfer functions
        tf_cb = tf_master*growth_cb/self.growth_k0
        tf_cbnu = tf_master*growth_cbnu/self.growth_k0
        return tf_cb

    # Initial Power Spectrum
    def Po(self,k):
        return k**self.n

    ### Transfer Function from Bond and Efstathiou
    def T_BE(self,k,a=6.4,b=3.,c=1.7, nu=1.13):
        Gamma = self.Omega_m*self.h
        q = k/Gamma
        #print 'using BE'
        return (1+(a*q+(b*q)**1.5+(c*q)**2)**nu)**(-1./nu)

    def P(self,k,T='BE'):
        if T=='BE':
            return self.T_BE(k)**2*self.Po(k)
        else:
            return self.T_EH(k,z=0)**2*self.Po(k)
    
    def b(self,z):
        return 1./(1.+z)

    def W(self,k,M):
        R = (3*M/(4*np.pi*self.rho_m))**(1/3.)
        u = k*R
        return 3*(np.sin(u)-u*np.cos(u))/u**3

    def sigma_func(self,k,M,T='BE'):
        return k**2*self.P(k,T)*self.W(k,M)**2

    ### redshift computed here, not in transfer function.
    # OK for BE, not for EH. Need to fix this.
    # Returns normalized sigma, not sigma^2
    def Sigma(self,M, z=0,T='BE'):
        if type(M) is np.ndarray:
            results = np.array([0.0]*M.size,dtype='float')
            for i in range(0,M.size):
                results[i] = quad(self.sigma_func, 0.0, np.inf, args=(M[i],T),limit=2000)[0]
        else:
            results = quad(self.sigma_func, 0.0, np.inf, args=(M,T),limit=2000)[0]
        M8 = 4./3*np.pi*8**3*self.rho_m
        s8 = self.b(0)/(2*np.pi**2)*quad(self.sigma_func, 0.0, np.inf, args=(M8,T),limit=2000)[0]
        results *= self.sigma8**2/s8 #normalize to sigma8 at redshift 0.
        return np.sqrt(self.b(z)*results/(2*np.pi**2))

    def N(self,M,z=0,T='BE'):
        s = self.Sigma(M,z,T)
        f = np.sqrt(2/np.pi)*self.delta_c/s*np.exp(-self.delta_c**2/(2*s**2))
        return f*self.rho_m/M

    def PSNofM(self,M, z=0,T='BE'):
        """
        Use Press Schecter Formalism to compute the mass function N(m),
        fractional number density of halos in mass range (m,m+dm). Should
        integrate to total number of halos/volume.
        @param M: x-range of masses in Msun/h
        @param z: red shift at which to evaluate mass function
        @return: [x-axis in log10(mass), y-axis in log10(N(m)), xlabel, ylabel]
        """
        s = self.Sigma(M,z,T)
        tck = interpolate.splrep(M,s)
        dsdm = -interpolate.splev(M, tck, der=1)
        return self.rho_m*self.delta_c/M*np.sqrt(2/np.pi)*np.exp(-self.delta_c**2/(2*s**2))/s**2*dsdm

    def PSdNdLogM(self,M, z=0,T='BE'):
        """
        Use Press Schecter Formalism to compute the mass function
        dn(m)/dlog(M).
        @param M: x-range of masses in Msun/h
        @param z: red shift at which to evaluate mass function
        @return: [x-axis in log10(mass), y-axis in log10(N(m)), xlabel, ylabel]
        """
        s = self.Sigma(M,z,T)
        tck = interpolate.splrep(np.log10(M),s)
        dsdlogm = -interpolate.splev(np.log10(M), tck, der=1)
        return self.rho_m*self.delta_c/M*np.sqrt(2/np.pi)*np.exp(-self.delta_c**2/(2*s**2))/s**2*dsdlogm

    def STdNdM(self,M, z=0,T='BE'):
        """
        Use Sheth & Tormen Formalism to compute the mass function
        dn(m)/d(M).
        @param M: x-range of masses in Msun/h
        @param z: red shift at which to evaluate mass function
        @return: [x-axis in log10(mass), y-axis in log10(N(m)), xlabel, ylabel]
        """
        s = self.Sigma(M,z,T)
        tck = interpolate.splrep(M,s)
        dsdm = -interpolate.splev(M, tck, der=1)
        A = 0.3222
        a = 0.707
        p = 0.3
        f = A*np.sqrt(2*a/np.pi)*(1+(s**2/a/self.delta_c**2)**p)*self.delta_c/s*np.exp(-.5*a*self.delta_c**2/s**2)
        
        return f/s*self.rho_m/M*dsdm

    def STdNdLogM(self,M, z=0, T='BE'):
        """
        Use Sheth & Tormen Formalism to compute the mass function
        dn(m)/dlog(M).
        @param M: x-range of masses in Msun/h
        @param z: red shift at which to evaluate mass function
        @return: [x-axis in log10(mass), y-axis in log10(N(m)), xlabel, ylabel]
        """
        s = self.Sigma(M,z,T)
        tck = interpolate.splrep(np.log10(M),s)
        dsdlogm = -interpolate.splev(np.log10(M), tck, der=1)
        A = 0.3222
        a = 0.707
        p = 0.3
        f = A*np.sqrt(2*a/np.pi)*(1+(s**2/a/self.delta_c**2)**p)*self.delta_c/s*np.exp(-.5*a*self.delta_c**2/s**2)
        return f/s*self.rho_m/M*dsdlogm

def PlotStuff():
    MF = AnalyticMF(h=0.7, Omega_m=0.27,sigma8=0.8,n=0.95) 
    kk = 10**np.arange(-2,2,.2)
    EH = [MF.T_EH(k,z=0.0) for k in kk]
    plt.plot(kk,EH)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

def PlotStuff2():
    MF = AnalyticMF(h=0.7, Omega_m=0.27,sigma8=0.8,n=0.95) 
    kk = 10**np.arange(-2,2,.2)
    BE = [MF.T_BE(k,a=6.4,b=3.,c=1.7, nu=1.13) for k in kk]
    plt.plot(kk,BE)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
