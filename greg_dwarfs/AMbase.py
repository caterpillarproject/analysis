import numpy as np
import reionization as reion
import AnalyticMF as MF
import reionization as reion
from scipy import interpolate
import DwarfMethods as dm

ML_ratio = 1.0
#N_iter = 20000 # make 20,000 later
N_iter = 2000 # make 20,000 later

# subhalo mass function next. Make alpha positive
def Nsubs_bigger(msub, Mhost,alpha,K):
    return K * Mhost /(alpha-1) * (msub**(1-alpha) - Mhost**(1-alpha))

# generate a list of halo masses from mlow=10**7.5 to mhigh = Mhost that follows the mean of the SHMF
def generate_shm_sample_static(Mhost,alpha,K, superPoisson=False): 
    mlow = 10**7.4 # lower limit of where stars can form above 10**3 Lsun, and not dark from reionization
    if superPoisson:
        return generate_shm_sample_superPoisson_static(Mhost, alpha, K, mlow)
    # as in sawala paper. Could go to 10**7 to be safer, but the luminous fraction is almost 0 there.
    mean = Nsubs_bigger(mlow,Mhost,alpha,K)
    nsubs = np.random.poisson(mean)
    randnums = np.random.rand(nsubs) # a number between 0 and 1
    masses = (mlow**(1-alpha) - (alpha-1)/(K*Mhost) * randnums*mean)**(1/(1-alpha))
    return masses

def generate_shm_sample_superPoisson_static(Mhost, alpha, K, mlow):
    # divide into intervals
    cutoffs = np.logspace(np.log10(mlow), np.log10(Mhost), 21) 
    sI = 0.18

    final_list = []
    for i in range(len(cutoffs)-1):
        mean = K * Mhost /(alpha-1) * (cutoffs[i]**(1-alpha) - cutoffs[i+1]**(1-alpha)) # tested to be correct
        mean_rest = K * Mhost /(alpha-1) * (cutoffs[i+1]**(1-alpha) - Mhost**(1-alpha)) # tested to be correct
        variance = mean + sI**2 * (mean**2 + 2*mean*mean_rest)
        p = mean / variance
        r = mean**2 / (variance - mean)
        nsubs = np.random.negative_binomial(r,p)   # sample_neg_binomial(mean, r, p)
        randnums = np.random.rand(nsubs)
        masses = (cutoffs[i]**(1-alpha) - (alpha-1)/(K*Mhost) * randnums*mean)**(1/(1-alpha)) # tested this is correct
        final_list = np.append(final_list, masses)
    return final_list



def Nvisible_From_Ngreater(min_lum, lums):
    if type(min_lum)==list or type(min_lum)==np.ndarray:
        N_visible = []; N_std = []; per10=[]; per90=[]
        for minlum in min_lum:
            Ngreater = np.array([np.sum(lum>minlum) for lum in lums])
            N_visible.append(np.mean(Ngreater))
            N_std.append(np.std(Ngreater))
            per10.append(np.percentile(Ngreater,20))
            per90.append(np.percentile(Ngreater,80))
        return np.array(N_visible), np.array(N_std), np.array(per10), np.array(per90)
    else:
        Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
        N_visible = np.mean(Ngreater)
        N_std = np.std(Ngreater)
        per10 = np.percentile(Ngreater,20)
        per90 = np.percentile(Ngreater,80)
        return N_visible, N_std, per10, per90

class AMbase(object):
    """
    methods that are standard to all AM models.
    methods that are unique should be over written by particular AM models
    """
    def __init__(self):
	return

    # given halo mass, what is stellar mass?
    def halo_to_stellar_mass(self,mvir,a=1):
	if type(mvir)==list:
            mvir = np.array(mvir)
        return self.getStellarMass(mvir,a)

    # given stellar mass/luminosity what is DM halo mass?
    # this assumes a scatter of 0
    def stellar_to_halo_mass(self, mstar, a=1):
	#USE INTERPOLATION     
 	from scipy import interpolate
 	mass_array = 10**np.arange(7,12,.07)
	mstar_array = self.halo_to_stellar_mass(mass_array,a)
        tck = interpolate.splrep(np.log10(mstar_array),np.log10(mass_array), k=1)
        return 10**interpolate.splev(np.log10(mstar),tck)


    def generate_shm_sample(self, Mhost, factor=1): 
        mlow = 10**7.4 # lower limit of where stars can form above 10**3 Lsun, and not dark from reionization
        if self.label[0:4] == 'Behr':
            mlow = 10**6.8
        if self.superPoissonian:
            return self.generate_shm_sample_superPoisson(Mhost, mlow, factor)
        # as in sawala paper. Could go to 10**7 to be safer, but the luminous fraction is almost 0 there.
        mean = Nsubs_bigger(mlow,Mhost,self.alpha,self.K)
        nsubs = np.random.poisson(mean*factor)  # only change the number of subs I want to generate with factor
        randnums = np.random.rand(nsubs) # a number between 0 and 1
        masses = (mlow**(1-self.alpha) - (self.alpha-1)/(self.K*Mhost) * randnums*mean)**(1/(1-self.alpha))
        return masses

    def generate_shm_sample_superPoisson(self, Mhost, mlow, factor=1):
        # divide into intervals
        cutoffs = np.logspace(np.log10(mlow), np.log10(Mhost), 21) 
        sI = 0.18 #0.14

        final_list = []
        for i in range(len(cutoffs)-1):
            mean = self.K * Mhost /(self.alpha-1) * (cutoffs[i]**(1-self.alpha) - cutoffs[i+1]**(1-self.alpha)) # tested to be correct
            mean_rest = self.K * Mhost /(self.alpha-1) * (cutoffs[i+1]**(1-self.alpha) - Mhost**(1-self.alpha)) # tested to be correct
            
            mean_mod=factor*mean
            mean_rest_mod=factor*mean_rest

            variance = mean_mod + sI**2 * (mean_mod**2 + 2*mean_mod*mean_rest_mod)
            p = mean_mod / variance  # leave this unaffected by multiplicative factor
            r = mean_mod**2 / (variance - mean_mod)
            nsubs = np.random.negative_binomial(r,p)   # sample_neg_binomial(mean, r, p)  # this is the only one I want affected
            randnums = np.random.rand(nsubs)
            masses = (cutoffs[i]**(1-self.alpha) - (self.alpha-1)/(self.K*Mhost) * randnums*mean)**(1/(1-self.alpha)) # tested this is correct
            final_list = np.append(final_list, masses)
        return final_list


    # assume to be at z=0
    # input should be Mhalo, the mass fitting the definition from the AM model.
    def mass_to_rvir(self, Mhalo):
        # first convert Mhalo to Mvir
        if type(Mhalo)==list or type(Mhalo)==np.ndarray:
            rvir = []
            for m in Mhalo:
                mvir, r = dm.convert_Mhalo_d_d(m, self.delta_crit, 103.86)
                rvir.append(r)
            return np.array(rvir)
        else:
            mvir, rvir = dm.convert_Mhalo_d_d(Mhalo, self.delta_crit, 103.86)
        return rvir
        #G = 4.157e-39  # kpc^3/Msun/s^2
        #H = dm.Hubble(a=1)
        #rho_crit = (3*H**2)/(8*np.pi*G)
        #rvir = (Mvir*3/(4*np.pi*self.delta_crit*rho_crit))**(1/3.)

    def mstar_to_mhalo(self,masses,a=1):
        if type(masses)==list or type(masses)==np.ndarray:
            return np.array([self.mstar_to_mhalo_single(Mstar,a) for Mstar in masses])
        else:
            return self.mstar_to_mhalo_single(masses,a)



    # use for the host halos. 
    def mstar_to_mhalo_single(self, Mstar, a=1):
        def gauss(masses, sigma, mean):
            mstar_dist = np.log10(mean) - np.log10(self.getStellarMass(masses,a))
            return 1./np.sqrt(2*np.pi * sigma**2)  * np.e**(-( mstar_dist**2 )/(2*sigma**2))
    
        Mhalo = self.stellar_to_halo_mass(Mstar, a=1)
        halo_masses = np.logspace(np.log10(Mhalo)-.8,  np.log10(Mhalo)+0.3, 30)
        sigmas = self.get_sigma(halo_masses)
        
        #sigmas = np.array([0.5]*len(halo_masses)) # only for experiment purposes

        #print sigmas
        weight1 = gauss(halo_masses, sigmas, Mstar)
         
        h0 = 0.6711; Om = 0.31; sig8 = 0.83; ns = 0.96
        
        Theoretical = MF.AnalyticMF(h=h0, Omega_m=Om,sigma8=sig8, n=ns)
        #STdndlogM = Theoretical.STdNdLogM(Mhalo, z=1.)
        weight2 = Theoretical.STdNdLogM(halo_masses,z=1)

        lum_tck, _, _ = reion.get_lum_fraction_mpeak()
        # shift lum frac by 0.1 dex
        shifted_masses = 10**(np.log10(halo_masses)-.16)
        lum_chances = reion.frac_luminous(shifted_masses, self.lum_tck) # double check this is between 0 and 1
        #print lum_chances, 'lum chances'
        weight3 = lum_chances
        #print weight3

        values = weight1*weight2*weight3
        values = values / np.sum(values)
        cum = np.cumsum(values)
        #print np.log10(cum)
    
        tck = interpolate.splrep(np.log10(cum)[3:-3] , np.log10(halo_masses)[3:-3], k=1)  
        mass = 10**interpolate.splev(np.log10(0.5), tck)  # THIS IS THE FINAL ANSWER
       

        # ESTIMATE typical uncertainty in total mass
        #m1siglow = 10**interpolate.splev(np.log10(0.1587), tck)
        #m1sighigh = 10**interpolate.splev(np.log10(0.8413), tck)
        #print m1siglow, mass, m1sighigh, 'mlow, median, hight'      
        #print (- mass + m1sighigh)/mass, '% change up'
        #print (mass - m1siglow)/mass, '% change down'

        #print mass/Mhalo, Mstar, 'ratio of new mass to 0 scatter mass'


        if mass/Mhalo < 0.5 or mass/Mhalo > 1.0:
            raise Exception('mstar to mhalo failure!!!!')
        if mass is None:
            raise Exception('mstar to mhalo failure!!!!')
        return mass


    def P_at_least_one(self, halo_masses, min_lum):
        if self.isPoisson:
            lowest_mass = self.stellar_to_halo_mass(ML_ratio*min_lum)
            mean = Nsubs_bigger(lowest_mass,halo_masses,self.alpha,self.K)
            Pgt1 = 1 - np.e**-mean
            return Pgt1
        else:
            frac_gte_one = []
            for halo_mass in halo_masses:
                samples = self.generate_stellar_mass_samples(halo_mass,N=N_iter) # generate list of lists of stellar masses
                lums = samples/ML_ratio
                Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
                frac_gte_one.append(np.sum(Ngreater>0) / float(len(Ngreater)))
            Pgt1 = np.array(frac_gte_one)
            return Pgt1


    def ngreater(self,halo_masses,min_lum, percentile=False):
        if self.isPoisson:
            lowest_mass = self.stellar_to_halo_mass(ML_ratio*min_lum)
            N_visible = Nsubs_bigger(lowest_mass, halo_masses,self.alpha,self.K)
            if percentile:
                print 'percentile not supported with poisson mode on'
                return None
            else:
                return N_visible, np.sqrt(N_visible)
        else:
            if type(halo_masses)==list or type(halo_masses)==np.ndarray:
                num_gte_min = []; spread = []; p10=[]; p90=[]
                for halo_mass in halo_masses:
                    samples = self.generate_stellar_mass_samples(halo_mass,N=N_iter) # generate list of lists of stellar masses
                    lums = samples/ML_ratio
                    Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
                    num_gte_min.append( np.mean(Ngreater))
                    spread.append( np.std(Ngreater))
                
                    if percentile:
                        # get 10 and 90th percentiles
                        p10.append(np.percentile(Ngreater,20))
                        p90.append(np.percentile(Ngreater,80))

                N_visible = np.array(num_gte_min); N_std = np.array(spread);
                if percentile:
                    per10 = np.array(p10); per90 = np.array(p90)
            else:
                samples = self.generate_stellar_mass_samples(halo_masses,N=N_iter) # generate list of lists of stellar masses
                lums = samples/ML_ratio
                N_visible, N_std, per10, per90 = Nvisible_From_Ngreater(min_lum, lums)
                #Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
                #N_visible = np.mean(Ngreater)
                #N_std = np.std(Ngreater)
                #if percentile:
                #    per10 = np.percentile(Ngreater,20)
                #    per90 = np.percentile(Ngreater,80)
            if not percentile:
                return N_visible, N_std
            else:
                return N_visible, N_std, per10, per90


    # generate list of lists of stellar masses
    # must have getStellarMass_random() implemented in subclass
    # multiplicative value applied to the mean as a way to handle subvolumes
    def generate_stellar_mass_samples(self,halo_mass,N=N_iter, factor=1):
        stellar_mass_lists=[]
    	for i in range(N):
            dm_masses = self.generate_shm_sample(halo_mass, factor)
            if self.reionization:
                lum_chances = reion.frac_luminous(dm_masses, self.lum_tck) # hopefully a list of fractions between 0 and 1
                randnums = np.random.rand(len(lum_chances))
                mask = randnums < lum_chances # chooses which halos will form stars
                dm_masses = dm_masses[mask] # dm masses that host stars
            stellar_mass_lists.append(self.getStellarMass_random(dm_masses,a=self.scale_factor))   # this is explicitly for subhalos
   	return np.array(stellar_mass_lists)



    def get_field_total(self,dm_masses,min_lum):
        N, _ = self.ngreater(dm_masses,min_lum)
        return np.sum(N)

    def get_field_distr(self,dm_masses,min_lum,N=N_iter):
        distr=np.array([0]*N)
        if type(dm_masses)==list or type(dm_masses)==np.ndarray:
            for halo_mass in dm_masses:
                samples = self.generate_stellar_mass_samples(halo_mass,N)
                lums = samples/ML_ratio
                Ngreater = np.array([np.sum(lum>min_lum) for lum in lums]) #should be 10,000 long
                distr+=Ngreater
            return distr
        else:
            samples = self.generate_stellar_mass_samples(dm_masses,N)
            lums = samples/ML_ratio
            Ngreater = np.array([np.sum(lum>min_lum) for lum in lums]) #should be 10,000 long
            return Ngreater







