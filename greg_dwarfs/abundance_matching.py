import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PlotParams import *
from AMbase import *
import infall_times
import random
#import DwarfMethods as dm
import reionization as reion
from matplotlib import colors
from scipy import interpolate

##### HOST SHMF DEFAULT VALUES MUST BE SET TO THE CORRECT VALUES AFTER THEY ARE COMPUTED!!!!!!######

ML_ratio = 1
N_iter = 1000  # 10000 for final plots, but for speed of testing things, make 2000


#### SHFM VALUES ######
m200_alpha_host = 1.884924; m200_k_host = .000989
m200inf_alpha_host = 1.837915; m200inf_k_host = .000854
mpeak_alpha_host = 1.870010; mpeak_k_host = .001878
m350peak_alpha_host = 1.867326; m350peak_k_host = .002003
mvir_alpha_host = 1.874085; mvir_k_host = .000708


m200_alpha_field = 1.931428; m200_k_field = .002578
m200inf_alpha_field = 1.810661 ; m200inf_k_field = .000635
mpeak_alpha_field = 1.823999 ; mpeak_k_field = .000892
m350peak_alpha_field = 1.811803; m350peak_k_field = .000765


#mvir_alpha_lmc =1.97255832244; mvir_k_lmc = 0.00505623635311
mvir_alpha_lmc =1.90; mvir_k_lmc = 0.00122202932021
mpeak_alpha_lmc = 1.925666; mpeak_k_lmc = .005879
m200inf_alpha_lmc = 1.880723; m200inf_k_lmc = .002194
m350peak_alpha_lmc = 1.881158; m350peak_k_lmc = .002820


########################




# assume a slope of -1.90 for the host at z=0. 
# r_ratio is r/r_vir, where r_vir depends on the host mass Mvir.
def radial_dependence(r_ratio):
    return -.00015+.00071*r_ratio + .00086*r_ratio**2 - .00031*r_ratio**3



# Not updated with any of the LMC stuff
class Sawala(AMbase):
    def __init__(self, hostSHMF=False, hostAlpha= m200_alpha_host, hostK= m200_k_host): # Host SHMF from the z=0 radial fit.
	self.slope1 = 2.2#2.00
	self.slope2 = 0.41#.546
	self.M1 = 10**10.347 / (10**(7.6143/self.slope1))  #34.7e5
	self.M2 = 10**8.0214 / (10**(4.9/self.slope2)) #0.11144

        # for Mpeak instead of z=0
	self.slope1 = 2.1#2.00
	self.M1 = 10**10.347 / (10**(7.6143/self.slope1))  #34.7e5
        self.slope2 = 0.10 #.546
	self.M2 = 10**8.0214 / (10**(4.98/self.slope2)) #0.11144

        # these values correspond to equation 1 in the Sawala abundance of not just dark matter haloes paper
	self.a=0.69
	self.Mt = 10**11.6
	self.w = 0.79
        self.b = 0.98
        self.ls = '-'

        # SHMF parameters. Slope alpha and constant, K
        self.alpha = m200_alpha_field 
        self.K = m200_k_field
        if hostSHMF:
            self.alpha = hostAlpha
            self.K = hostK
        self.superPoissonian = False

        self.color = 'cyan'
        self.label = 'Sawala'

        self.sigma = 0.5  ## the lognormal scatter in the SHMF relationship. Estimated from plot
        self.isPoisson=False

        self.delta_crit = 200
        self.scale_factor = 1.0

    # Infered from figure 4 lower panel of Bent by Baryons
    # I should add in some scatter too
    # high mass end needs to be adjusted - force it to converge to moster for instance?????
    # yes, at 10**11, make it switch over to Moster et al.
    def getStellarMass_single(self,M,a=1):
        if M < 10**11:
            return (M/(self.M2))**self.slope2 + (M/(self.M1))**self.slope1
        else:
            model=Moster()
            return model.getStellarMass(M,a)

    def getStellarMass(self,M,a=1):
        if type(M)==list or type(M)==np.ndarray:
            stellar=[-1]*len(M)
            for i,mass in enumerate(M):
                stellar[i]=self.getStellarMass_single(mass)
            return np.array(stellar)
        else:
            return self.getStellarMass_single(M)

    def getStellarMass_up1sig(self,M,a=1):
        return 10**( np.log10(self.getStellarMass(M)) + self.sigma)
    
    def getStellarMass_down1sig(self,M,a=1):
        return 10**( np.log10(self.getStellarMass(M)) - self.sigma)

    # this is for subhalos, not centrals
    # from equation 2 in "The abundance of not just dark matter haloes"
    def reduceDMmassEQ2(self,M):
	return M* (.65+(M/10**11.4)**.51)/(1+ (M/10**11.4)**.51)

    # from equation 1 in "The abundance of not just dark matter haloes"
    def reduceDMmass(self,M):
        return M* self.b * (self.a/self.b + (M/self.Mt)**self.w ) / (1+ (M/self.Mt)**self.w)

    # Taken from The Chosen Few Figure 3
    # given halo mass, what fraction are luminous?
    def frac_luminous(self,M,a=1):
        from scipy import interpolate
        frac_array = np.array([0, 0.01389365137212395, 0.03227229587377156,  0.0402324326871204, 0.09594819785818098,   0.17944972724833808,0.3697539417495146, 0.7910273213219392, 0.8458747813611525, 1  ])
        log_mass = np.array([7.283067472982055, 7.998634558928878,8.247742935525217,  8.49887928058564, 8.747863613591717, 8.996755635088757,9.247286185871147, 9.497049397234692, 9.746036614975427, 10])
        tck = interpolate.splrep(log_mass,frac_array, k=1)
        if len(M)==0:
            return np.array([])
        return interpolate.splev(np.log10(M),tck)

    # should only input reduced mass into this function
    def getStellarMass_random(self,M,a=1):
        # get random scatter, perturbation
        if type(M)==list or type(M)==np.ndarray:
            perturbation = np.array([random.gauss(0,self.sigma) for i in range(len(M))])
        else:
            perturbation = random.gauss(0,self.sigma)
        return 10**( np.log10(((M/(self.M2))**self.slope2 + (M/(self.M1))**self.slope1)) +perturbation )

    # !!! This function includes all of the new reionization and DM mass suppression stuff!! 
    def generate_stellar_mass_samples(self,halo_mass,N=N_iter,factor=1):
        stellar_mass_lists=[]
    	for i in range(N):
            # generate dm_masses
            dm_masses = generate_shm_sample(self, halo_mass,factor)
            # reduce dm_masses by certain amount
            dm_masses = self.reduceDMmass(dm_masses)
            lum_chances = self.frac_luminous(dm_masses) # hopefully a list of fractions between 0 and 1
            randnums = np.random.rand(len(lum_chances))
            mask = randnums < lum_chances # chooses which halos will form stars
            dm_masses = dm_masses[mask] # dm masses that host stars

            stellar_mass_lists.append(self.getStellarMass_random(dm_masses,a=1.0))
   	return np.array(stellar_mass_lists)


    # must implement this differently
    def P_at_least_one(self, halo_masses, min_lum):
        frac_gte_one = []
        for halo_mass in halo_masses:
            samples = self.generate_stellar_mass_samples(halo_mass,N=N_iter) # generate list of lists of stellar masses
            lums = samples/ML_ratio
            Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
            frac_gte_one.append(np.sum(Ngreater>0) / float(len(Ngreater)))
        Pgt1 = np.array(frac_gte_one)
        return Pgt1


    def ngreater(self,halo_masses,min_lum):
        if type(halo_masses)==list or type(halo_masses)==np.ndarray: 
            num_gte_min = []
            spread = []
            for halo_mass in halo_masses:
                samples = self.generate_stellar_mass_samples(halo_mass,N=N_iter) # generate list of lists of stellar masses
                lums = samples/ML_ratio
                Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
                num_gte_min.append( np.mean(Ngreater) ) # this is taking the average
                spread.append( np.std(Ngreater))
            N_visible = np.array(num_gte_min); N_std = np.array(spread)
        else:
            samples = self.generate_stellar_mass_samples(halo_masses,N=N_iter) # generate list of lists of stellar masses
            lums = samples/ML_ratio
            Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
            N_visible = np.mean(Ngreater)
            N_std = np.std(Ngreater)
        return N_visible, N_std

    def get_sigma(self, M=None):
        if M is None:
            return self.sigma
        if type(M)==list or type(M)==np.ndarray:
            return np.array([self.sigma]*len(M))
        else:
            return self.sigma



"""        
# data points for interpolating
10, 1
9.746036614975427, 0.8458747813611525
9.497049397234692, 0.7910273213219392
9.247286185871147, 0.3697539417495146
8.996755635088757, 0.17944972724833808
8.747863613591717, 0.09594819785818098
8.49887928058564, 0.0402324326871204
8.247742935525217, 0.03227229587377156
7.998634558928878, 0.01389365137212395
7.283067472982055, -0.00041347863421270503
"""
	

class Brook(AMbase):
    def __init__(self,hostSHMF=False,hostAlpha=m350peak_alpha_host, hostK=m350peak_k_host, reionization=False,z=13, lmcSHMF=False,lmcAlpha=m350peak_alpha_lmc, lmcK=m350peak_k_lmc, catreion=False, vmax_ach=9.48535156 , vmax_filt= 23.54248047): # host SHMF should match plot in the hostfigs folder, and value written in paper
        super(Brook,self).__init__()

        self.alphaAM = 3.1
        self.M0 = 79.6  # 63.1 #63.1   #79.6    # for the high mass transition it no longer works too well.... MUST FIX!!

        # SHMF parameters. Slope alpha and constant, K of the field values
        self.alpha = m350peak_alpha_field
        self.K = m350peak_k_field
        self.color = colors.hex2color('#e41a1c')     #'red'
        self.ls = '--'
        self.label = 'Brook'
        self.isPoisson=True
        self.z=z
        self.delta_crit = 350
        self.scale_factor = 1.0
        self.moster_model=Moster(no_load_reionization = True)
        self.superPoissonian = False
        self.catreion = catreion

        if hostSHMF:
            self.alpha = hostAlpha
            self.K = hostK
            self.superPoissonian = True

        if lmcSHMF:
            self.alpha = lmcAlpha
            self.K = lmcK
            self.superPoissonian = True
            self.isPoisson=False

        """
        self.reionization = reionization
        if reionization:
            #self.color = 'midnightblue'
            self.label = 'Brook + Reion'
            #self.lum_tck, _, _ = reion.get_lum_fraction_Barber()
            #self.lum_tck, _, _,_,_,_ = reion.get_lum_fraction_caterpillar(mdef='m350',z=self.z)
            self.lum_tck, _, _ = reion.get_lum_fraction_m350()

            self.isPoisson=False
            self.ls = '-'
        """

        self.reionization = reionization
        if reionization:
            self.label = 'Brook + Reion'
            #self.lum_tck, _, _ = reion.get_lum_fraction_Barber()
            if catreion:
                self.lum_tck, _, _ = reion.get_lum_fraction_caterpillar(mdef='m350',z=self.z,vmax_ach=vmax_ach, vmax_filt=vmax_filt)
            else:
                self.lum_tck, _, _ = reion.get_lum_fraction_m350()

            self.isPoisson=False
            self.ls = '-'
        else:
            self.lum_tck, _, _ = reion.get_lum_fraction_m350()






    def getStellarMass_single(self,M,a=1):
        Mcrit = self.M0*10**6 * (10**8)**(1/3.1)   # 10.48 for M0=79.6, 10.38 for M0 = 63.1
        Msmooth =  164850728228.91016   # for 3*10**9
        #Msmooth = model.stellar_to_halo_mass(3*10**9)  # this is very slow. do it once, then hard code the value
        #print Msmooth, 'moster mass at 10**9'

        if M <= Mcrit:   # how did we get to 10.48?
            return (M/(self.M0*10**6))**self.alphaAM
        elif M > Mcrit and M < Msmooth:
            # linearly interpolate between brook at 10**8 and moster at 10**9
            return 10**( 8 + ( np.log10(3.e9)-8.0)/( np.log10(Msmooth) - np.log10(Mcrit)) *(np.log10(M) - np.log10(Mcrit)))

        else:   # M < Msmooth
            frac = 1 #((Mcrit/(self.M0*10**6))**self.alphaAM) / model.getStellarMass(Mcrit,a=1)
            #print frac, Mcrit, 'frac and Mcrit in Brook'
            return frac*self.moster_model.getStellarMass(M,a)

    def getStellarMass(self,M,a=1):
        if type(M)==list or type(M)==np.ndarray:
            stellar=[-1]*len(M)
            for i,mass in enumerate(M):
                stellar[i]=self.getStellarMass_single(mass)
            return np.array(stellar)
        else:
            return self.getStellarMass_single(M)

    def getStellarMass_random(self,M,a=1):
        return self.getStellarMass(M,a)
                
    def get_sigma(self, M=None):
        if M is None:
            return 0.2
        if type(M)==list or type(M)==np.ndarray:
            return np.array([0.2]*len(M))
        else:
            return 0.2

class GarrisonKimmel(AMbase):
    def __init__(self,hostSHMF=False,hostAlpha=mpeak_alpha_host, hostK= mpeak_k_host, reionization=False,z=13, lmcSHMF=False,lmcAlpha=mpeak_alpha_lmc, lmcK= mpeak_k_lmc, catreion=False, vmax_ach=9.48535156 , vmax_filt= 23.54248047): # host SHMF should match plot in the hostfigs folder, and value written in paper
        self.epsilon0=-1.777
        self.M10 = 11.514
        self.alphaAM0=-1.92
        self.delta0=3.508
        self.gamma0=0.316
        self.ls = '--'
        self.z=z
        self.delta_crit = 103.86
        # SHMF parameters. Slope alpha and constant, K
        #self.alpha = 1.867  # these values at z=0 mass
        #self.K = .000703
        self.alpha = mpeak_alpha_field    # these values for Mpeak, bryan & norman
        self.K =  mpeak_k_field  
        self.scale_factor = 1.0
        self.superPoissonian = False
        self.catreion = catreion

        if hostSHMF:
            self.alpha= hostAlpha
            self.K = hostK
            self.superPoissonian = True

        if lmcSHMF:
            self.alpha = lmcAlpha
            self.K = lmcK
            self.superPoissonian = True
            self.isPoisson=False

        self.color = colors.hex2color('#984ea3')   #'cyan'
        self.label = 'GK14'
        self.isPoisson=True


        self.reionization = reionization
        if reionization:
            self.label = 'GK14 + Reion'
            #self.lum_tck, _, _ = reion.get_lum_fraction_Barber()
            if catreion:
                self.lum_tck, _, _ = reion.get_lum_fraction_caterpillar(mdef='mpeak',z=self.z,vmax_ach=vmax_ach, vmax_filt=vmax_filt)
            else:
                self.lum_tck, _, _ = reion.get_lum_fraction_mpeak()

            self.isPoisson=False
            self.ls = '-'
        else:
            self.lum_tck, _, _ = reion.get_lum_fraction_mpeak()

    def f(self,x):
        return -np.log10(10**(x*self.alphaAM())+1)+ self.delta()*(np.log10(1+np.e**x))**self.gamma() / (1+np.e**(10**-x))

    def getStellarMass(self,M,a=1):
        return 10**(self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0))

    def getStellarMass_random(self,M,a=1):
        return self.getStellarMass(M,a)

    def M1(self):
        return self.M10
    def epsilon(self):
        return self.epsilon0
    def alphaAM(self):
        return self.alphaAM0
    def delta(self):
        return self.delta0
    def gamma(self):
        return self.gamma0

    def get_sigma(self, M=None):
        if M is None:
            return 0.2
        if type(M)==list or type(M)==np.ndarray:
            return np.array([0.2]*len(M))
        else:
            return 0.2



class GarrisonKimmel16(AMbase):
    # default values for hostK and hostAlpha are not corrrect for peak values!!!
    def __init__(self,sigma=0.8,panel='satellites',variable=False,hostSHMF=False,hostAlpha= mpeak_alpha_host, hostK= mpeak_k_host ,reionization=False,z=13, lmcSHMF=False,lmcAlpha= mpeak_alpha_lmc, lmcK= mpeak_k_lmc):
        if panel=='satellites':
            self.alphaAM0 = -(0.14* sigma**2 + .14*sigma+1.79)
        if panel=='field':
            self.alphaAM0 = -(0.24*sigma**2+0.16*sigma+1.99)        

        self.epsilon0=-1.777
        self.M10 = 11.514
        self.delta0=3.508
        self.gamma0=0.316
        self.sigma = sigma
        self.ls = '--'
        self.z=z
        self.delta_crit = 103.86
        # SHMF parameters. Slope alpha and constant, K
        self.alpha = mpeak_alpha_field
        self.K = mpeak_k_field
        self.scale_factor = 1.0
        self.superPoissonian = False

        if hostSHMF:
            self.alpha= hostAlpha
            self.K = hostK
            self.superPoissonian = True

        if lmcSHMF:
            self.alpha = lmcAlpha
            self.K = lmcK
            self.superPoissonian = True
            self.isPoisson=False

        self.color = 'orange'
        self.label = 'GK16_const'
        self.isPoisson=False

        self.reionization = reionization
        if reionization:
            #self.color = 'midnightblue'
            self.label = 'GK16_const + Reion'
            #self.lum_tck, _, _ = reion.get_lum_fraction_Barber()
            #self.lum_tck, _, _ = reion.get_lum_fraction_caterpillar(mdef='mpeak',z=self.z)
            self.lum_tck, _, _ = reion.get_lum_fraction_mpeak()
            self.isPoisson=False
            self.ls = '-'

    def f(self,x):
        return -np.log10(10**(x*self.alphaAM())+1)+ self.delta()*(np.log10(1+np.e**x))**self.gamma() / (1+np.e**(10**-x))

    def getStellarMass(self,M,a=1):
        return 10**(self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0))

    def getStellarMass_up1sig(self,M,a=1):
        return 10**( (self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0)) +self.sigma)
    
    def getStellarMass_down1sig(self,M,a=1):
        return 10**( (self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0)) -self.sigma)

    def getStellarMass_random(self,M,a=1):
        if type(M)==list or type(M)==np.ndarray:
            perturbation = np.array([random.gauss(0,self.sigma) for i in range(len(M))])
        else:
            perturbation = random.gauss(0,self.sigma)
        return 10**( (self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0)) +perturbation )

    def M1(self):
        return self.M10
    def epsilon(self):
        return self.epsilon0
    def alphaAM(self):
        return self.alphaAM0
    def delta(self):
        return self.delta0
    def gamma(self):
        return self.gamma0


    def get_sigma(self, M=None):
        if M is None:
            return self.sigma
        if type(M)==list or type(M)==np.ndarray:
            return np.array([self.sigma]*len(M))
        else:
            return self.sigma



class GK16_grow(AMbase):
    def __init__(self,plotgamma=-0.2,panel='satellites',variable=False,hostSHMF=False,hostAlpha=mpeak_alpha_host, hostK=mpeak_k_host, reionization=False,z=13, lmcSHMF=False,lmcAlpha=mpeak_alpha_lmc, lmcK=mpeak_k_lmc, catreion=False, vmax_ach=9.48535156 , vmax_filt= 23.54248047, bent = False, bentsigdex = 0.4, color= colors.hex2color('#ff7f00'), label = None):
        # plotgamma is from x axis of figure 4 in the paper
        if panel=='satellites':
            self.alphaAM0 = -(0.25*plotgamma**2 - 1.37*plotgamma + 1.69)
            #print self.alphaAM0, 'should be close to 1.8'
        if panel=='field':
            self.alphaAM0 = -(0.47*plotgamma**2 - 1.48*plotgamma + 1.81)

        self.plotgamma = plotgamma
        self.epsilon0=-1.777
        self.M10 = 11.514
        self.delta0=3.508
        self.gamma0=0.316
        self.ls = '--'
        self.z=z
        self.delta_crit = 103.86
        #self.sigma = sigma
        
        # SHMF parameters. Slope alpha and constant, K
        self.alpha = mpeak_alpha_field
        self.K = mpeak_k_field
        self.color = color # colors.hex2color('#ff7f00')   #'purple'
        self.label = 'GK16'
        self.isPoisson=False
        self.scale_factor = 1.0
        self.superPoissonian = False
        self.catreion = catreion
        self.bent = bent
        self.bentsigdex = bentsigdex

        if hostSHMF:
            self.alpha= hostAlpha
            self.K = hostK
            self.superPoissonian = True

        if lmcSHMF:
            self.alpha = lmcAlpha
            self.K = lmcK
            self.superPoissonian = True
            self.isPoisson=False

        self.reionization = reionization
        if reionization:
            #self.color = 'midnightblue'
            self.label = 'GK16 + Reion'
            #self.lum_tck, _, _ = reion.get_lum_fraction_Barber()
            #self.lum_tck, _, _ = reion.get_lum_fraction_caterpillar(mdef='mpeak',z=self.z)
            if catreion:
                self.lum_tck, _, _ = reion.get_lum_fraction_caterpillar(mdef='mpeak',z=self.z,vmax_ach=vmax_ach, vmax_filt=vmax_filt)
            else:
                self.lum_tck, _, _ = reion.get_lum_fraction_mpeak()
            self.isPoisson = False
            self.ls = '-'
        else:
            self.lum_tck, _, _ = reion.get_lum_fraction_mpeak()

        if label != None:
            self.label = label

    def get_sigma(self,M):
        if self.bent:
            mstarcrit = 1.5*10**3
            Mcrit = self.find_mcrit(mstarcrit)
            Mcrit = 10**11  # make it constant variation

        # above Mpeak, sigma = 0.2, then it grows linearly.
        # plotgamma negative, so this should really grow!!!
        if type(M)==list or type(M)==np.ndarray:
            sigs = np.array([0.2]*len(M))                                    
            sigs[M<10**11]+=self.plotgamma * (np.log10(M[M<10**11]) - self.M1())
            if self.bent:
                sigs[M < Mcrit] = self.bentsigdex
        else:
            if self.bent:
                if M < Mcrit:
                    return self.bentsigdex
            if M>10**11:
                return 0.2
            else:
                return self.plotgamma * (np.log10(M) - self.M1())
        return np.array(sigs)


    def f(self,x):
        return -np.log10(10**(x*self.alphaAM())+1)+ self.delta()*(np.log10(1+np.e**x))**self.gamma() / (1+np.e**(10**-x))


    def find_mcrit(self, mstar, a=1):
   	mass_array = 10**np.arange(7,10,.07)
	mstar_array = 10**(self.epsilon()+self.M1() + self.f(np.log10(mass_array) -self.M1()) - self.f(0))
        tck = interpolate.splrep(np.log10(mstar_array),np.log10(mass_array), k=1)
        return 10**interpolate.splev(np.log10(mstar),tck)

    def getStellarMass_single(self,M,a=1):
        mstarcrit = 1.5*10**3
        Mcrit = self.find_mcrit(mstarcrit)  #  what mass is mstar = mstarcrit

        if M <= Mcrit:  
            const = mstarcrit / (Mcrit**0.2)
            return const * M ** 0.2
        else:
            return 10**(self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0))

    def getStellarMass_bent(self,M,a=1):
        if type(M)==list or type(M)==np.ndarray:
            stellar=[-1]*len(M)
            for i,mass in enumerate(M):
                stellar[i]=self.getStellarMass_single(mass)
            return np.array(stellar)
        else:
            return self.getStellarMass_single(M)
   


    def getStellarMass(self,M,a=1):
        if self.bent:
            return self.getStellarMass_bent(M,a)
        else:
            return 10**(self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0))

    def getStellarMass_up1sig(self,M,a=1):
        return 10**(np.log10(self.getStellarMass(M,a))+self.get_sigma(M))
        #return 10**( (self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0)) +self.get_sigma(M))
    
    def getStellarMass_down1sig(self,M,a=1):
        return 10**(np.log10(self.getStellarMass(M,a))-self.get_sigma(M))
        #return 10**( (self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0)) -self.get_sigma(M))

    def getStellarMass_random(self,M,a=1):
        if type(M)==list or type(M)==np.ndarray:
            perturbation = np.array([random.gauss(0,self.get_sigma(mpeak)) for i,mpeak in zip(range(len(M)), M) ])
        else:
            perturbation = random.gauss(0,self.get_sigma(M))
        return 10**(np.log10(self.getStellarMass(M,a))+ perturbation)
        #return 10**( (self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0)) +perturbation )

    def M1(self):
        return self.M10
    def epsilon(self):
        return self.epsilon0
    def alphaAM(self):
        return self.alphaAM0
    def delta(self):
        return self.delta0
    def gamma(self): 
        return self.gamma0

# not updated with the LMC stuff
class Behroozi(AMbase):
    def __init__(self,hostSHMF=False,hostAlpha=mpeak_alpha_host, hostK=mpeak_k_host,reionization=False,z=13, lmcSHMF=False,lmcAlpha=mpeak_alpha_lmc, lmcK=mpeak_k_lmc, scale_factor=0.53):
        self.epsilon0=-1.777
        self.epsilona=-0.006
        self.epsilonz=0.00
        self.epsilona2=-0.119
        self.M10 = 11.514
        self.M1a=-1.793
        self.M1z=-0.251
        self.alphaAM0=-1.412
        self.alphaAMa=0.731
        self.delta0=3.508
        self.deltaa=2.608
        self.deltaz=-0.043
        self.gamma0=0.316
        self.gammaa=1.319
        self.gammaz=0.279
 
        # SHMF parameters. Slope alpha and constant, K
        self.alpha = mpeak_alpha_field
        self.K = mpeak_k_field
        self.color =  colors.hex2color('#4daf4a')  #'orange'
        self.label = 'Behroozi'
        self.isPoisson=True
        self.ls = '--'
        self.z=z
        self.delta_crit = 103.86
        self.scale_factor = scale_factor  # 0.53  # z=0.9 is approx t infall = 7.5  gyr ago so a = .53
        self.superPoissonian = False

        if hostSHMF:
            self.alpha= hostAlpha
            self.K = hostK
            self.superPoissonian = True

        if lmcSHMF:
            self.alpha = lmcAlpha
            self.K = lmcK
            self.superPoissonian = True
            self.isPoisson=False

        self.reionization = reionization
        if reionization:
            #self.color = 'midnightblue'
            self.label = 'Behroozi + Reion'
            #self.lum_tck, _, _ = reion.get_lum_fraction_Barber()
            #self.lum_tck, _, _ = reion.get_lum_fraction_caterpillar(mdef='mpeak',z=self.z)
            self.lum_tck, _, _ = reion.get_lum_fraction_mpeak()
            self.isPoisson=False
            self.ls = '-'
        else:
            self.lum_tck, _, _ = reion.get_lum_fraction_mpeak()
 
 
    def f(self,x,a):
        return -np.log10(10**(x*self.alphaAM(a))+1)+ self.delta(a)*(np.log10(1+np.e**x))**self.gamma(a) / (1+np.e**(10**-x))
 
    def getStellarMass(self,M,a=1):
        return 10**(self.epsilon(a)+self.M1(a) + self.f(np.log10(M) -self.M1(a),a) - self.f(0,a))
 
    def getStellarMass_random(self,M,a=1):
        return self.getStellarMass(M,a)
 
    def nu(self,a):
        return np.e**(-4*a**2)
    def M1(self,a):
        z=1/float(a)-1
        return self.M10 + (self.M1a*(a-1)+self.M1z*z) * self.nu(a)
    def epsilon(self,a):
        z=1/float(a)-1
        return self.epsilon0 + (self.epsilona*(a-1)+self.epsilonz*z)*self.nu(a)+self.epsilona2*(a-1)
 
    def alphaAM(self,a):
        return self.alphaAM0+(self.alphaAMa*(a-1))*self.nu(a)
    def delta(self,a):
        z=1/float(a)-1
        return self.delta0+(self.deltaa*(a-1)+self.deltaz*z)*self.nu(a)
    def gamma(self,a):
        z=1/float(a)-1
        return self.gamma0+(self.gammaa*(a-1)+self.gammaz*z)*self.nu(a)
     

    def get_sigma(self, M=None):
        if M is None:
            return 0.2
        if type(M)==list or type(M)==np.ndarray:
            return np.array([0.2]*len(M))
        else:
            return 0.2


#vmax_ach = 9.48535156 
#vmax_filt = 23.54248047 
class Moster(AMbase):
    def __init__(self,hostSHMF=False,hostAlpha=m200inf_alpha_host, hostK=m200inf_k_host, reionization=False,z=13,t_offset=0, catreion=False, vmax_ach=9.48535156 , vmax_filt= 23.54248047,  lmcSHMF=False,lmcAlpha=m200inf_alpha_lmc, lmcK=m200inf_k_lmc, no_load_reionization = False):
        self.M10 = 11.590
        self.M11 = 1.195
        self.N10 = .0351
        self.N11 = -0.0247
        self.B10 = 1.376
        self.B11 = -0.826
        self.G10 = 0.608
        self.G11 = 0.329

        self.delta_crit = 200

        # SHMF parameters. Slope alpha and constant, K
        self.alpha = m200inf_alpha_field
        self.K = m200inf_k_field
        self.color = colors.hex2color('#377eb8')   # 'blue'
        self.label = 'Moster'
        self.isPoisson=False
        self.ls = '--'
        self.z=z
        self.t_offset = t_offset
        self.superPoissonian = False
        self.catreion = catreion

        if hostSHMF:
            self.alpha= hostAlpha
            self.K = hostK
            self.superPoissonian = True

        if lmcSHMF:
            self.alpha = lmcAlpha
            self.K = lmcK
            self.superPoissonian = True


        self.infall_tck = infall_times.get_infall_tck(field = not hostSHMF)
        self.lum_tck = None

        self.reionization = reionization
        if reionization:
            #self.color = 'midnightblue'
            self.label = 'Moster + Reion'

            if catreion:
                self.lum_tck, _, _ = reion.get_lum_fraction_caterpillar(mdef='m200',z=self.z,vmax_ach=vmax_ach, vmax_filt=vmax_filt)
            else:
                self.lum_tck, _, _ = reion.get_lum_fraction_Barber()
            self.isPoisson=False
            self.ls = '-'
        elif not no_load_reionization:
            self.lum_tck, _, _ = reion.get_lum_fraction_Barber()  # need to load it anyway for mstar_to_mhalo

    def getStellarMass(self,M, a=1):
        return M*(2*self.N(a)*( (M/self.M1(a) )**-self.beta(a) + ( M/self.M1(a) )**self.gamma(a) )**-1)

    def M1(self,a):
        return 10**(self.M10+self.M11*(1-a))

    def N(self,a):
        return self.N10 + self.N11*(1-a)

    def beta(self,a):
        return self.B10 + self.B11*(1-a)

    def gamma(self,a):
        return self.G10 + self.G11*(1-a)


    # should only input reduced mass into this function
    def getStellarMass_random(self,M,a=1):
        # get random scatter, perturbation
        if type(M)==list or type(M)==np.ndarray:
            perturbation = np.array([random.gauss(0,self.sigma) for i in range(len(M))])
        else:
            perturbation = random.gauss(0,self.sigma)
        return 10**( ((M/(self.M2))**self.slope2 + (M/(self.M1))**self.slope1) + perturbation)


    ### NEW FUNCTION TO GENERATE ###
    # why is it so slow????
    # is it the append?
    # is it get stellar mass function?
    # generate dm_masses?
    # scale_factors??
    def generate_stellar_mass_samples(self,halo_mass,N=N_iter,factor=1):
        #import time
        #ts=time.time()
        #tlum=0; tinfall=0
        stellar_mass_lists=[]
    	for i in range(N):
            # generate dm_masses
            dm_masses = self.generate_shm_sample(halo_mass, factor)
            if self.reionization:  # reduce the sample
                #if i%100==0:
                #    print i
                #    print tlum, tinfall, 'lum and infall time'
                #t0=time.time()
                lum_chances = reion.frac_luminous(dm_masses, self.lum_tck) # hopefully a list of fractions between 0 and 1
                #t1=time.time()
                #tlum+=t1-t0
                randnums = np.random.rand(len(lum_chances))
                mask = randnums < lum_chances # chooses which halos will form stars
                dm_masses = dm_masses[mask] # dm masses that host stars
            
            #t0=time.time()
            scale_factors = infall_times.random_infall_scale(self.infall_tck, len(dm_masses), offset=self.t_offset) ## the luminous ones should be shifted towards earlier times
            #t1=time.time()
            #tinfall +=t1-t0
            
            stellar_mass_lists.append(self.getStellarMass(dm_masses,scale_factors))
        #tf=time.time()
        #print tf-ts, 'total time'
        #print tlum, tinfall, 'lum and infall time'
   	return np.array(stellar_mass_lists)


    def get_sigma(self, M=None):
        if M is None:
            return 0.2
        if type(M)==list or type(M)==np.ndarray:
            return np.array([0.2]*len(M))
        else:
            return 0.2


## where are the lower end extrapolation points of these models?
## Brook - 10^7 - 10^8 Mstar (confirmed)
# GK14   - 10^8 Mstar    (confirmed)
# GK16   - 4.5 * 10^5 Mstar   (confirmed)
# Moster - just above 10^7 Mstar??  From figure 3 and figure 5, extrapolation starts at 10^7.4
# Behroozi = 3e8 Mstar or 8*10^10 Mhalo for mpeak
def plotMstar_v_Mhalo():
    a = 1.0
    M = 10**np.arange(7.5,12.5,.1)
    #AMmodels=[Moster(), Behroozi(), GarrisonKimmel(), Brook(), GarrisonKimmel16(), Sawala()]
    AMmodels=[Behroozi(), GarrisonKimmel(), GK16_grow(), Moster(), Brook()]
    for model in AMmodels:
        if isinstance(model, Brook):
            # 10^7 to 10^8 Mstar
            Mlow = np.log10(model.stellar_to_halo_mass(10**7))
            Mhigh = np.log10(model.stellar_to_halo_mass(10**8))
            print Mlow, Mhigh, 'mlow and high'
            Mtmp = 10**np.arange(Mlow,Mhigh,.1)
            mstar = model.getStellarMass(Mtmp,a)
            plt.plot(Mtmp, mstar, label=model.label,lw=linewidth, color=model.color, ls='-')
            
            Mtmp = 10**np.arange(7.5,Mlow,.1)
            mstar = model.getStellarMass(Mtmp,a)
            plt.plot(Mtmp, mstar, lw=linewidth, color=model.color, ls='--')
            
            Mtmp = 10**np.arange(Mhigh,12.5,.1)
            mstar = model.getStellarMass(Mtmp,a)
            plt.plot(Mtmp, mstar, lw=linewidth, color=model.color, ls='--')

        else:

            if isinstance(model, GK16_grow):  #or isinstance(model, Sawala):
                mstar_up = model.getStellarMass_up1sig(M,a)
                mstar_down = model.getStellarMass_down1sig(M,a)
                plt.fill_between(M,mstar_up,mstar_down, facecolor=model.color, alpha=0.25)
                #mstar_rand = model.getStellarMass_random(M,a)
                #plt.scatter(M,mstar_rand,color=model.color)
                Mlow = np.log10(model.stellar_to_halo_mass(4.5*10**5))

            if isinstance(model, GarrisonKimmel):
                Mlow = np.log10(model.stellar_to_halo_mass(10**8))
            if isinstance(model, Behroozi):
                Mlow = np.log10(8*10**10)
            if isinstance(model, Moster):
                Mlow = np.log10(model.stellar_to_halo_mass(10**7.4))

            # solid line
            Mtmp = 10**np.arange(Mlow,12.5,.1)
            mstar = model.getStellarMass(Mtmp,a)
            plt.plot(Mtmp, mstar, label=model.label,lw=linewidth, color=model.color)
        
            #extrapolation
            Mtmp = 10**np.arange(7.5,Mlow,.1)
            mstar = model.getStellarMass(Mtmp,a)
            plt.plot(Mtmp, mstar,lw=linewidth, color=model.color, ls='--')          


        #else:
        #    mstar = model.getStellarMass(M,a)
        #    plt.plot(M, mstar, label=model.label,lw=linewidth, color=model.color)



    plt.legend(loc='lower right',fontsize=legend_size,frameon=False)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('$\mathrm{M_{halo}} \, \mathrm{[M_\odot]}$',fontsize=label_font)
    plt.ylabel('$\mathrm{M_*} \, \mathrm{[M_\odot]}$',fontsize=label_font)
    plt.ylim((10**3,10**11))
    plt.xlim((10**7.5,10**12.5))
    plt.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig('smhm_relation.pdf') # compare to figure 5 in paper
    plt.close()







def plotBentbyBaryons():
    #f,ax = plt.subplots(ncols=1)
    a = 1.0
    M = 10**np.arange(7.5,np.log10(5e10),.15)
    AMmodels= [Moster(), Sawala()]
    for model in AMmodels:
        mstar = model.getStellarMass(M,a)
        plt.plot(M, mstar, label=model.label,linewidth=5, color=model.color)

    plt.legend(loc='lower right',frameon=False,fontsize=26)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('$M_{halo} \mathrm{[M_\odot]}$',fontsize=30)
    plt.ylabel('$M_* \mathrm{[M_\odot]}$',fontsize=30)
    plt.ylim((10**3,10**8))
    plt.xlim((10**8,5e10))
    plt.tick_params(axis='both', which='major', labelsize=21)
    #plt.gca().tight_layout()	
    plt.gcf().subplots_adjust(bottom=0.17)
    plt.gcf().subplots_adjust(left=0.17)
    plt.savefig('LMCPlots/BentByBaryons_peak') # compare to figure 4 in bent by baryons
    plt.close()

def plotBentbyBaryons2():
    #f,ax = plt.subplots(ncols=1)
    pg = -1.0
    a = 1.0
    M = 10**np.arange(7.5,np.log10(5e10),.05)
    AMmodels= [Sawala(), GK16_grow(bent=True, plotgamma = pg)]
    AMlabels = ['Sawala AM model', 'GK16 with\nslope = 3.31 and bend']
    for model, label in zip(AMmodels, AMlabels):
        mstar = model.getStellarMass(M,a)
        plt.plot(M, mstar, label=label,linewidth=5, color=model.color)
        if isinstance(model, GK16_grow):  #or isinstance(model, Sawala):
            mstar_up = model.getStellarMass_up1sig(M,a)
            mstar_down = model.getStellarMass_down1sig(M,a)
            plt.fill_between(M,mstar_up,mstar_down, facecolor=model.color, alpha=0.25)

    plt.legend(loc='upper left',frameon=False,fontsize=legend_size)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('$M_{\mathrm{peak}} \mathrm{[M_\odot]}$',fontsize=label_font)
    plt.ylabel('$M_* \mathrm{[M_\odot]}$',fontsize=label_font)
    plt.ylim((10**2,10**8))
    plt.xlim((10**7.5,4e10))
    plt.tick_params(axis='both', which='major', labelsize=tick_size)
    #plt.gca().tight_layout()	
    plt.gcf().subplots_adjust(bottom=0.14)
    plt.gcf().subplots_adjust(left=0.14)
    plt.savefig('LMCPlots/BentByBaryons2') # compare to figure 4 in bent by baryons
    plt.close()



# Figure 1 Plot
#plotMstar_v_Mhalo()
