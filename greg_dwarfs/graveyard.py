#### I SHOULD MAKE THIS AN inherited method ####

# want functions:
# given halo mass, what is stellar mass? 
'moster, behroozi, garrison-kimmel'
def halo_to_stellar_mass(mvir,model,a=1):
    if type(mvir)==list:
        mvir = np.array(mvir)
    return model.getStellarMass(mvir,a)


# given stellar mass/luminosity what is DM halo mass?
def stellar_to_halo_mass(mstar,model,a=1):
    #USE INTERPOLATION     
    from scipy import interpolate
    mass_array = 10**np.arange(7,12,.07)
    mstar_array = halo_to_stellar_mass(mass_array,model,a)
    tck = interpolate.splrep(np.log10(mstar_array),np.log10(mass_array), k=1)
    return 10**interpolate.splev(np.log10(mstar),tck)




#### I SHOULD MAKE THIS AN inherited method that gets over ridden with
# special AM models - garrison kimmel and sawala  ####
def P_at_least_one(halo_masses,min_lum,model):
    if isinstance(model, am.GarrisonKimmel16):
        frac_gte_one = []
        for halo_mass in halo_masses:
            samples = generate_stellar_mass_samples(halo_mass,model,N=10000) # generate list of lists of stellar masses
            lums = samples/ML_ratio
            Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
            frac_gte_one.append(np.sum(Ngreater>0) / float(len(Ngreater)))
        Pgt1 = np.array(frac_gte_one)
    else:
        lowest_mass = am.stellar_to_halo_mass(ML_ratio*min_lum,model=model)
        mean = Nsubs_bigger(lowest_mass,halo_masses,model.alpha,model.K)
        Pgt1 = 1 - np.e**-mean
    return Pgt1



##### MOVE THIS TO SUPER CLASS ######
# generate list of lists of stellar masses
def generate_stellar_mass_samples(halo_mass,model,N):
    stellar_mass_lists=[]
    for i in range(N):
        dm_masses = generate_shm_sample(halo_mass,model.alpha,model.K)
        stellar_mass_lists.append(model.getStellarMass_random(dm_masses,a=1.0))
    return np.array(stellar_mass_lists)


#### I SHOULD MAKE THIS AN inherited method that gets over ridden with
# special AM models - garrison kimmel and sawala  ####
def ngreater(halo_masses,min_lum,model):
    if isinstance(model, am.GarrisonKimmel16):
        num_gte_one = []
        for halo_mass in halo_masses:
            samples = generate_stellar_mass_samples(halo_mass,model,N=10000) # generate list of lists of stellar masses
            lums = samples/ML_ratio
            Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
            num_gte_one.append(np.sum(Ngreater) / float(len(Ngreater)))
        N_visible = np.array(num_gte_one)

    
    else:
        lowest_mass = am.stellar_to_halo_mass(ML_ratio*min_lum,model=model)
        N_visible = Nsubs_bigger(lowest_mass, halo_masses,model.alpha,model.K)
    return N_visible


#### SHOULD BE ABLE TO REMOVE LATER ####
# subhalo mass function next. Make alpha positive
def Nsubs_bigger(msub, Mhost,alpha,K):
    return K * Mhost /(alpha-1) * (msub**(1-alpha) - Mhost**(1-alpha))



### MOVE THIS TO SUPER CLASS ###
# generate a list of halo masses from mlow=10**8 to mhigh = Mhost that follows the mean of the SHMF
def generate_shm_sample(Mhost,alpha,K): # lower limit of where stars can form. minihalo analysis
    mlow = 10**8.0
    mean = Nsubs_bigger(mlow,Mhost,alpha,K)
    nsubs = np.random.poisson(mean)
    randnums = np.random.rand(nsubs)
    masses = (mlow**(1-alpha) - (alpha-1)/(K*Mhost) * randnums*mean)**(1/(1-alpha))
    return masses
