import numpy as np
import haloutils
from scipy.integrate import quad
import MTanalysis_field as mtaf
import MTanalysis3 as mta
import MTaddition as mtadd
import MTaddition_field as mtaddf
import pandas
from scipy.optimize import fsolve 

import abundance_matching as am

# return contamination distance in Mpc/h
def get_contam_dist(hpath):
    # create dictionary
    hpaths, mindists = np.load('/nfs/blank/h4231/gdooley/Dropbox/DwarfsOfDwarfs/code/contam_dists.npy')
    look_up = dict(zip(hpaths, mindists.astype('float')))
    return look_up[hpath]

def generate_contam_dist():
    import haloutils
    #hpaths = haloutils.get_all_halo_paths_lx(lx=14)
    catnums = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,29,31,33,36,37,40,53]
    hpaths = [haloutils.catnum_hpath(i,14) for i in catnums]   
    mindists=[]
    for hpath in hpaths:
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        posns = haloutils.load_partblock(hpath,snap_z0,'POS ',parttype=2)
        hostID = haloutils.load_zoomid(hpath)
        cat=haloutils.load_rscat(hpath,snap_z0,rmaxcut=True)
        dists = distance(np.array(cat.ix[hostID][['posX','posY','posZ']]), posns)
        mindists.append(np.min(dists))
        print np.min(dists)
    print hpaths
    print mindists
    np.save('contam_dists',np.array([hpaths,mindists]))

# return virial radii in kpc
def get_rvir(hpath):
    # create dictionary
    hpaths, mindists = np.load('/nfs/blank/h4231/gdooley/Dropbox/DwarfsOfDwarfs/code/rvirs.npy')
    look_up = dict(zip(hpaths, mindists.astype('float')))
    return look_up[hpath]

def generate_host_rvirs():
    rvirs=[]
    hpaths = haloutils.get_all_halo_paths_lx(lx=14)[0:30]
    for hpath in hpaths:
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        hostID = haloutils.load_zoomid(hpath)
        cat=haloutils.load_rscat(hpath,snap_z0,rmaxcut=True)
        rvirs.append(cat.ix[hostID]['rvir']/cat.h0)
        print cat.ix[hostID]['rvir']/cat.h0
    np.save('rvirs',np.array([hpaths,rvirs]))


def contam_ratios():
    hpaths = haloutils.get_all_halo_paths_lx(lx=14)[0:30]
    ratios = []
    for hpath in hpaths:
        contam = get_contam_dist(hpath)*1000/0.67
        rvir = get_rvir(hpath)
        ratios.append(contam/rvir)
    return np.array(ratios)




# compute distance from posA to posB.       
 # posA can be an array. boxsize must be in same units as positions. 
def distance(posA, posB,boxsize=100.):
    dist = abs(posA-posB)
    tmp = dist > boxsize/2.0
    dist[tmp] = boxsize-dist[tmp]
    if dist.shape == (3,):
        return np.sqrt(np.sum(dist**2))
    else:
        return np.sqrt(np.sum(dist**2,axis=1))

def getMidpoints(bins):
    """               
    Given a range of cut-off values, find the midpoints.
    @return: array of length len(bins)-1 
    """
    spacing = bins[1:]-bins[:-1]
    return bins[:-1]+spacing/2.0



# load dataE but eliminate cases of halos that accreted into the host
# on their own, but then merged with a subhalo that was already in the
# host.
"""
these are halos that fell into the host on their own
then entered the virial radius of a subhalo
then merged with the subhalo which is now an extant halo
and they are perfectly good halos to tag, but I can only find them traversing the extant tree
so unless you care about particle tagging at infall, you should add to your dataE reader a clause to mask out any halos with ['depth'] > 0
or ['backsnap'] < 319
or ['sub_rank'] < 0
they all give info on these halos, and the repeated ['rsid'] allows you to know what these guys merged into
"""

# appears halo 10 does not have up-to-date extant data run on it.
# halo 18 does not have any field halos
def get_extant_data(hpath,field=False):
    if not field:
        AE = mta.AllExtantData()
        dataE = AE.read(hpath)
        #mask = dataE['depth'] == 0
        #dataE = dataE[mask]
        hostID = haloutils.load_zoomid(hpath)
        pidmask = dataE['pid'] == hostID
        #print np.sum(pidmask), np.sum(~pidmask), float(np.sum(pidmask))/len(pidmask)
        
        # cat-7 is messed up!! double check the hostID on this one. only 17% are subs
        # more are subhalos of 45758 than of 45347. must have been a recent accretion
        # do not include in the host halo sample!
        # also 40 and 53

        E = mtadd.ExtantDataReionization()
        data_reion = E.read(hpath)
        #maskR = data_reion['depth_reion'] == 0 
        #data_reion = data_reion[maskR]
        return pandas.concat((dataE[pidmask] ,data_reion[pidmask] ),axis=1) 

    if field:
        #FieldData = mtaf.FieldHaloSubstructureFirstPass()
        #dataE = FieldData.read(hpath)
        AEF = mtaf.AllExtantFieldData()
        dataE = AEF.read(hpath)

        mask = dataE['depth'] == 0 
        E = mtaddf.ExtantDataReionizationField()
        data_reion = E.read(hpath)
        return pandas.concat((dataE[mask] ,data_reion[mask] ),axis=1) 
        #return dataE[mask]


#9.4 and up, make it 100%. 6.8 and below, make it 0


# cosmological values from caterpillar
def Hubble(a,Om=.3175,Ol=.6825,h=0.6711):
    to_inv_sec = 3.241*10**-20
    Ok = 1-Om-Ol
    return to_inv_sec * 100*h*(Om*(1/a)**3 + Ok*(1/a)**2 + Ol)**.5 # value in 1/sec

"""
def H(z,Om,Ol,h):
    Ok = 1-Om-Ol
    return (Om*(1+z)**3 + Ok*(1+z)**2 + Ol)**.5
"""

def lookback_func(z,Om,Ol):
    return 1./((1+z)*np.sqrt(Om*(1+z)**3 + (1-Om-Ol)*(1+z)**2 + Ol))

def lookback(a,Om=0.3175,Ol=0.6825,h=0.6711):
    #HubbleConst = 3.421 * 10**-18 #inverse seconds for h = 1.0
    z = 1./a - 1.
    times = np.zeros((len(z)))
    for i in range(len(z)):
        times[i] = quad(lambda x: lookback_func(x,Om,Ol),0,z[i])[0]
    return times/(Hubble(1))*3.171e-17


def time_to_scale(time,Om=0.3175,Ol=0.6825,h=0.6711):
    def func(z):
        return time - quad(lambda x: lookback_func(x,Om,Ol),0,z)[0]/(Hubble(1))*3.171e-17  
    z = fsolve(func,1)[0]
    a = 1./(1.+z)
    return a



# cat-7 is messed up!! double check the hostID on this one. only 17% are subs
# more are subhalos of 45758 than of 45347. must have been a recent accretion
# do not include in the host halo sample!
# also 40 and 53

def get_hpaths(field, lx=14):
    if field:
        # remove halo 18 - contaminated
        catnums = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,19,20,21,22,23,24,25,26,27,29,31,33,36,37,40,53]  
        # 34 has no rockstar
        # 11 is messed up
        # cat-7 is messed up!! only 17% are subs
        # more are subhalos of 45758 than of 45347. must have been a recent accretion
        # do not include in the host halo sample!
        # also 40 and 53
    else:
        # remove halo 11 temporarily
        catnums = [1,2,3,4,5,6,8,9,10,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,29,31,33,36,37,40,53]
        #catnums = [36,37,53]  # 34 not run, 40 not run, 53 not run. 34 no rockstar
        # 11, 40, 53 running, should be added later
    hpaths = [haloutils.catnum_hpath(i,lx) for i in catnums]
    return hpaths



def rho_enlcNWF(r,rs,rho_0):
    mencl = rho_0*4*np.pi*rs**3 *( np.log((rs+r)/rs) - r/(r+rs) )
    return mencl / (4/3. * np.pi * r**3)


# convert from bryan and norman z=0 mass to mcrit_delta.
# standard delta is 103.86 for Bryan and Norman
def convert_Mhalo_z0(Mvir, delta):  # convert to the stated delta
    G = 4.157e-39  # kpc^3/Msun/s^2
    H = Hubble(a=1) # in 1/s
    rho_crit = (3*H**2)/(8*np.pi*G) # in Msun/kpc^3
    rvir = (Mvir*3/(4*np.pi*103.86*rho_crit))**(1/3.) # in kpc
    #conc = 7.5
    conc = 90.1 * Mvir**(-1./11) # just an ok guess
    rs = rvir/conc
    rho_0 = Mvir/ (4*np.pi*rs**3 *( np.log((rs+rvir)/rs) - rvir/(rvir+rs) )) # msun/kpc^3
    
    
    from scipy.optimize import fsolve    
    def funcNFW(x):
        return rho_enlcNWF(x,rs,rho_0)/(rho_crit) - delta
    rvirNFW = fsolve(funcNFW,1)
    mvirNFW = 4./3 * np.pi * rvirNFW**3 * rho_crit * delta
    #print rvir, rvirNFW, 'rvir350 and rvir350 NFW'
    return mvirNFW[0]




# convert from delta1 to delta2. delta1 is like 103.86
# delta2 is like delta in the previous function
def convert_Mhalo_d_d(Mvir, delta1, delta2):  # convert to the stated delta
    G = 4.157e-39  # kpc^3/Msun/s^2
    H = Hubble(a=1) # in 1/s
    rho_crit = (3*H**2)/(8*np.pi*G) # in Msun/kpc^3
    rvir = (Mvir*3/(4*np.pi*delta1*rho_crit))**(1/3.) # in kpc

    if delta1 == delta2:
        return Mvir, rvir

    #conc = 7.5
    conc = 90.1 * Mvir**(-1./11) # just an ok guess
    rs = rvir/conc
    rho_0 = Mvir/ (4*np.pi*rs**3 *( np.log((rs+rvir)/rs) - rvir/(rvir+rs) )) # msun/kpc^3
        
    from scipy.optimize import fsolve    
    def funcNFW(x, rscale, rho):
        return rho_enlcNWF(x,rscale,rho)/(rho_crit) - delta2

    if type(Mvir)==list or type(Mvir)==np.ndarray:
        rvirNFW = np.array([fsolve(funcNFW,1, args=(rscale, rho))[0] for rscale, rho in zip(rs, rho_0) ])
        mvirNFW = 4./3 * np.pi * rvirNFW**3 * rho_crit * delta2
        return mvirNFW, rvirNFW    
    else:
        rvirNFW = fsolve(funcNFW,1,args=(rs, rho_0))
        mvirNFW = 4./3 * np.pi * rvirNFW**3 * rho_crit * delta2
        #print rvir, 'rvir original'
        #print rvirNFW[0], 'rvir new'
        #print rvir, rvirNFW, 'rvir350 and rvir350 NFW'
        return mvirNFW[0], rvirNFW[0]
    


def load_nearby_gx(v=2):
    if v==1:
        data = np.loadtxt('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/LMCPlots/nearby_gx_short_dMCs(v1).dat', dtype='str')
    else:
        data = np.loadtxt('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/LMCPlots/nearby_gx_short_dMCs.dat', dtype='str')
    dtype = [('Name', np.str_,17), ('RA', np.float), ('Dec', np.float), ('dist', np.float), ('dist_GC', np.float), ('dist_LMC', np.float), ('dist_SMC', np.float), ('mstar', np.float) ]
    arr = np.ndarray((len(data)-2,),dtype=dtype)
    arr['Name'] = np.array(data[2:,0], dtype='string')
    arr['RA'] = np.array(data[2:,1],dtype='float')
    arr['Dec'] = np.array(data[2:,2],dtype='float')
    arr['dist'] = np.array(data[2:,3],dtype='float')
    arr['dist_GC'] = np.array(data[2:,4],dtype='float')
    arr['dist_LMC'] = np.array(data[2:,5],dtype='float')
    arr['dist_SMC'] = np.array(data[2:,6],dtype='float')
    arr['mstar'] = np.array(data[2:,7], dtype='float')
    NoCanis = arr['Name'] != '*Canis_Major'
    return arr[NoCanis]   # exclude Canis Major because it is highly disputed. Likely not a galaxy.
 
    

# for plotting a step plot of y vs. x
def step_plotX(x):
    return np.repeat(x,2)[1:]
def step_plotY(y):
    return np.repeat(y,2)[0:-1]


def convert_Mvir(Mvir, model):
    if isinstance(model, am.Moster) or isinstance(model, am.Sawala):
        return convert_Mhalo_z0(Mvir, 200)
    if isinstance(model, am.GarrisonKimmel16) or isinstance(model, am.GarrisonKimmel) or isinstance(model, am.Behroozi) or isinstance(model, am.GK16_grow):
        return Mvir
    if isinstance(model, am.Brook):
        return convert_Mhalo_z0(Mvir, 350)
