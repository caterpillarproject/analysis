import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import abundance_matching as am
import AMbase as ambase
import DwarfMethods as dm
import haloutils as htils

# he uses mvir at z=0.
# I need to translate from mvir at z=0 to mpeak, m350 peak, and minfall.
# just like that other time.
# given msub/mhost in mvir z=0 form, how do I translate it to m200 infall, etc?
# do that sampling... take all m200 infall, compute mvir z=0 / mhost z=0 for them.
# there was huge scatter, right?

#model = am.GK16_grow(lmcSHMF=True)
#mu = np.logspace(-4,-1,20)
#mhost = 2*10**11
#msub = mu*Mhost

def get_sig_ratio(Msub, Mhost, mvir_alpha =1.90, mvir_k = 0.00122202932021, plotit = False):
    #msub = mu*Mhost
    mu = Msub/Mhost
    #print mu, 'mu in get ratio'

    N_mean = ambase.Nsubs_bigger(Msub, Mhost,mvir_alpha,  mvir_k)    #model.alpha, model.K)
    s1 = 0.14  #0.18
    sigmaI = N_mean * s1
    sigmaP2 = N_mean
    sigma2 = sigmaP2 + sigmaI**2
    sigma_ratio = np.sqrt(sigma2)/np.sqrt(sigmaP2)
    
    sigma_ratio2 = np.sqrt( N_mean + s1**2*N_mean**2 ) / np.sqrt(N_mean)
    #if not (sigma_ratio == sigma_ratio2).all():
    #    print 'mismatched sigmas', sigma_ratio, sigma_ratio2

    #print sigma_ratio, 'sig ratio for mu'
    if plotit:
        plt.plot(mu, sigma_ratio)
        plt.xscale('log')
        plt.ylabel('$\sigma / \sigma_{POISSON}$')
        plt.xlabel('$\mu \equiv M_{sub}/M_{vir}$')
        plt.savefig('LMCPlots/variance.png')
        plt.close()
    return sigma_ratio


### MUST MODIFY THIS TO MY VARIANCE NEEDS ###
def convert_masses_for_variance():
    mvirz0 = [];  ng_mvirz0=[]
    m350=[];  ng_m350=[]
    mpeak = [];   ng_mpeak = []
    m200 = []; ng_m200 = []
    sig_ratio = []
    hpaths = dm.get_hpaths(False)
    mu = np.logspace(-4,-1,20)

    for hpath in hpaths:
        snap_z0 = htils.get_numsnaps(hpath)-1
        cat=htils.load_rscat(hpath,snap_z0,rmaxcut=False)
        data = dm.get_extant_data(hpath,False)
        z0mass = np.array(cat.ix[data['rsid']]['mgrav']/cat.h0)
        hostID = htils.load_zoomid(hpath)
        mhost = float(cat.ix[hostID]['mgrav']/cat.h0 )
        sig_ratio.append(  get_sig_ratio(z0mass, mhost) )

        m350.append( np.array(data['max_mass350NFW']/cat.h0) / dm.convert_Mhalo_z0(mhost, 350) )  # needs to be divided by m350 of the host
        mpeak.append( np.array(data['max_mass']/cat.h0) / mhost )
        m200.append( np.array(data['infall_mass200']/cat.h0 / dm.convert_Mhalo_z0(mhost, 200)))
        mvirz0.append(z0mass / mhost)

        ng_m350.append( np.array([np.sum(m350[-1] > mu[i]) for i in range(len(mu)) ]) )
        ng_mpeak.append( np.array([np.sum(mpeak[-1] > mu[i]) for i in range(len(mu)) ]) )
        ng_m200.append( np.array([np.sum(m200[-1] > mu[i]) for i in range(len(mu)) ]) )
        ng_mvirz0.append( np.array([np.sum(mvirz0[-1] > mu[i]) for i in range(len(mu)) ]) )
        print 'done with cat', htils.hpath_catnum(hpath)



    host_mass = 1.4*10**12
    hm = dm.convert_Mhalo_z0(host_mass, 350)
    # this fit is not very good. for m350 shmf.  max_mass350NFW is good for both. the host mass must be in m350 mass
    m350_mean = ambase.Nsubs_bigger(mu*hm,hm,am.m350peak_alpha_host, am.m350peak_k_host)
    m350_mean_data = np.mean(ng_m350,axis=0)
    m350_std = np.std(ng_m350,axis=0)
    m350_ratio = m350_std / np.sqrt(m350_mean_data)

    # this now fits wonderfully
    mpeak_mean = ambase.Nsubs_bigger(mu*host_mass,host_mass,am.mpeak_alpha_host, am.mpeak_k_host)
    mpeak_mean_data = np.mean(ng_mpeak,axis=0)
    mpeak_std = np.std(ng_mpeak,axis=0)
    mpeak_ratio = mpeak_std / np.sqrt(mpeak_mean_data)

    hm = dm.convert_Mhalo_z0(host_mass, 200)
    m200_mean = ambase.Nsubs_bigger(mu*hm,hm,am.m200inf_alpha_host, am.m200inf_k_host)
    m200_mean_data = np.mean(ng_m200,axis=0)
    m200_std = np.std(ng_m200,axis=0)
    m200_ratio = m200_std / np.sqrt(m200_mean_data)

    mvirz0_mean = ambase.Nsubs_bigger(mu*host_mass,host_mass,am.mvir_alpha_host, am.mvir_k_host)
    mvirz0_mean_data = np.mean(ng_mvirz0,axis=0)
    mvirz0_std = np.std(ng_mvirz0,axis=0)
    mvirz0_ratio = mvirz0_std / np.sqrt(mvirz0_mean_data)


    m350 = np.array([item for arr in m350 for item in arr])
    mpeak = np.array([item for arr in mpeak for item in arr])
    m200 = np.array([item for arr in m200 for item in arr])
    mvirz0 = np.array([item for arr in mvirz0 for item in arr])
    sig_ratio = np.array([item for arr in sig_ratio for item in arr])  # this is now completely messed up

    np.save('m350_variance', m350)
    np.save('mpeak_variance', mpeak)
    np.save('m200_variance', m200)
    np.save('mvirz0_variance', mvirz0)
    np.save('sig_ratio_variance', sig_ratio)

    np.save('m350_ratio', m350_ratio)
    np.save('mpeak_ratio', mpeak_ratio)
    np.save('m200_ratio', m200_ratio)
    np.save('mvirz0_ratio', mvirz0_ratio)



def plot_things():
    mu = np.logspace(-4,-1,20)
    m350_mu = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/m350_variance.npy')
    mpeak_mu = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/mpeak_variance.npy')
    m200_mu = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/m200_variance.npy')
    mvirz0_mu = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/mvirz0_variance.npy')
    sig_ratio = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/sig_ratio_variance.npy')

    m350_ratio = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/m350_ratio.npy')
    mpeak_ratio = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/mpeak_ratio.npy')
    m200_ratio = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/m200_ratio.npy')
    mvirz0_ratio = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/mvirz0_ratio.npy')

    midpoints = 10**dm.getMidpoints(np.log10(mu))

    """
    masks = [ (m350_mu > mu[i])&(m350_mu < mu[i+1]) for i in range(len(mu)-1)]
    y_values = [np.mean(sig_ratio[mask]) for mask in masks]
    plt.plot(midpoints, y_values, linestyle='-', label='m350', color='green')

    masks = [ (mvirz0_mu > mu[i])&(mvirz0_mu < mu[i+1]) for i in range(len(mu)-1)]
    y_values = [np.mean(sig_ratio[mask]) for mask in masks]
    plt.plot(midpoints, y_values, linestyle='-', label='mvir z0', color='red')

    masks = [ (mpeak_mu > mu[i])&(mpeak_mu < mu[i+1]) for i in range(len(mu)-1)]
    y_values = [np.mean(sig_ratio[mask]) for mask in masks]
    plt.plot(midpoints, y_values, linestyle='-', label='mpeak', color='blue')

    masks = [ (m200_mu > mu[i])&(m200_mu < mu[i+1]) for i in range(len(mu)-1)]
    y_values = [np.mean(sig_ratio[mask]) for mask in masks]
    plt.plot(midpoints, y_values, linestyle='-', label='m200', color='orange')
    """

    
    mhost = 1.4*10**12
    plt.plot(mu,  get_sig_ratio(mhost*mu, mhost, plotit = False), linestyle='--', label='mvir z0 analytic mw', color='red')
    model = am.GK16_grow(hostSHMF=True)    
    plt.plot(mu,  get_sig_ratio(mhost*mu, mhost, mvir_alpha = model.alpha, mvir_k=model.K, plotit = False), linestyle='--', label='mpeak analytic mw', color='blue')
    model = am.Moster(hostSHMF=True)
    hm = dm.convert_Mhalo_z0(mhost, 200)    
    plt.plot(mu,  get_sig_ratio(hm*mu, hm, mvir_alpha = model.alpha, mvir_k=model.K, plotit = False), linestyle='--', label='m200 infall analytic mw', color='orange')
    model = am.Brook(hostSHMF=True)
    hm = dm.convert_Mhalo_z0(mhost, 350)
    plt.plot(mu,  get_sig_ratio(hm*mu, hm, mvir_alpha = model.alpha, mvir_k=model.K, plotit = False), linestyle='--', label='m350 peak analytic mw', color='green')


    plt.scatter(mu, m350_ratio, color='green',marker='x')
    plt.scatter(mu, mpeak_ratio, color='blue',marker='x')
    plt.scatter(mu, m200_ratio, color='orange',marker='x')   
    plt.scatter(mu, mvirz0_ratio, color='red',marker='x')

    plt.legend()
    plt.xscale('log')
    plt.xlim((10**-4, 10**-1))
    #plt.ylim((1,4.0))
    plt.savefig('LMCPlots/test_mass_variance')
    plt.close()



def test_generate_dm_masses_brook():
    m350_ratio = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/m350_ratio.npy')
    mhost = 1.4*10**12
    model = am.Brook(hostSHMF=True)
    halo_mass = dm.convert_Mhalo_z0(mhost, 350)
    mu = np.logspace(-4,-1,20)
    ng = []; ng_super=[]
    for i in range(4000):    
        dm_masses = ambase.generate_shm_sample_static(halo_mass,model.alpha,model.K, superPoisson=False)
        ng.append( np.array([np.sum(dm_masses/halo_mass > mu[i]) for i in range(len(mu)) ]) )
        dm_masses2 = ambase.generate_shm_sample_static(halo_mass,model.alpha,model.K, superPoisson=True)
        ng_super.append( np.array([np.sum(dm_masses2/halo_mass > mu[i]) for i in range(len(mu)) ]) )
    mean_data = np.mean(ng,axis=0)
    std = np.std(ng,axis=0)
    ratio = std / np.sqrt(mean_data)


    mean_data2 = np.mean(ng_super,axis=0)
    std2 = np.std(ng_super,axis=0)
    ratio2 = std2 / np.sqrt(mean_data2)

    plt.plot(mu, ratio, label='poisson generated masses')
    plt.plot(mu, ratio2, label='NBD generated masses')

    plt.scatter(mu, m350_ratio, color='green',marker='x')
    plt.plot(mu,  get_sig_ratio(halo_mass*mu, halo_mass, mvir_alpha = model.alpha, mvir_k=model.K, plotit = False), linestyle='--', label='m350 peak analytic mw', color='green')


    plt.legend()
    plt.xscale('log')
    plt.xlim((10**-4, 10**-1))

    plt.savefig('LMCPlots/generate_dm_masses_test')
    plt.close()


def test_generate_dm_masses_GK16():
    mpeak_ratio = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/mpeak_ratio.npy')
    mhost = 1.4*10**12
    model = am.GK16_grow(hostSHMF=True)
    halo_mass = mhost
    mu = np.logspace(-4,-1,20)
    ng = []; ng_super=[]
    for i in range(4000):    
        #dm_masses = ambase.generate_shm_sample_static(halo_mass,model.alpha,model.K, superPoisson=False)
        #ng.append( np.array([np.sum(dm_masses/halo_mass > mu[i]) for i in range(len(mu)) ]) )
        dm_masses2 = ambase.generate_shm_sample_static(halo_mass,model.alpha,model.K, superPoisson=True)
        ng_super.append( np.array([np.sum(dm_masses2/halo_mass > mu[i]) for i in range(len(mu)) ]) )
    #mean_data = np.mean(ng,axis=0)
    #std = np.std(ng,axis=0)
    #ratio = std / np.sqrt(mean_data)


    mean_data2 = np.mean(ng_super,axis=0)
    std2 = np.std(ng_super,axis=0)
    ratio2 = std2 / np.sqrt(mean_data2)

    #plt.plot(mu, ratio, label='poisson generated masses')
    plt.plot(mu, ratio2, label='NBD generated masses')

    plt.scatter(mu, mpeak_ratio, color='blue',marker='x')
    plt.plot(mu,  get_sig_ratio(halo_mass*mu, halo_mass, mvir_alpha = model.alpha, mvir_k=model.K, plotit = False), linestyle='--', label='mpeak peak analytic mw', color='blue')


    plt.legend()
    plt.xscale('log')
    plt.xlim((10**-4, 10**-1))

    plt.savefig('LMCPlots/generate_dm_masses_test_GK16')
    plt.close()
