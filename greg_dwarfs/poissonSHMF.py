import numpy as np
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import haloutils as htils  # must be imported after matplotlib without xterm on
import scipy.misc
from field_shmf import get_field_halos


def dndM(Msub,Mhost,alpha=1.895046,K=0.00102):
    return K*Msub**(-alpha) * Mhost

def N(mlow,mhigh,Mhost,alpha=1.895046,K=0.00102):
    return K*Mhost/(alpha-1) * (mlow**(1-alpha) - mhigh**(1-alpha) )

def get_mlow(Mhost, Mhost2,f,mlow1,alpha=1.895046,K=0.00102):
    return ( Mhost/Mhost2 * (mlow1**(1-alpha) - (f*Mhost)**(1-alpha)) + (f*Mhost2)**(1-alpha) )**(1/(1-alpha))


# now want to test poisson distribution
K = .000703; alpha = 1.866887
f=0.1; mlow = 10**9.5; Mhost=10**12
mean = N(mlow,f*Mhost,Mhost,alpha,K)
print mean, 'mean number of occurences'
plt.figure()
ax = plt.subplot(111)

lx=14
hpaths = htils.get_all_halo_paths_lx(lx)[0:30]
nsubs_list=[]
for hpath in hpaths:   # I will have 30 data points. enough for poisson distribution?
    snap_z0 = htils.get_numsnaps(hpath)-1
    cat=htils.load_rscat(hpath,snap_z0,rmaxcut=True)
    hostID = htils.load_zoomid(hpath)
    field_halos = get_field_halos(cat,hpath,hostID,mlow=10,mhigh=11.5)
    rsid_list = np.array(field_halos['id'])
    # loop over field halos here too. by rsid, concatenated with host rsid
    for rsid in rsid_list:
        masses_sub = np.array(cat.get_subhalos_within_halo(rsid)['mgrav']/cat.h0)
        Mhost_cat = cat.ix[rsid]['mgrav']/cat.h0
        mhigh_cat = f*Mhost_cat
        mlow_cat = get_mlow(Mhost, Mhost_cat,f,mlow,alpha,K)
        nsubs = np.sum((masses_sub < mhigh_cat)&(masses_sub>mlow_cat))
        print nsubs
        nsubs_list.append(nsubs)
    print 'done with cat', htils.hpath_catnum(hpath)

print np.mean(nsubs_list), 'mean of nsubs caterpillar'
print len(nsubs_list), 'number of samples'
plt.hist(np.array(nsubs_list),bins=np.arange(0,np.max(nsubs_list)+2), normed=False,histtype='step',color='red',linewidth=1.5)


# plot pure poisson distribution
nocc = np.arange(0,mean*3+1)
prob = (mean**nocc * np.e**-mean /scipy.misc.factorial(nocc))*len(nsubs_list)
tmp=np.repeat(nocc,2)
ax.plot(np.append(tmp[1:],tmp[-1]+1), np.repeat(prob,2),linewidth=1.5)
plt.vlines(mean,0,np.max(prob))
plt.ylabel('Probability')
plt.xlabel('Num Occurences')


plt.savefig('poisson_%d' %(int(mean)))
plt.close()



"""
Mhost1 = 10**12
Mhost2 = 10**11
f=0.1
mhigh1=f*Mhost1; mhigh2=f*Mhost2
mlow1 = 10**10
print N(mlow1,mhigh1,Mhost1), 'N interval 1'
mlow2 = get_mlow(Mhost1,Mhost2,f,mlow1)
print np.log10(mlow2), 'mlow2', np.log10(mlow1), 'low1'
print N(mlow2,mhigh2,Mhost2), 'N interval 2'
"""
