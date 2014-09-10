import numpy as np
import pylab as plt
import haloutils

def _radial_profile(hpos,subs,rbins,rvir,getmeandens=False):
    hpos = np.array(hpos)
    subpos = np.array(subs[['posX','posY','posZ']])
    dist = np.sqrt(np.sum((subpos-hpos)**2,axis=1))
    h,x = np.histogram(dist,rbins)
    dV = (4./3. * np.pi * ((rbins[1:])**3 - (rbins[:-1])**3))
    n_density = h/dV
    mean_ndens = np.sum(h)/(4./3. * np.pi*rvir**3)
    if getmeandens: return n_density/mean_ndens,mean_ndens
    return n_density/mean_ndens
    
def radial_profile(rsid,rscat,rbins=None,nbins=10):
    hrow = rscat.ix[rsid]
    hpos = hrow[['posX','posY','posZ']]
    subs = rscat.get_subhalos_from_halo(rsid)
    rvir = hrow['rvir']/1000.
    if rbins==None:
        rbins = np.linspace(0,1,nbins+1)*rvir
    rbinmid = (rbins[1:]+rbins[:-1])/2.0/rvir
    n_density_rel = _radial_profile(hpos,subs,rbins,rvir)
    return rbinmid,n_density_rel
    
def radial_profile_bin(rsid,rscat,Mbins,rbins=None,nbins=10,usetotnorm=False):
    Mbins = np.array(Mbins)
    hrow = rscat.ix[rsid]
    hpos = hrow[['posX','posY','posZ']]
    subs = rscat.get_subhalos_from_halo(rsid)
    rvir = hrow['rvir']/1000.
    if rbins==None:
        rbins = np.linspace(0,1,nbins+1)*rvir
    rbinmid = (rbins[1:]+rbins[:-1])/2.0
    submass = np.array(subs['mvir']/rscat.h0)
    nlist = []
    for i in range(len(Mbins)-1):
        mask = (submass >= Mbins[i]) & (submass < Mbins[i+1])
        nlist.append(_radial_profile(hpos,subs[mask],rbins,rvir,getmeandens=usetotnorm))
    if usetotnorm:
        totndens = 0
        nlistnew = []
        for n,mean_ndens in nlist:
            totndens += mean_ndens
            nlistnew.append(n*mean_ndens)
        nlist = []
        for n in nlistnew:
            nlist.append(n/totndens)
    return rbinmid/rvir,nlist

if __name__=="__main__":
    """
    Example of how to use these functions
    """
    lx = 14; nv = 4; ictype = "BB"
    haloids = [1327707]
    rsids = [51132]
    Mbins = [10**5,10**6,10**7,10**8,10**9,10**10,10**11,10**12]

    # Normalizing by <n>
    fig,ax=plt.subplots()
    labellist = [r"$10^{%.1f}-10^{%.1f} M_\odot$" % (np.log10(M1),np.log10(M2)) for M1,M2 in zip(Mbins[:-1],Mbins[1:])]
    for haloid,rsid in zip(haloids,rsids):
        hpath = haloutils.get_hpath(haloid,ictype,lx,nv)
        rscat = haloutils.load_rscat(hpath,255)
        xtot,ntot = radial_profile(rsid,rscat)
        x,nlist = radial_profile_bin(rsid,rscat,Mbins,usetotnorm=True)
        for i,n in enumerate(nlist):
            ax.plot(x,n,label=labellist[i],drawstyle='steps-mid')
        ax.plot(xtot,ntot,label='Total',drawstyle='steps-mid',color='k',lw=2)
    ax.legend(loc='best',fontsize=10)
    ax.set_xlabel(r'$r/r_{vir}$')
    ax.set_ylabel(r'$n_M(r)/<n(r)>$')

    # Normalizing by <n_M>
    fig,ax=plt.subplots()
    labellist = [r"$10^{%.1f}-10^{%.1f} M_\odot$" % (np.log10(M1),np.log10(M2)) for M1,M2 in zip(Mbins[:-1],Mbins[1:])]
    for haloid,rsid in zip(haloids,rsids):
        hpath = haloutils.get_hpath(haloid,ictype,lx,nv)
        rscat = haloutils.load_rscat(hpath,255)
        xtot,ntot = radial_profile(rsid,rscat)
        x,nlist = radial_profile_bin(rsid,rscat,Mbins,usetotnorm=False)
        for i,n in enumerate(nlist):
            ax.plot(x,n,label=labellist[i],drawstyle='steps-mid')
        ax.plot(xtot,ntot,label='Total',drawstyle='steps-mid',color='k',lw=2)
    ax.legend(loc='best',fontsize=10)
    ax.set_xlabel(r'$r/r_{vir}$')
    ax.set_ylabel(r'$n_M(r)/<n_M(r)>$')
    plt.show()
