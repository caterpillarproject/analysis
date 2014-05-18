import numpy as np
import pylab as plt
import asciitable

from profilefit import NFWprofile,fitNFW

from optparse import OptionParser

def plot_profile(ax,outpath,
                 subfind=False,rockstar_get_all=False,
                 plotp03=True,plotrvir=False,plotNFW=False,**kwargs):
    # First line is p03r(Mpc), rvir (kpc)
    # Subsequent lines are r(Mpc), rho (10^10 Msun/Mpc^3)
    if subfind:
        filename = outpath+'/subf-halo-profile.dat'
    elif rockstar_get_all:
        filename = outpath+'/rs-halo-profile-allpart.dat'
    else:
        filename = outpath+'/rs-halo-profile.dat'
    data = np.array(asciitable.read(filename,delimiter=" ",data_start=1))
    r = 1000.*data['col1']; rho = data['col2']
    f = open(filename,'r')
    p03r,rvir,r200c,halomass = f.readline().split(" ")
    p03r = 1000.*float(p03r); rvir = float(rvir); r200c = float(r200c)
    halomass = float(halomass)
    f.close()

    ax.plot(r,(r/1000.)**2 * rho,**kwargs)
    ymin = 10**-1.5; ymax = 10**2.5
    if plotp03:
        ax.plot([p03r,p03r],[ymin,ymax],ls='--',**kwargs)
    if plotrvir:
        ax.plot([rvir,rvir],[ymin,ymax],'k-.')
        ax.plot([r200c,r200c],[ymin,ymax],'k:')
    if plotNFW:
        good_r = (r > p03r) & (r < rvir)
        rs,rhos = fitNFW(r[good_r],rho[good_r],x0=[50,5],bounds=[(1,3000),(4,8)])
        ax.plot(r,(r/1000.)**2 * NFWprofile(r,rs,rhos),ls=':',lw='0.5',**kwargs)
    plt.xlabel(r'r [$h^{-1}$ kpc]')
    plt.ylabel(r'$r^2 \rho(r)$ [$10^{10} M_\odot$ Mpc$^{-1}$]')
    plt.xlim([10**-2,10**3])
    plt.ylim([ymin,ymax])
    if plotNFW:
        return rs
    else:
        return None

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("--subfind",action="store_true",dest="subfind",default=False,
                      help="Plot subfind")
    parser.add_option("--rockstarall",action="store_true",dest="rockstar_get_all",
                      default=False,
                      help="Plot rockstar with subhalo particles")
    parser.add_option("--save",action="store_true",dest="save",default=False,
                      help="Save a png")
    options,args = parser.parse_args()
    subfind = options.subfind
    rockstar_get_all = options.rockstar_get_all
    if subfind and rockstar_get_all:
        exit("ERROR: don't specify both subfind and rockstar_get_all, only do one")

    fig = plt.figure(figsize=(10,10))
    plt.subplots_adjust(wspace=.3,hspace=.3)
    haloidlist = ['H268422','H121869','H241932','H21047']
    #haloidlist = ['H268422']
    base = "/bigbang/data/AnnaGroup/caterpillar/halos"
    for i,haloid in enumerate(haloidlist):
        ax = plt.subplot(2,2,i+1)
        rs11 = plot_profile(ax,base+'/'+haloid+'/'+haloid+"_BB_Z127_P7_LN7_LX11_O4_NV3",
                            color='blue',plotp03=True,plotrvir=True,plotNFW=True,
                            subfind=subfind,rockstar_get_all=rockstar_get_all)
        if rs11 != None:
            ax.text(10**-1.8,10**1.9,r"LX11 $r_s=$%3.2f kpc" % (rs11),color='blue',fontsize='x-small')
        rs12 = plot_profile(ax,base+'/'+haloid+'/'+haloid+"_BB_Z127_P7_LN7_LX12_O4_NV3",
                            color='red',plotp03=True,plotrvir=True,plotNFW=True,
                            subfind=subfind,rockstar_get_all=rockstar_get_all)
        if rs12 != None:
            ax.text(10**-1.8,10**1.7,r"LX12 $r_s=$%3.2f kpc" % (rs12),color='red', fontsize='x-small')
        if haloid=='H241932' and subfind:
            rs13 = plot_profile(ax,base+'/'+haloid+'/'+haloid+"_BB_Z127_P7_LN7_LX13_O4_NV3",
                                color='green',plotp03=True,plotrvir=True,plotNFW=True,
                                subfind=subfind,rockstar_get_all=rockstar_get_all)
            if rs13 != None:
                ax.text(10**-1.8,10**1.5,r"LX13 $r_s=$%3.2f kpc" % (rs13),color='green', fontsize='x-small')
        if haloid=='H268422' and subfind:
            rs13 = plot_profile(ax,base+'/'+haloid+'/'+haloid+"_BB_Z127_P7_LN7_LX13_O4_NV3",
                                color='green',plotp03=True,plotrvir=True,plotNFW=True,
                                subfind=subfind,rockstar_get_all=rockstar_get_all)
            if rs13 != None:
                ax.text(10**-1.8,10**1.5,r"LX13 $r_s=$%3.2f kpc" % (rs13),color='green', fontsize='x-small')
            rs14 = plot_profile(ax,base+'/'+haloid+'/'+haloid+"_BB_Z127_P7_LN7_LX14_O4_NV3",
                                color='cyan',plotp03=True,plotrvir=True,plotNFW=True,
                                subfind=subfind,rockstar_get_all=rockstar_get_all)
            if rs14 != None:
                ax.text(10**-1.8,10**1.3,r"LX14 $r_s=$%3.2f kpc" % (rs14),color='cyan', fontsize='x-small')

        if i==0:
            if subfind:
                plotlabel = 'Subfind'
            elif rockstar_get_all:
                plotlabel = 'Rockstar (subs)'
            else:
                plotlabel = 'Rockstar'
            ax.text(10**-1.8,10**2.2,plotlabel,color='black',fontsize='small')
        
        plt.title(haloid)
        plt.loglog()
    if options.save:
        if subfind:
            extra = 'subfind'
        elif rockstar_get_all:
            extra = 'rockstarsubs'
        else:
            extra = 'rockstar'
        plt.savefig('four_profiles_'+extra+'.png',bbox_tight=True)
    else:
        plt.show()
