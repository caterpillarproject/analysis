import numpy as np
import pylab as plt
import asciitable
import os

from profilefit import NFWprofile,fitNFW
from haloutils import get_foldername

from optparse import OptionParser

def plot_profile(ax,outpath,
                 subfind=False,subfindradius=False,rockstar_get_all=False,
                 plotp03=True,plotrvir=False,plotNFW=False,**kwargs):
    # First line is p03r(Mpc), rvir (kpc)
    # Subsequent lines are r(Mpc), rho (10^10 Msun/Mpc^3)
    if subfind:
        filename = outpath+'/subf-halo-profile.dat'
    elif subfindradius:
        filename = outpath+'/subf-halo-profile-radius.dat'
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
    ax.set_xlabel(r'r [$h^{-1}$ kpc]')
    ax.set_ylabel(r'$r^2 \rho(r)$ [$10^{10} M_\odot$ Mpc$^{-1}$]')
    ax.set_xlim([10**-2,10**3])
    ax.set_ylim([ymin,ymax])
    if plotNFW:
        return rs
    else:
        return None

def find_profilefiles(haloid,subfind,subfindradius,rockstar_get_all,
                  levellist=[11,12,13,14],nrvirlist=[4],
                  basepath="/bigbang/data/AnnaGroup/caterpillar/halos"):
    if subfind:
        profilefile = 'subf-halo-profile.dat'
    elif subfindradius:
        profilefile = 'subf-halo-profile-radius.dat'
    elif rockstar_get_all:
        profilefile = 'rs-halo-profile-allpart.dat'
    else:
        profilefile = 'rs-halo-profile.dat'

    filelist = []
    hidstr = 'H'+str(haloid)
    for lx in levellist:
        for nv in nrvirlist:
            path = basepath+'/'+hidstr+'/'+hidstr+"_BB_Z127_P7_LN7_LX"+str(lx)+"_O4_NV"+str(nv)
            filename = path+'/'+profilefile
            if os.path.exists(filename): filelist.append((lx,path))
    return filelist

if __name__=="__main__":
    PLOTP03 = True
    PLOTRVIR= True
    PLOTNFW = True

    parser = OptionParser()
    parser.add_option("--subfind",action="store_true",dest="subfind",default=False,
                      help="Plot subfind")
    parser.add_option("--subfindradius",action="store_true",dest="subfindradius",default=False,
                      help="Plot subfind (all within radius): TODO (not implemented yet)")
    parser.add_option("--rockstarall",action="store_true",dest="rockstar_get_all",
                      default=False,
                      help="Plot rockstar with subhalo particles")
    parser.add_option("--save",action="store_true",dest="save",default=False,
                      help="Save a png")
    options,args = parser.parse_args()
    sheet = int(args[0])
    print "Plotting sheet %i" % (sheet)

    subfind = options.subfind
    subfindradius = options.subfindradius
    rockstar_get_all = options.rockstar_get_all
    if sum([subfind,subfindradius,rockstar_get_all]) > 1:
        exit("ERROR: don't specify more than one of subfind, subfindradius, rockstar_get_all")

    #fig = plt.figure(figsize=(15,15))
    fig, axes = plt.subplots(3,3,sharex=True,sharey=True,figsize=(15,15),
                             subplot_kw=dict(xscale='log',yscale='log')) 
    plt.subplots_adjust(wspace=0,hspace=0)

    if sheet==1:
        haloidlist = [1194083,1292049,1327707,
                      1232333,1725139,230667,
                      4847,   649524, 706754]
    elif sheet==2:
        haloidlist = [1269360, 134911, 1506656,
                      1665066, 889027, 263605, 
                      1079498, 1326950, 1542569]
    elif sheet==3:
        haloidlist = [1129405, 917969, 299792,
                      299792, 299792, 299792,
                      299792, 299792, 299792]
    elif sheet==10:
        haloidlist = [1764135, 1353966, 889079,
                      795050,  794721, 327580,
                      1476079,135990,795187]
    else:
        exit("Invalid sheet number")
    assert len(haloidlist) == 9
    
    colordict = {11:'blue',12:'red',13:'green',14:'cyan'}
    for i,haloid in enumerate(haloidlist):
        filelist = find_profilefiles(haloid,subfind,subfindradius,rockstar_get_all)
        if len(filelist)==0: 
            print 'ERROR: H%i does not have any profiles (skipping)' % (haloid)
            continue
        else:
            print 'Plotting H%i' % (haloid)
        
        #ax = plt.subplot(3,3,i+1)
        ax = axes[i/3,i%3] 
        for lx,path in filelist:
            rs = plot_profile(ax,path,
                              color=colordict[lx],
                              plotp03=PLOTP03,plotrvir=PLOTRVIR,plotNFW=PLOTNFW,
                              subfind=subfind,subfindradius=subfindradius,
                              rockstar_get_all=rockstar_get_all)
            if rs != None:
                yexponent = 1.3 + 0.2*(14-lx)
                ax.text(10**-1.8,10**yexponent,r"LX%i $r_s=$%3.2f kpc" % (lx,rs),
                        color=colordict[lx],fontsize='x-small')
        plotlabel = 'H'+str(haloid)
        if i==0:
            if subfind:
                plotlabel += ' subf'
            elif subfindradius:
                plotlabel += ' subfradius'
            elif rockstar_get_all:
                plotlabel += ' rockstarsubs'
            else:
                plotlabel += ' rockstar'
        ax.text(10**-1.8,10**2.2,plotlabel,color='black',fontsize='small')
    for ax in np.ravel(axes[0:2,:]):
        ax.set_xlabel('')
    for ax in np.ravel(axes[:,1:3]):
        ax.set_ylabel('')
        
        #plt.title(haloid)
        #plt.loglog()

    #plt.tight_layout()
    if options.save:
        if subfind:
            extra = 'subfind'
        elif subfindradius:
            extra = 'subfindradius'
        elif rockstar_get_all:
            extra = 'rockstarsubs'
        else:
            extra = 'rockstar'
        plt.savefig('cp_profiles_s'+str(sheet)+'_'+extra+'.png',bbox_inches='tight')
    else:
        plt.show()
