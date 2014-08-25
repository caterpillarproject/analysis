# core
import numpy as np
import os, sys, platform, asciitable
from optparse import OptionParser
from operator import itemgetter

# caterpillar
import readhalos.RSDataReader as RDR
import readhalos.readsubf as rsf
import readsnapshots.readsnapHDF5 as rs
import readsnapshots.readsnapHDF5_greg as rsg
from brendanlib.grifflib import determinebasepath

# analysis
import findhalos.haloutils as haloutils

h0 = .6711

def contamdist(outpath,pzindex,checknum=5,verbose=False):
    haloid = haloutils.get_parent_hid(outpath)
    ictype,lx,nv = haloutils.get_zoom_params(outpath)
    #if haloid not in pzindex['parentid']:
    #print str(haloid)+" not in parent zoom index, finding halo with subfind"

    if verbose:
        print "---using subfind for "+haloutils.get_foldername(outpath)
        print "most massive 5 groups:"
    scat = haloutils.load_scat(outpath)
    mindistlist = []
    rvirlist = []
    contamlist = []
    for i in range(checknum):
        gpos = scat.group_pos[i]
        mindist = 1e5
        for parttype in [2,3,4,5]:
            ppos = rsg.read_block(outpath+'/outputs/snapdir_255/snap_255',"POS ",parttype=2)
            dist = np.sqrt(np.sum((ppos-gpos)**2,axis=1))
            thismindist = np.min(dist)*1000./h0
            if thismindist < mindist:
                mindist = thismindist
                if parttype != 2: print "closest contam is parttype "+str(parttype)
        mindistlist.append(mindist)
        rvir = scat.group_r_crit200[i]*1000./h0
        rvirlist.append(rvir)
        contamlist.append(scat.group_contamination_count[i])
        if verbose:
            print "  mass     rvir  rcontam ccount cmass partdiff"
            print "  %3.2e %4.1f %7.1f %6i %3.2e" % \
                (scat.group_m_crit200[i]*10**10/h0,
                 rvir,
                 mindist,
                 scat.group_contamination_count[i],
                 scat.group_contamination_mass[i]*10**10/h0,
                 scat.group_len)

    bestgroup = min(enumerate(scat.group_contamination_count[0:checknum]),key=itemgetter(1))[0]
    if verbose: print "bestgroup: %i" % (bestgroup)
    #spos = scat.group_pos[bestgroup]
    return bestgroup,rvirlist[bestgroup],contamlist[bestgroup],mindistlist[bestgroup]

if __name__=="__main__":
    basepath = determinebasepath(platform.node())
    halobase = basepath+'/caterpillar/halos'
    pzindex = haloutils.get_parent_zoom_index()

    parser = OptionParser()
    options,args = parser.parse_args()

    
    haloids = ['H1194625']
    ictypes = ['BC','BE','BB']
    lxlist = [11]
    nvlist = [1,2,3,4,5,6]

    distarr = np.zeros((len(haloids),len(ictypes),len(lxlist),len(nvlist)))-1
    rvirarr = np.zeros((len(haloids),len(ictypes),len(lxlist),len(nvlist)))-1
    bestgrouparr = np.zeros((len(haloids),len(ictypes),len(lxlist),len(nvlist)))-1
    contamarr = np.zeros((len(haloids),len(ictypes),len(lxlist),len(nvlist)))-1

    for i0,haloid in enumerate(haloids):
        for i1,ictype in enumerate(ictypes):
            for i2,lx in enumerate(lxlist):
                for i3,nv in enumerate(nvlist):
                    outpath = haloutils.get_outpath(haloid,ictype,lx,nv,halobase)
                    if os.path.exists(outpath+'/outputs/snapdir_255/snap_255.0.hdf5'):
                        bestgroup,rvir,contamcount,dist = contamdist(outpath,pzindex)
                        distarr[i0,i1,i2,i3] = dist
                        rvirarr[i0,i1,i2,i3] = rvir
                        contamarr[i0,i1,i2,i3] = contamcount
                        bestgrouparr[i0,i1,i2,i3] = bestgroup
    distarr[distarr==-1] = np.nan
    rvirarr[rvirarr==-1] = np.nan
    bestgrouparr[bestgrouparr==-1] = np.nan
    contamarr[contamarr==-1] = np.nan
    print bestgrouparr[0,0,0,:]
    print bestgrouparr[0,1,0,:]
    import pylab as plt
    plt.figure()
    plt.plot(nvlist,distarr[0,0,0,:],'bo-',label='BC')
    plt.plot(nvlist,distarr[0,1,0,:],'go-',label='BE')
    plt.legend(loc='best',fontsize=18)
    for i,nv in enumerate(nvlist):
        rvir = rvirarr[0,0,0,i]
        plt.plot([i+1-.25,i+1+.25],[rvir,rvir],'b--')
    for i,nv in enumerate(nvlist):
        rvir = rvirarr[0,1,0,i]
        plt.plot([i+1-.25,i+1+.25],[rvir,rvir],'g--')
    plt.xlabel('NV')
    plt.ylabel('Closest Contamination Particle (kpc)')
    plt.title(haloids[0])
    plt.savefig(str(haloids[0])+'_contamdist.png',bbox_inches='tight')

    plt.figure()
    plt.plot(nvlist,contamarr[0,0,0,:],'bo-',label='BC')
    plt.plot(nvlist,contamarr[0,1,0,:],'go-',label='BE')
    plt.xlabel('NV')
    plt.ylabel('Contamination Count')
    plt.title(haloids[0])
    plt.savefig(str(haloids[0])+'_contamcount.png',bbox_inches='tight')

    #plt.savefig()
