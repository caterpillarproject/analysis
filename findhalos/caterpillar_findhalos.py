import sys
import numpy as np
import asciitable
import glob
from optparse import OptionParser
from operator import itemgetter

import haloutils

def get_zoom_id(parenthid,hcat,scat,pcat,
                cutoffratio=0.07,checknum=10,maxdist=.1,
                hubble = 0.6711,
                nofix=False,verbose=False,retflag=False):
    hosts = hcat.get_hosts()
    mostpartid = hosts.index[hosts['npart'].argmax()]
    mostpart = np.max(hosts['npart'])
    mostpos = hosts.ix[mostpartid][['posX','posY','posZ']]
    badsubfflag = False

    bestgroup = min(enumerate(scat.group_contamination_count[0:checknum]),key=itemgetter(1))[0]
    if bestgroup != 0: 
        print "     WARNING: subfind bestgroup != 0 (is %i)" % (bestgroup)
        badsubfflag = True
        #print "group masses:",scat.group_mass[0:(bestgroup+1)]*10**10/header.hubble
        #print "group npart:",scat.group_len[0:(bestgroup+1)]
        #print "group contamination:",scat.group_contamination_count[0:(bestgroup+1)]
    #substart = scat.group_firstsub[bestgroup]
    #subnum = scat.group_nsubs[bestgroup]
    #print "     R200c mass of bestgroup %i: %3.2e" % (bestgroup, scat.group_m_crit200[bestgroup] * 10**10/hubble)
    scatpos = scat.group_pos[bestgroup]
    distdiff = np.sqrt(np.sum((mostpos-scatpos)**2))
    if distdiff > maxdist:
        print "     WARNING: mostpart and bestgroup pos differ by %3.2f" % (distdiff)
        badsubfflag = True

    parentmass = pcat.ix[parenthid]['mvir']
    zoommass = hosts.ix[mostpartid]['mvir']
    if (np.abs(zoommass/parentmass - 1) <= cutoffratio) or options.nofix:
        if retflag: return mostpartid,False,badsubfflag
        return mostpartid #this one is probably the right one
    #otherwise, check the top several and use the subhalo particles too

    print "     fixing %i... mostpart (id %i): %i ratiodiff %3.2f (using subhalo parts too)" % (parenthid,mostpartid,mostpart,np.abs(zoommass/parentmass - 1))
    idlist = []; npartlist = []
    while len(idlist) <= 0:
        largestids = hosts.index[np.argsort(np.array(hosts['npart']))[-checknum-1:-1]] #sorts ascending
        for thisid in largestids:
            thisratio = np.abs(hosts.ix[thisid]['mvir']/parentmass-1)
            if thisratio <= cutoffratio:
                thispart = hcat.get_all_num_particles_from_halo(thisid) #np.array(hosts.ix[thisid]['npart'])
                thispos  = np.array(hosts.ix[thisid][['posX','posY','posZ']])
                thisdist = np.sqrt(np.sum((mostpos-thispos)**2)) * 1000.
                groupdist = np.sqrt(np.sum((thispos-scatpos)**2)) * 1000.
                print "     id %8i, ratiodiff %3.2f, partnum: %8i, distance: %7.2f kpc/h groupdist: %7.2f kpc/h" % (thisid,thisratio,thispart,thisdist,groupdist) #float(mostpart-thispart)/float(mostpart)
                idlist.append(thisid); npartlist.append(thispart)
        if len(idlist)==0:
            print "     didn't find matches in the first %i at cutoff %3.2f, increasing both" % (checknum, cutoffratio)
            checknum = checknum + 5; cutoffratio = cutoffratio + 0.01
    idlist = np.array(idlist); npartlist = np.array(npartlist)
    if retflag: return idlist[np.argmax(npartlist)],True,badsubfflag
    return idlist[np.argmax(npartlist)]

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("-w","--write-summary",dest="write_summary",
                      action="store_true",default=False,
                      help="flag to write a summary file in /bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index.txt (using asciitable FixedWidth format)")
    parser.add_option("-n","--no-fix",dest="nofix",
                      action="store_true",default=False,
                      help="flag to use simple zoomid (summary file in /bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index_nofix.txt)")
    parser.add_option("-o","--include-old",dest="includeold",
                      action="store_true",default=False,
                      help="flag to include oldhalos")
    (options, args) = parser.parse_args()

    levellist = [11,12,13,14]
    nrvirlist = [3,4]

    print "Reading parent..."
    pcat = haloutils.load_pcatz0()
    hlist = haloutils.find_halo_paths(require_rockstar=True,require_subfind=True,
                                      require_sorted=True,
                                      nrvirlist=nrvirlist,levellist=levellist)
    if options.includeold:
        pcatold = haloutils.load_pcatz0(old=True)
        hlist = hlist+haloutils.find_halo_paths(require_rockstar=True,require_subfind=True,
                                                nrvirlist=nrvirlist,levellist=levellist,
                                                basepath="/bigbang/data/AnnaGroup/caterpillar/halos/oldhalos",hdf5=False)
    print "Number of halos with rockstar and subfind:",len(hlist)

    hindex = []
    
    for hpath in hlist:
        oldflag = 'oldhalos' in hpath
        parenthid = haloutils.get_parent_hid(hpath)
        ictype,lx,nv = haloutils.get_zoom_params(hpath)
        lastsnap = haloutils.get_numsnaps(hpath) - 1
        hcat = haloutils.load_rscat(hpath,lastsnap)
        scat = haloutils.load_scat(hpath)
        if oldflag:
            print "--old halo"
            zoomid,badhaloflag,badsubfflag = get_zoom_id(parenthid,hcat,scat,pcatold,verbose=True,retflag=True,nofix=options.nofix)
        else:
            zoomid,badhaloflag,badsubfflag = get_zoom_id(parenthid,hcat,scat,pcat,verbose=True,retflag=True,nofix=options.nofix)
        zoommass = hcat.ix[zoomid]['mvir']/hcat.h0
        zoomrvir = hcat.ix[zoomid]['rvir']/hcat.h0
        zoomx,zoomy,zoomz = hcat.ix[zoomid][['posX','posY','posZ']]
        print "Halo %7i LX %2i has ID %6i: M = %3.2e R = %4.1f (%4.2f,%4.2f,%4.2f)" % (parenthid, lx, zoomid, zoommass,zoomrvir, zoomx,zoomy,zoomz)
        hindex.append([parenthid,ictype,lx,nv,zoomid,
                       zoomx,zoomy,zoomz,zoommass,zoomrvir,int(badhaloflag),int(oldflag),int(badsubfflag)])
        sys.stdout.flush()
        
    if options.write_summary:
        if options.nofix:
            outname = "/bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index_nofix.txt"
        else:
            outname = "/bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index.txt"
        asciitable.write(hindex,outname,
                         Writer=asciitable.FixedWidth,
                         names=['parentid','ictype','LX','NV','zoomid','x','y','z','mvir','rvir','badflag','oldflag','badsubf'],
                         formats={'x': '%0.3f','y': '%0.3f','z': '%0.3f',
                                  'mvir':'%3.2e','rvir': '%0.1f'})
