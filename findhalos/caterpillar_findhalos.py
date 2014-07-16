import sys
import numpy as np
import asciitable
import glob
from optparse import OptionParser

import haloutils

def get_zoom_id(parenthid,hcat,pcat,cutoffratio=0.07,checknum=10,verbose=False,retflag=False):
    hosts = hcat.get_hosts()
    mostpartid = hosts.index[hosts['npart'].argmax()]
    mostpart = np.max(hosts['npart'])
    mostpos = hosts.ix[mostpartid][['posX','posY','posZ']]

    parentmass = pcat.ix[parenthid]['mvir']
    zoommass = hosts.ix[mostpartid]['mvir']
    if np.abs(zoommass/parentmass - 1) <= cutoffratio:
        if retflag: return mostpartid,False
        return mostpartid #this one is probably the right one
    #otherwise, check the top several and use the subhalo particles too

    print "     fixing %i... mostpart: %i ratiodiff %3.2f (using subhalo parts too)" % (parenthid,mostpart,np.abs(zoommass/parentmass - 1))
    idlist = []; npartlist = []
    while len(idlist) <= 0:
        largestids = hosts.index[np.argsort(np.array(hosts['npart']))[-checknum-1:-1]] #sorts ascending
        for thisid in largestids:
            thisratio = np.abs(hosts.ix[thisid]['mvir']/parentmass-1)
            if thisratio <= cutoffratio:
                thispart = hcat.get_all_num_particles_from_halo(thisid) #np.array(hosts.ix[thisid]['npart'])
                thispos  = np.array(hosts.ix[thisid][['posX','posY','posZ']])
                thisdist = np.sqrt(np.sum((mostpos-thispos)**2)) * 1000.
                print "     id %8i, ratiodiff %3.2f, partnum: %8i, distance: %7.2f kpc/h" % (thisid,thisratio,thispart,thisdist) #float(mostpart-thispart)/float(mostpart)
                idlist.append(thisid); npartlist.append(thispart)
        if len(idlist)==0:
            print "     didn't find matches in the first %i at cutoff %3.2f, increasing both" % (checknum, cutoffratio)
            checknum = checknum + 5; cutoffratio = cutoffratio + 0.01
    idlist = np.array(idlist); npartlist = np.array(npartlist)
    if retflag: return idlist[np.argmax(npartlist)],True
    return idlist[np.argmax(npartlist)]

    #zoommassratios = hcat['mvir']/parentmass
    #mostpartids = 

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("-w","--write-summary",dest="write_summary",
                      action="store_true",default=False,
                      help="flag to write a summary file in /bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index.txt (using asciitable FixedWidth format)")
    (options, args) = parser.parse_args()

    levellist = [11,12,13,14]
    nrvirlist = [3,4]

    print "Reading parent..."
    pcat = haloutils.load_pcatz0()
    pcatold = haloutils.load_pcatz0(old=True)
    hlist = haloutils.find_halo_paths(require_rockstar=True,
                                      nrvirlist=nrvirlist,levellist=levellist)
    hlist = hlist+haloutils.find_halo_paths(require_rockstar=True,
                                            nrvirlist=nrvirlist,levellist=levellist,
                                            basepath="/bigbang/data/AnnaGroup/caterpillar/halos/oldhalos",hdf5=False)
    print "Number of halos with rockstar:",len(hlist)

    hindex = []
    
    for hpath in hlist:
        parenthid = haloutils.get_parent_hid(hpath)
        ictype,lx,nv = haloutils.get_zoom_params(hpath)
        lastsnap = haloutils.get_numsnaps(hpath) - 1
        hcat = haloutils.load_rscat(hpath,lastsnap)
        if 'oldhalos' in hpath:
            zoomid,badhaloflag = get_zoom_id(parenthid,hcat,pcatold,verbose=True,retflag=True)
        else:
            zoomid,badhaloflag = get_zoom_id(parenthid,hcat,pcat,verbose=True,retflag=True)
        zoommass = hcat.ix[zoomid]['mvir']/hcat.h0
        zoomrvir = hcat.ix[zoomid]['rvir']/hcat.h0
        zoomx,zoomy,zoomz = hcat.ix[zoomid][['posX','posY','posZ']]
        print "Halo %7i LX %2i has ID %6i: M = %3.2e R = %4.1f (%4.2f,%4.2f,%4.2f)" % (parenthid, lx, zoomid, zoommass,zoomrvir, zoomx,zoomy,zoomz)
        hindex.append([parenthid,ictype,lx,nv,zoomid,
                       zoomx,zoomy,zoomz,zoommass,zoomrvir,int(badhaloflag)])
        sys.stdout.flush()
        
    if options.write_summary:
        asciitable.write(hindex,"/bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index.txt",
                         Writer=asciitable.FixedWidth,
                         names=['parentid','ictype','LX','NV','zoomid','x','y','z','mvir','rvir','badflag'],
                         formats={'x': '%0.3f','y': '%0.3f','z': '%0.3f',
                                  'mvir':'%3.2e','rvir': '%0.1f'})
