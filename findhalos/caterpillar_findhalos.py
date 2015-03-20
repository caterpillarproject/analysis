import sys
import numpy as np
import asciitable
import glob,os
from optparse import OptionParser
from operator import itemgetter

import haloutils

def get_zoom_id(parenthid,hcat,scat,pcat,
                cutoffratio=0.07,checknum=10,maxdist=.1,
                hubble = 0.6711,
                nofix=False,verbose=False,retflag=False):
    hosts = hcat.get_hosts()
    mostpartid = hosts.index[np.array(hosts['num_cp']).argmax()]
    mostpart = np.max(hosts['num_cp'])
    mostpos = hosts.ix[mostpartid][['posX','posY','posZ']]
    badsubfflag = False

    if scat != None:
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
    elif scat == None:
        badsubfflag=True

    parentmass = pcat.ix[parenthid]['mgrav']
    zoommass = hosts.ix[mostpartid]['mgrav']
    if (np.abs(zoommass/parentmass - 1) <= cutoffratio) or nofix:
        if retflag: return mostpartid,False,badsubfflag
        return mostpartid #this one is probably the right one
    #otherwise, check the top several and use the subhalo particles too

    print "     fixing %i... mostpart (id %i): %i ratiodiff %3.2f, parentmass %3.2e" % (parenthid,mostpartid,mostpart,np.abs(zoommass/parentmass - 1),parentmass/hcat.h0)
    idlist = []; npartlist = []
    while len(idlist) <= 0:
        largestids = hosts.index[np.argsort(np.array(hosts['num_cp']))[-checknum:]] #sorts ascending
        for thisid in largestids:
            thisratio = np.abs(hosts.ix[thisid]['mgrav']/parentmass-1)
            #print "         %i %3.2f" % (thisid,thisratio)
            if thisratio <= cutoffratio:
                thispart = hcat.get_all_num_particles_from_halo(thisid) #np.array(hosts.ix[thisid]['npart'])
                thispos  = np.array(hosts.ix[thisid][['posX','posY','posZ']])
                thisdist = np.sqrt(np.sum((mostpos-thispos)**2)) * 1000.
                if scat != None:
                    groupdist = np.sqrt(np.sum((thispos-scatpos)**2)) * 1000.
                    print "     id %8i, ratiodiff %3.2f, partnum: %8i, distance: %7.2f kpc/h groupdist: %7.2f kpc/h" % (thisid,thisratio,thispart,thisdist,groupdist) #float(mostpart-thispart)/float(mostpart)
                else:
                    print "     id %8i, ratiodiff %3.2f, partnum: %8i, distance: %7.2f kpc/h (no groupdist)" % (thisid,thisratio,thispart,thisdist) #float(mostpart-thispart)/float(mostpart)
                idlist.append(thisid); npartlist.append(thispart)
        if len(idlist)==0:
            print "     didn't find matches in the first %i at cutoff %3.2f, increasing both" % (checknum, cutoffratio)
            checknum = checknum + 5; cutoffratio = cutoffratio + 0.01
    idlist = np.array(idlist); npartlist = np.array(npartlist)
    if retflag: return idlist[np.argmax(npartlist)],True,badsubfflag
    return idlist[np.argmax(npartlist)]

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("-f","--force",dest="force",
                      action="store_true",default=False,
                      help="flag to re-find all halos (as opposed to just updating the table)")
    parser.add_option("-w","--write-summary",dest="write_summary",
                      action="store_true",default=False,
                      help="flag to write a summary file in /bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index.txt (using asciitable FixedWidth format)")
    parser.add_option("-n","--no-fix",dest="nofix",
                      action="store_true",default=False,
                      help="flag to use simple zoomid (summary file in /bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index_nofix.txt)")
    parser.add_option("-c","--contam",dest="contam",
                      action="store",type="int",default=0,
                      help="flag to run index on contamination suite: 0 (default) for normal, 1 for low mass, 2 for med mass, 3 for high mass")
    (options, args) = parser.parse_args()

    levellist = [11,12,13,14]
    nrvirlist = [3,4,5,6,7]


    print "Reading parent..."
    pcat = haloutils.load_pcatz0()
    if options.contam != 0:
        if options.contam==1: contampath=haloutils.global_halobase+'/low_mass_halos'
        elif options.contam==2: contampath=haloutils.global_halobase+'/middle_mass_halos'
        elif options.contam==3: contampath=haloutils.global_halobase+'/high_mass_halos'
        else: raise ValueError("contam argument must be 0, 1, 2, or 3")
        hlist = haloutils.find_halo_paths(basepath=contampath,
                                          require_rockstar=True,require_subfind=True,
                                          contamsuite=True,onlychecklastsnap=True,
                                          nrvirlist=nrvirlist,levellist=levellist,
                                          use_fullbin_rockstar=True)
        outname = contampath+"/contam_zoom_index.txt"
    else:
        hlist = haloutils.find_halo_paths(require_rockstar=True,require_subfind=False,
                                          onlychecklastsnap=True,
                                          nrvirlist=nrvirlist,levellist=levellist,verbose=True)
        if options.nofix:
            outname = "/bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index_nofix.txt"
        else:
            outname = "/bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index.txt"
    print "Number of halos with rockstar and subfind:",len(hlist)
    hindex = []
    
    if not options.force: #what halos have already been added to the index?
        currentindextable = asciitable.read(outname,Reader=asciitable.FixedWidth)
        currentindex = []
        existinglines = {}
        for line in currentindextable:
            currentindex.append(list(line))
            key = haloutils.hidstr(line[0])+'_'+line[1]+'_'+str(line[2])+'_'+str(line[3])
            existinglines[key]=line

    for hpath in hlist:
        parenthid = haloutils.get_parent_hid(hpath)
        ictype,lx,nv = haloutils.get_zoom_params(hpath)

        if not options.force: #skip if already in index
            key = haloutils.hidstr(parenthid)+'_'+ictype+'_'+str(lx)+'_'+str(nv)
            if key in existinglines:
                hindex.append(existinglines[key])
                print "Already in index: "+hpath
                continue

        lastsnap = haloutils.get_numsnaps(hpath) - 1
        if options.contam != 0:
            try:
                hcat = haloutils.load_rscat(hpath,lastsnap,halodir='halos',version=8,rmaxcut=False)
            except:
                hcat = haloutils.load_rscat(hpath,lastsnap,halodir='halos',version=7,rmaxcut=False)
        else:
            hcat = haloutils.load_rscat(hpath,lastsnap,rmaxcut=False)

        try:
            scat = haloutils.load_scat(hpath)
        except IOError:
            scat = None
        zoomid,badhaloflag,badsubfflag = get_zoom_id(parenthid,hcat,scat,pcat,verbose=True,retflag=True,nofix=options.nofix)
        zoommass = hcat.ix[zoomid]['mgrav']/hcat.h0
        zoomrvir = hcat.ix[zoomid]['rvir']/hcat.h0
        zoomx,zoomy,zoomz = hcat.ix[zoomid][['posX','posY','posZ']]
        hpos = np.array([zoomx,zoomy,zoomz])
        print hpath
        print "Halo %7i LX %2i has ID %6i: M = %4.3e R = %4.1f (%4.2f,%4.2f,%4.2f)" % (parenthid, lx, zoomid, zoommass,zoomrvir, zoomx,zoomy,zoomz)
        if options.contam != 0:
            icsize = np.sum([os.path.getsize(icfile) for icfile in glob.glob(hpath+'/ics.*')])/10.**6 #MB
            #hpos = np.array([zoomx,zoomy,zoomz])#np.array(hcat.ix[zoomid][['posX','posY','posZ']])
            drarr = []
            for parttype in [2,3,4,5]:
                ppos = haloutils.load_partblock(hpath,255,"POS ",parttype=parttype)
                drarr.append(np.min(np.sqrt(np.sum((ppos-hpos)**2,1))))
            hindex.append([parenthid,ictype,lx,nv,zoomid,
                           icsize,drarr[0],drarr[1],drarr[2],drarr[3],
                           zoomx,zoomy,zoomz,zoommass,zoomrvir,
                           int(badhaloflag),int(badsubfflag)])
        else:
            allsnapsthere = True
            for snap in xrange(256):
                snapstr = str(snap).zfill(3)
                snappath = hpath+"/outputs/snapdir_"+snapstr+"/snap_"+snapstr+".0.hdf5"
                if (not os.path.exists(snappath)):
                    allsnapsthere = False
                    break
            ppos = haloutils.load_partblock(hpath,255,"POS ",parttype=2)
            min2 = np.min(np.sqrt(np.sum((ppos-hpos)**2,1)))
            hindex.append([parenthid,ictype,lx,nv,zoomid,min2,
                           zoomx,zoomy,zoomz,zoommass,zoomrvir,
                           int(badhaloflag),int(badsubfflag),int(allsnapsthere)])
        sys.stdout.flush()

    if options.write_summary:
        if options.contam != 0:
            asciitable.write(hindex,outname,
                             Writer=asciitable.FixedWidth,
                             names=['parentid','ictype','LX','NV','zoomid',
                                    'icsize','min2','min3','min4','min5',
                                    'x','y','z','mgrav','rvir','badflag','badsubf'],
                             formats={'x': '%0.3f','y': '%0.3f','z': '%0.3f',
                                      'mgrav':'%4.3e','rvir': '%0.1f',
                                      'min2':'%0.6f','min3':'%0.6f','min4':'%0.6f','min5':'%0.6f',})
        else:
            asciitable.write(hindex,outname,
                             Writer=asciitable.FixedWidth,
                             names=['parentid','ictype','LX','NV','zoomid','min2',
                                    'x','y','z','mgrav','rvir',
                                    'badflag','badsubf','allsnaps'],
                             formats={'x': '%0.3f','y': '%0.3f','z': '%0.3f',
                                      'mgrav':'%4.3e','rvir': '%0.1f','min2':'%0.6f'})
