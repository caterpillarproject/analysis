import glob
import haloutils
import os, sys, time
import cPickle as pickle
from six import iteritems
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd

import MTaddition as mtadd
import MTanalysis3 as mta
AE = mta.AllExtantData()
E = mtadd.ExtantDataReionization()

## Greg's abundance matching
sys.path.insert(0,"/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code")
import abundance_matching as greg_am
## Added infall_times.py to my local directory

#logMbins = np.arange(6,12.1,.1)
#Htot = np.zeros(len(logMbins)-1)
numgreg = []
numrelics = []

#suffix = ""
#suffix = "_7585"
suffix = "_7685"
#suffix = "_7785"
#suffix = "_7885"
#suffix = "_vmax12gd"

## TODO gotta do a better way for grabbing the relevant files
for fnames in glob.glob("UFDSEARCHTMP/*"+suffix+".p"):
    if "MBs" in fnames: continue
    if suffix=="":
        if "_7585" in fnames or "_7885" in fnames or "_vmax12gd" in fnames: continue
        hid = os.path.basename(fnames)[:-2]
        _suffix = "_8085"
    elif suffix=="_vmax12gd":
        hid = os.path.basename(fnames)[:-11]
        _suffix = suffix
    else:
        hid = os.path.basename(fnames)[:-7]
        _suffix = suffix

    hid = haloutils.hidint(hid)
    hpath = haloutils.get_hpath_lx(hid,14)
    print hid, fnames
    sys.stdout.flush()
    
    rscat = haloutils.load_rscat(hpath,haloutils.get_lastsnap(hpath))
    zoomid = haloutils.load_zoomid(hpath)
    host = rscat.ix[zoomid]
    subs = rscat.get_all_subhalos_within_halo(zoomid)
    subids = set(subs['id'])
    

    print "Loading MTC"
    start = time.time()
    mtc = haloutils.load_mtc(hpath,indexbyrsid=True)
    print "{:.2f}".format(time.time()-start)
    sys.stdout.flush()

    start = time.time()
    with open(fnames,"r") as f:
        data = pickle.load(f)
    relic_indices, relic_rows, merged_relic_flags, merged_relic_flags2 = data

    mtkeys = merged_relic_flags.keys()
    in_main_host = [mtkey in subids for mtkey in mtkeys]
    in_main_host = dict(zip(mtkeys, in_main_host))

    # Get out relics that did not merge
    subids_of_relics = []
    subids_of_relics_in_host = []
    num_total_relics = 0
    relics = {}
    for mtkey, flags in iteritems(merged_relic_flags):
        mt = mtc[mtkey]
        relic_cands = mt[relic_rows[mtkey]]
        if len(relic_cands[~flags]) > 0:
            relics[mtkey] = relic_cands[~flags]
            num_total_relics += len(relics[mtkey])
            subids_of_relics.append(mtkey)
            if in_main_host[mtkey]: subids_of_relics_in_host.append(mtkey)
    # Pull out z=0 relics
    #rscat_relics = rscat.ix[subids_of_relics]
    # Filter by those in the host
    sub_relics = rscat.ix[subids_of_relics_in_host]
    num_relics_in_host = len(sub_relics)
    
    print "{}/{} relics are in main host".format(num_relics_in_host, num_total_relics)
    sys.stdout.flush()
    
    ### Get Greg's extant data (based on DwarfMethods.py get_extant_data())
    dataE = AE.read(hpath)
    data_reion = E.read(hpath)
    extant = pd.concat((dataE, data_reion), axis=1)
    ## Somehow this cut misses some halos in the virial radius of the host
    #pidmask = dataE['pid'] == zoomid
    #extant = pd.concat((dataE[pidmask], data_reion[pidmask]), axis=1)
    extant.index = extant['rsid']
    ## Cut depth==0 which is only the actually extant objects
    ## depth is recursively searching for things that have merged
    extant = extant[extant['depth']==0]

    greg_rsids = set(extant.index)
    if not greg_rsids.issuperset(subids_of_relics_in_host):
        print "ERROR: not all relic subs in Greg's extant catalog. Here are the bad IDs:"
        badIDs = set(subids_of_relics_in_host).difference(greg_rsids)
        print badIDs
        #continue 
        print "Continuing anyway removing badIDs from my relics"
        for badID in badIDs:
            sub_relics = sub_relics[sub_relics.index != badID]

    ##### Apply selection to find surviving satellites from Greg's criteria
    ## These criteria are chosen to reproduce the statistics of the 
    ## Barber et al. 2014 luminous fraction (Fig 2 of Greg's paper)
    vmax_ach = 9.48535156
    vmax_filt= 23.54248047
    mask_pre = np.array(extant['vmax_12'] > vmax_ach)
    mask_post= np.array(extant['peak_vmax'] >= vmax_filt)
    
    ## Use abundance matching to grab UFDs
    ## 1000 to 2e5 Msun
    ## logM ~ 8.2 to 9.4
    am_model = greg_am.GarrisonKimmel()
    mpeak_upper = am_model.stellar_to_halo_mass(2e5)
    mpeak_lower = am_model.stellar_to_halo_mass(1e3)
    
    ## This is Greg's satellites
    in_range = mask_pre | mask_post
    ii_ufd = np.logical_and(np.array(extant[in_range]['max_mass'] < mpeak_upper),
                            np.array(extant[in_range]['max_mass'] > mpeak_lower))
    ufds = extant[in_range][ii_ufd]
    
    numgreg.append(len(ufds))

    ## This is my satellites
    sub_relics_max_mass = extant.ix[sub_relics.index]['max_mass']
    in_range = np.logical_and(sub_relics_max_mass < mpeak_upper,
                              sub_relics_max_mass > mpeak_lower)
    sub_relics = sub_relics[in_range]
    numrelics.append(len(sub_relics))
    
    print "We have {} vs Greg has {}".format(numrelics[-1],numgreg[-1])
    sys.stdout.flush()

    #h,x = np.histogram(np.log10(ufds['m200_8']),bins=logMbins)
    #Htot += h
    
    ## TODO Properties of our satellites not in Greg's catalog if any
    ## There are none in my test example so just count how many there are for now

    ## Save main branches
    aj_mbs = []
    gd_mbs = []
    print "    ----Going through AJ----"
    for subid in sub_relics.index:
        try:
            mt = mtc[subid]
        except:
            print "    {} not in mtc".format(subid)
        else:
            if mt[0]['origid'] != subid:
                print "    {} does not match".format(subid)
            else:
                aj_mbs.append(mt.getMainBranch(0))
    
    print "    ----Going through GD----"
    for subid in ufds.index:
        try:
            mt = mtc[subid]
        except:
            print "    {} not in mtc".format(subid)
        else:
            if mt[0]['origid'] != subid:
                print "    {} does not match".format(subid)
            else:
                gd_mbs.append(mt.getMainBranch(0))

    with open("UFDSEARCHTMP/H{}_MBs{}.p".format(hid,_suffix),"w") as f:
        pickle.dump([aj_mbs,gd_mbs],f)

    ## End Loop
    sys.stdout.flush()
    
#print Htot
print numrelics
print numgreg
print np.array(numrelics).astype(float)/np.array(numgreg).astype(float)

