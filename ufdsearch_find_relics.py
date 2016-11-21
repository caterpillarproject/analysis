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

for fnames in glob.glob("UFDSEARCHTMP/*.p"):
    hid = os.path.basename(fnames)[:-2]
    hid = haloutils.hidint(hid)
    hpath = haloutils.get_hpath_lx(hid,14)
    print hid, fnames
    
    rscat = haloutils.load_rscat(hpath,haloutils.get_lastsnap(hpath))
    zoomid = haloutils.load_zoomid(hpath)
    host = rscat.ix[zoomid]
    subs = get_all_subhalos_within_halo(zoomid)
    subids = set(subs['id'])
    

    print "Loading MTC"
    start = time.time()
    mtc = haloutils.load_mtc(hpath,indexbyrsid=True)
    print "{:.2f}".format(time.time()-start)

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
    
    # Get Greg's extant data
    dataE = AE.read(hpath)
    dataE.index = dataE['rsid']
    
    greg_rsids = set(dataE.index)
    if not greg_rsids.issuperset(subids_of_relics_in_host):
        print "ERROR: not all relic subs in Greg's extant catalog. Here are the bad IDs:"
        print set(subids_of_relics_in_host).difference(greg_rsids)
        continue 

    ## TODO get MTaddition for reionization quantities

    ## TODO Apply selection to find surviving satellites from Greg's criteria

    ## TODO Properties of Greg's satellites we did not get in our selection

    ## TODO Properties of our satellites not in Greg's catalog
    
    break
