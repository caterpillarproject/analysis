import glob
import haloutils
import os, sys, time
import cPickle as pickle
from six import iteritems
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd

## Greg's abundance matching
sys.path.insert(0,"/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code")
import abundance_matching as greg_am
## Added infall_times.py to my local directory

if __name__=="__main__":
    # Take all objects within main host
    # How many merged into host long ago?
    # How many merged into host recently?
    # How many are extant in satellites?

    for fnames in glob.glob("UFDSEARCHTMP/*_7585.p"):
        if "MBs" in fnames: continue
        hid = os.path.basename(fnames)[:-7]
        hid = haloutils.hidint(hid)
        hpath = haloutils.get_hpath_lx(hid,14)
        print hid, fnames
        sys.stdout.flush()
        
        rscat = haloutils.load_rscat(hpath,haloutils.get_lastsnap(hpath))
        zoomid = haloutils.load_zoomid(hpath)
        host = rscat.ix[zoomid]
        subs = rscat.get_all_subhalos_within_halo(zoomid)
        subids = set(subs['id'])
        
        with open(fnames,"r") as f:
            data = pickle.load(f)
        relic_indices, relic_rows, merged_relic_flags, merged_relic_flags2 = data
        
        mtkeys = merged_relic_flags.keys()
        in_main_host = [mtkey in subids for mtkey in mtkeys]
        in_main_host = dict(zip(mtkeys, in_main_host))
        
        for mtkey, flags in iteritems(merged_relic_flags):
            if in_main_host[mtkey]: pass
            pass
        
        pass

