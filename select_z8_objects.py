import numpy as np
import haloutils
import time,sys
import pandas as pd
#import matplotlib.pyplot as plt
import cPickle as pickle

import MTanalysis3 as mta
import MTaddition as mtadd
from FindMiniHalos import mcrit

import sys
sys.path.append("./greg_dwarfs")
import DwarfMethods as dm
import abundance_matching

def select_z8_objects(hpath):
    hid = haloutils.get_parent_hid(hpath)
    
    h0 = 0.6711
    z_r = 8.346
    snap_r = 67

    start = time.time()
    mtc = haloutils.load_zoom_mtc(hpath, indexbyrsid=True)
    zoomid = haloutils.load_zoomid(hpath) #host halo id at z=0
    rscat = haloutils.load_rscat(hpath, haloutils.get_lastsnap(hpath), rmaxcut=False)
    subids = rscat.get_all_subhalos_within_halo(zoomid)
    print "Time to load zoom mtc and rscat: {:.1f}".format(time.time()-start)
    
    ## Extract all z=8 objects in the history of host and all its subs
    big_table = []
    for mtkey, mt in mtc.Trees.iteritems():
        ii_zr  = np.array(mt['snap'] == snap_r)
        ii_big = np.array(mt['mvir']/h0 >= 10**6.5)
        if mt[0]['mvir']/h0 > 1e7:
            print mtkey, np.log10(mt[0]['mvir']/h0), np.sum(ii_zr), np.sum(ii_zr & ii_big)
        df = pd.DataFrame(mt[ii_zr & ii_big])
        df['mtkey'] = mtkey
        big_table.append(df)
    final_table = pd.concat(big_table, ignore_index=True)
    np.save("UFDSEARCH_Z0/{}_z8halos.npy".format(haloutils.hidstr(hid)), final_table.to_records(index=False))
    
if __name__=="__main__":
    hpaths = dm.get_hpaths(field=False, lx=14)
    for hpath in hpaths:
        select_z8_objects(hpath)
    #hpath = hpaths[0]  # or whatever hpath you want
    #hpath = hpaths[4]
    #hpath = hpaths[7]
    #print hpath

