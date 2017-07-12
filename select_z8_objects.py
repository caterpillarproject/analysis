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

zin2snapr = {
    0:319,
    1:171,
    2:137,
    3:119,
    4:108,
    5:100,
    6:90,
    7:78,
    8:67,
    9:59,
    10:53,
    11:47,
    12:42,
    13:38,
    14:34
}
all_zin = np.sort(zin2snapr.keys())

def zin_to_zr_snapr(zin,verbose=True):
    #### From MTaddition.py:
    #### zsnaps = [90,78,67,59,53,47,42,38,34]  # corresponds to z = 6.33, 7.26, 8.346, 9.33, 10.22, 11.28, 12.33, 13.31, 14.44
    ## Get snap_r
    assert zin in all_zin, zin
    snap_r = zin2snapr[zin]
    # Hardcode a random LX14 halo
    hpath = haloutils.get_hpath_lx(5320,14)
    z_r = haloutils.get_z_snap(hpath,snap_r)[0]
    #if zin == 4:
    #    z_r = 4.26
    #    snap_r = 108
    #elif zin == 6:
    #    z_r = 6.33
    #    snap_r = 90
    #elif zin == 8:
    #    z_r = 8.346
    #    snap_r = 67
    #elif zin == 10:
    #    z_r = 10.22
    #    snap_r = 53
    #elif zin == 12:
    #    z_r = 12.33
    #    snap_r = 42
    if verbose:
        print "z = {}".format(zin)
    return z_r, snap_r

def select_z8_objects(hpath,zin=8):

    hid = haloutils.get_parent_hid(hpath)
    
    h0 = 0.6711
    z_r, snap_r = zin_to_zr_snapr(zin)
    #z_r = 8.346
    #snap_r = 67

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
    np.save("UFDSEARCH_Z0/{}_z{}halos.npy".format(haloutils.hidstr(hid),zin), final_table.to_records(index=False))
    
if __name__=="__main__":
    zin = 4
    hpaths = dm.get_hpaths(field=False, lx=14)
    for hpath in hpaths:
        select_z8_objects(hpath,zin=zin)
    #hpath = hpaths[0]  # or whatever hpath you want
    #hpath = hpaths[4]
    #hpath = hpaths[7]
    #print hpath

