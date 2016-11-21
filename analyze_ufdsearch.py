import glob
import haloutils
import os, sys
import cPickle as pickle
from six import iteritems
import numpy as np
import matplotlib.pyplot as plt

import pandas as pd

allNtot = 0; allN = 0
Nhalos = 0
halos_with_mergers = 0

Ntotarr = []
Narr = []
current_merged_halos = {}
all_current_halos = {}

MWhalos =  {}

num_merged_map = {}
num_cands_map = {}

def make_str(hid,key):
    return "{}_{}".format(hid,key)
def make_tup(strkey):
    hid,key = strkey.split("_")
    return (int(hid),int(key))

num_in_host_tree = []
all_subs_with_relics = {}
#logmass_bins = np.arange(4,13,.2)
#H = np.zeros(len(logmass_bins)-1)

for fnames in glob.glob("UFDSEARCHTMP/*.p"):
    hid = os.path.basename(fnames)[:-2]
    hid = haloutils.hidint(hid)
    hpath = haloutils.get_hpath_lx(hid,14)
    rscat = haloutils.load_rscat(hpath,haloutils.get_lastsnap(hpath))
    zoomid = haloutils.load_zoomid(hpath)
    rscat['parentsim_hid'] = np.zeros(len(rscat),dtype=int)+hid
    new_keys = [make_str(hid,key) for key in rscat.index]
    rscat.data.index = new_keys
    MWhalos[hid] = rscat.ix[make_str(hid,zoomid)]

    subids = rscat['id'][rscat['hostID']==zoomid]
    subids_with_relics = set([])

    has_mergers = 0
    with open(fnames,"r") as f:
        data = pickle.load(f)
    relic_indices, relic_rows, merged_relic_flags, merged_relic_flags2 = data
    Ntot=0; N=0
    print fnames
    for mtkey, flags in iteritems(merged_relic_flags):
        N += np.sum(flags); Ntot += len(flags)
        all_current_halos[(hid,mtkey)] = rscat.ix[make_str(hid,mtkey)]
        num_cands_map[(hid,mtkey)] = len(flags)
        num_merged_map[(hid,mtkey)] = np.sum(flags)
        if np.sum(flags) > 0: 
            print "   ",mtkey,np.sum(flags)
            current_merged_halos[(hid,mtkey)] = rscat.ix[make_str(hid,mtkey)]
            has_mergers = 1
        if mtkey == zoomid: num_in_host_tree.append(np.sum(flags))
        if mtkey in subids:
            subids_with_relics.add(mtkey)
    print N,Ntot
    Ntotarr.append(Ntot)
    Narr.append(N)
    allNtot += Ntot; allN += N
    Nhalos += 1
    halos_with_mergers += has_mergers
    
    subs_with_relics = rscat.ix[np.array(subids_with_relics)]
    all_subs_with_relics[hid] = subs_with_relics

print allN, allNtot, halos_with_mergers,Nhalos

dfMWhalos = pd.DataFrame(MWhalos.values())
dfhalos = pd.DataFrame(current_merged_halos.values())
dfallhalos = pd.DataFrame(all_current_halos.values())
#dfhalos['num_merged'] = [
dfhalos['logmgrav'] = np.log10(dfhalos['mgrav'])
dfhalos['logmvir'] = np.log10(dfhalos['mvir'])
dfallhalos['logmgrav'] = np.log10(dfhalos['mgrav'])
dfallhalos['logmvir'] = np.log10(dfhalos['mvir'])

frac_merged = np.array(Narr).astype(float)/np.array(Ntotarr)
dx = .01
plt.hist(frac_merged,bins=np.arange(0.8,1+dx,dx))
plt.show()

