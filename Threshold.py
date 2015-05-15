import numpy as np
from caterpillaranalysis import *
import readsnapshots.readsnapHDF5_greg as rsg
import haloutils
import os
from caterpillarplot import *

# input
param = 'vmax'
thresh = 5.0
redshift = 6.0


lx=14
halo_paths = haloutils.find_halo_paths(levellist=[lx],require_mergertree=True,require_subfind=False,verbose=False) 
hpath = halo_paths[0]
# end input


mtc = haloutils.load_mtc(hpath) #,haloids=[hostID])
trees = mtc.Trees
firsthalos = []

for tree in trees:
    # go back to redshift 6 halos
    a = 1/(1+redshift)
    a = mtc.scale_list[np.where(mtc.scale_list>a)[0][0]]
    epsilon = (mtc.scale_list[-1]-mtc.scale_list[-2])/10
    mask = np.logical_and(tree['scale']<a+epsilon, tree['scale']>a-epsilon)
    halos = np.where(mask)[0] # all z=6 halos

    for halo in halos:
        get_halos_recurse(halo, pre_value=None, pre_row=None, mmp=False)

def get_halos_recurse(halo, pre_value,pre_row,mmp):    
    cur_value = tree[halo][param]
    print halo, cur_value, tree[halo]['scale']
    if cur_value < thresh:
        if np.sum(tree.getMainBranch(halo)[param]>thresh) == 0:
            if pre_value == None or pre_value < thresh: # should never get to 2nd clause
                return # always too small
            elif mmp: 
                print 'found a halo'
                # pre_row is the line of the halo that crossed the threshold
                # to get this halo back out of the MT, load the Tree with "tree.rockstar_id", and go to line "pre_row"
                firsthalos.append([tree.rockstar_id,pre_row, tree[pre_row]['origid'], tree[pre_row]['snap']]) # pre_value, cur_value
                return
    mmp = tree.getMMP(halo)
    if mmp!=None:
        get_halos_recurse(mmp, cur_value, halo, mmp=True)
        progenitors = tree.getNonMMPprogenitors(halo)
        for progenitor in progenitors:
            get_halos_recurse(progenitor, cur_value, halo, mmp=False)
    
