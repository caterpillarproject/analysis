import haloutils as htils
import sys,os
import numpy as np
import pylab as plt
import pandas as pd
import subprocess as sub

lx=14 
suffix = 'all_firstgal'

def get_suffix(suffix):
    if suffix == 'all_chabrier':
        return 'all_lw_minihalo_array_chabrier.npy'
    if suffix == 'all_firstgal':
        return 'all_firstgal_array.npy'
    if suffix == 'all_kroupa':
        return 'all_lw_minihalo_array_kroupa.npy'
    if suffix == 'all_salpeter':
        return 'all_lw_minihalo_array_salpeter.npy'
    else:
        return ''

def run_get_next(suffix, lx):
    #hpaths = htils.get_paper_paths_lx(lx)
    hpaths = htils.get_all_halo_paths_lx(lx)
    N=20
    for hpath in hpaths:
        hid = htils.get_parent_hid(hpath)
        print hid
        mhalos_path = get_suffix(suffix)
        try:
            print hpath+"/analysis/"+mhalos_path
            minihalos = np.load(hpath+"/analysis/"+mhalos_path)
        except:
            print 'no minihalo data for', hid, 'skipping it'
            continue
        mtc = htils.load_mtc(hpath, indexbyrsid=True) #haloids=[]
        print 'loaded mtc', htils.hpath_name(hpath)
        # create mapping to trees
        tree_id_map = {}
        for value,key in enumerate(minihalos['base_rsid']):
            tree_id_map.setdefault(key, []).append(value)

        next_ids = [[]]*len(minihalos)
        merger = [-1]*len(minihalos)
        form_snap=[-1]*len(minihalos)
        merge_snap=[-1]*len(minihalos)
        for i,key in enumerate(tree_id_map.iterkeys()):
            tree=mtc.Trees[key]
            desc_map = tree.get_desc_map()
            for minirow in tree_id_map[key]:
                minihalo = minihalos[minirow]
                form_snap[minirow] = minihalo['snap']
                desc=tree.getDescBranch(minihalo['row'],desc_map)
                next_ids[minirow] = desc['origid'][1:1+N]
                mask = desc['mmp']!=1
                if np.sum(mask)==0:
                    merge_snap[minirow]=desc[-1]['snap']+1
                else:
                    merge_snap[minirow] = desc[mask][0]['snap']+1
            if i%500==0:
                print i
        np.save(hpath+'/analysis/minihalo_descendants_'+suffix, next_ids)
        np.save(hpath+'/analysis/form_snap_'+suffix, form_snap)
        np.save(hpath+'/analysis/merge_snap_'+suffix, merge_snap)
        print 'done'


def load_next_halos(hpath,suffix):
    mhalos_path = get_suffix(suffix)
    minihalos = np.load(hpath+"/analysis/"+mhalos_path)
    next_ids = np.load(hpath+"/analysis/minihalo_descendants_"+suffix+'.npy')
    form_snap = np.load(hpath+'/analysis/form_snap_'+suffix+'.npy')
    merge_snap = np.load(hpath+'/analysis/merge_snap_'+suffix+'.npy')
    return minihalos, next_ids, form_snap, merge_snap
#newdata=np.load(hpath+'/analysis/minihalo_descendants.npy')





"""
for minihalo,i in zip(minihalos,range(len(minihalos))):
    print get_next_N_ids(mtc,minihalo,N=20)
"""

#next_ids = []
#next_ids = np.ndarray((len(minihalos),20),dtype=np.int32)
#next_ids = [[-1]*20]*len(minihalos)

#def get_next_N_ids(mtc, minihalo, N,desc_map=None):
#    tree = mtc.Trees[minihalo['base_rsid']]
#    row = minihalo['row']
#    return tree.getDescBranch(row,desc_map)['origid'][1:1+N]  
# make this function faster by only going N steps in getDescBranch

