import haloutils as htils
import sys,os
import numpy as np
import pylab as plt
import pandas as pd
import subprocess as sub

lx=14 
suffix = 'all_minihalos'

def get_suffix(suffix):
    if suffix == 'all_minihalos':
        return 'all_minihalo_array.npy'
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

def run_get_next(suffix, lx,rerun=False):
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

        try:
            if rerun: # guarantee an exception, so it triggers re-running
                next_ids = np.load("/failure/minihalo.npy")
            else:
                next_ids = np.load(hpath+"/analysis/minihalo_descendants_"+suffix+'.npy')
        except:
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


def run_get_accretion(suffix, lx,rerun=False):
    import MTanalysis3 as mta
    #hpaths = htils.get_paper_paths_lx(lx)
    hpaths = htils.get_all_halo_paths_lx(lx)
    for hpath in hpaths:
        hid = htils.get_parent_hid(hpath)
        print hid
        sys.stdout.flush()
        mhalos_path = get_suffix(suffix)
        try:
            print hpath+"/analysis/"+mhalos_path
            minihalos = np.load(hpath+"/analysis/"+mhalos_path)
        except:
            print 'no minihalo data for', hid, 'skipping it'
            continue

        try:
            if rerun: # guarantee an exception, so it triggers re-running
                next_ids = np.load("/failure/minihalo.npy")
            else:
                next_ids = np.load(hpath+"/analysis/accrete_snap_"+suffix+'.npy')
        except:
            mtc = htils.load_mtc(hpath, indexbyrsid=True) #haloids=[]
            print 'loaded mtc', htils.hpath_name(hpath)
            sys.stdout.flush()
            # create mapping to trees
            tree_id_map = {}
            for value,key in enumerate(minihalos['base_rsid']):
                tree_id_map.setdefault(key, []).append(value)

            infall_snap = [-1]*len(minihalos)
            host_merge_snap = [-1]*len(minihalos)
            
            hostkey= htils.load_zoomid(hpath)
            hosttree=mtc.Trees[hostkey]
            mmp_map = hosttree.get_mmp_map()
            host_mb = hosttree.getMainBranch(0, mmp_map) # main branch of host starting from z=0 snapshot
            for i,key in enumerate(tree_id_map.iterkeys()):
                tree=mtc.Trees[key]
                desc_map = tree.get_desc_map()
                for minirow in tree_id_map[key]:
                    minihalo = minihalos[minirow]
                    desc=tree.getDescBranch(minihalo['row'],desc_map)
                    iloc,isnap = mta.getInfall(desc[::-1],host_mb) 
                    if isnap is None:
                        infall_snap[minirow] = -2
                    else:
                        infall_snap[minirow] = isnap

                    # if merged with host, desc['id'] == host_mb['id']
                    branch_len = min(len(desc),len(host_mb))
                    inhost = np.where(desc[::-1]['id'][0:branch_len]==host_mb['id'][0:branch_len])[0]
                    if len(inhost)==0:
                        host_merge_snap[minirow] = -2
                    else:
                        host_merge_snap[minirow] = host_mb['snap'][inhost[-1]]

                if i%500==0:
                    print i
                    sys.stdout.flush()
            np.save(hpath+'/analysis/accrete_snap_'+suffix, infall_snap)
            np.save(hpath+'/analysis/host_merge_snap_'+suffix, host_merge_snap)
            print 'done'


def load_next_halos(hpath,suffix):
    mhalos_path = get_suffix(suffix)
    minihalos = np.load(hpath+"/analysis/"+mhalos_path)
    next_ids = np.load(hpath+"/analysis/minihalo_descendants_"+suffix+'.npy')
    form_snap = np.load(hpath+'/analysis/form_snap_'+suffix+'.npy')
    merge_snap = np.load(hpath+'/analysis/merge_snap_'+suffix+'.npy')
    host_merge_snap = np.load(hpath+'/analysis/host_merge_snap_'+suffix+'.npy')
    accrete_snap = np.load(hpath+'/analysis/accrete_snap_'+suffix+'.npy')
    return minihalos, next_ids, np.array(form_snap), np.array(merge_snap),np.array(host_merge_snap),np.array(accrete_snap)

"""
run_get_accretion('all_chabrier', lx)
print 'done with chabrier'
run_get_accretion('all_firstgal', lx)
print 'done with firstgal'
run_get_accretion('all_kroupa', lx)
print 'done with kroupa'
run_get_accretion('all_salpeter', lx)
print 'done with salpeter'
"""

#run_get_next(suffix, lx)
#run_get_next('all_chabrier', lx)
#print 'done with chabrier'
#run_get_next('all_firstgal', lx)
#print 'done with firstgal'
#run_get_next('all_kroupa', lx)
#print 'done with kroupa'
#run_get_next('all_salpeter', lx)
#print 'done with salpeter'


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

