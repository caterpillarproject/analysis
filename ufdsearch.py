import numpy as np
import haloutils
import time,sys
from caterpillaranalysis import PluginBase
import cPickle as pickle

import readhalos.RSDataReader as RDR
import mergertrees.MTCatalogue as MTC
import readsnapshots.readsnapHDF5_greg as rsg

import multiprocessing


# Py 2 and 3 dictionary iteration
from six import iteritems

def search_mt_for_relics(relics,mt,snap):
    """
    Given a RSDataReader dataframe relics (from snap),
    search a given merger tree for relics with their IDs.
    Return the row in that tree.
    """

    ii = mt['snap']==snap
    indices_into_mt = np.where(ii)[0]
    _mt = mt[ii]
    
    mtrows = np.ones(len(relics))*-1
    if len(_mt['origid']) == 0: return mtrows

    for i,rsid in enumerate(np.array(relics['id']).astype(int)):
        diffs = np.abs(_mt['origid'] - rsid)
        ix = np.argmin(diffs)
        if diffs[ix] == 0:
            mtrows[i] = indices_into_mt[ix]
    return mtrows

def search_mtc_for_relics(relics, mtc, snap):
    """
    Iterate over mt in a merger tree catalog, looking for row of matches
    """
    relic_indices = {}
    for j,mt in iteritems(mtc.Trees):
        mtrows = search_mt_for_relics(relics,mt,snap)
        if np.any(mtrows != -1):
            relic_indices[j] = mtrows
    return relic_indices

def get_id2row_map(mt, idtype="tree"):
    """
    Return dict for tree or rsid to a row in the table.
    Sentinel value of -1 returns -1.
    """
    idtypes = {"tree":"id",
               "rsid":"origid"}
    assert idtype in idtypes
    idcol = idtypes[idtype]
    out = dict(zip(mt[idcol], np.arange(len(mt[idcol]))))
    out[-1] = -1
    return out
              
def trace_descendants(row, mt, id2row=None):
    """
    Given a row in a mt, get the rows of the descendants.
    """
    if id2row is None: id2row = get_id2row_map(mt)
    
    rowlist = [row]
    next_desc_row = id2row[mt[row]['desc_id']]
    while next_desc_row != -1:
        row = next_desc_row
        rowlist.append(row)
        next_desc_row = id2row[mt[row]['desc_id']]
    return rowlist

def get_relic_rows(relic_indices):
    relic_rows = {}
    for mtkey,relics_in_tree in iteritems(relic_indices):
        this_relic_rows = relics_in_tree[relics_in_tree != -1].astype(int)
        relic_rows[mtkey] = this_relic_rows
    return relic_rows

def check_if_relic_row_merged(row, mt, id2row):
    rowlist = trace_descendants(row, mt, id2row=id2row)
    rowlist = rowlist[1:]
    rows_where_experienced_merger = np.where(mt[rowlist]['num_prog'] > 1)[0]
    # Check each node with mergers to see if any mergers are larger than MX
    # In other words, if any non-mmp have M > MX, then there is a merger
    for merger_row in rows_where_experienced_merger:
        ii_prog = mt['desc_id'] == mt[merger_row]['id']
        progs = mt[ii_prog]
        assert np.sum(progs['mmp'])==1 or len(progs)==0,np.sum(progs['mmp'])
        # If has a large merger:
        if np.any(progs[progs['mmp']==0]['mvir'] > MX):
            # Only need one large merger to disqualify this guy
            return True #this_merged_flag[j] = True
    return False
    
def tag_relics_as_merged(relic_rows, mtc, MX=None, verbose=False, nthreads=1):
    if MX is None: MX = 10**9.0
    merged_relics = {}
    for mtkey,this_relic_rows in iteritems(relic_rows):
        start = time.time()
        mt = mtc[mtkey]
        this_merged_flag = np.zeros(len(this_relic_rows),dtype=bool)
        # Loop over all relics in this tree
        id2row = get_id2row_map(mt)
        for j,row in enumerate(this_relic_rows):
            rowlist = trace_descendants(row, mt, id2row=id2row)
            lastid = mt[rowlist[0]]['id']
            rowlist = rowlist[1:]
            rows_where_experienced_merger = np.where(mt[rowlist]['num_prog'] > 1)[0]
            # Check each node with mergers to see if any mergers are larger than MX
            # In other words, if any non-mmp have M > MX, then there is a merger
            start2 = time.time()
            for _merger_row in rows_where_experienced_merger:
                merger_row = rowlist[_merger_row]
                if verbose:
                    print "    row={} merger_row={}".format(row,merger_row)
                ii_prog = mt['desc_id'] == mt[merger_row]['id']
                progs = mt[ii_prog]
                #assert np.sum(progs['mmp'])==1 or len(progs)==0,np.sum(progs['mmp'])
                if len(progs) != mt[merger_row]['num_prog'] or len(progs)==0:
                    print "      MTERROR mtkey={} row={} merger_row={} progs={} num_prog={}".format(mtkey,row,merger_row,len(progs),mt[merger_row]['num_prog'])

                # If has a large merger:
                max_mass = np.max(progs[progs['id'] != lastid]['mvir'])
                if max_mass > MX:
                    # Only need one large merger to disqualify this guy
                    this_merged_flag[j] =True
                    if verbose:
                        print "      found merger: a={} M={}".format(mt[merger_row]['scale'],max_mass)
                    break
                lastid = mt[merger_row]['id']
            if verbose:
                print "    {:.2f} {}".format(time.time()-start2,len(rows_where_experienced_merger))
        merged_relics[mtkey] = this_merged_flag
        if verbose:
            print mtkey,"{:.2f} {}".format(time.time()-start,len(this_relic_rows))
    return merged_relics

def ufdsearch(hid,lx,z_r,MX):
#if __name__=="__main__":
#    hid=1387186;lx=14;z_r=8;MX=1e9
#    hid=5320;lx=14;z_r=8;MX=1e9

    print "Running on {} LX{}".format(hid,lx)
    sys.stdout.flush()

    hpath = haloutils.get_hpath_lx(hid, lx)
    reion_redshifts = [z_r]
    reion_snaps = haloutils.get_snap_z(hpath, reion_redshifts)
    snap_r = reion_snaps[0]
    #z_r = reion_redshifts[0]

    start = time.time()
    rscat_reion = haloutils.load_rscat(hpath, snap_r)
    print "Load rscat {:.2f}".format(time.time()-start)
    sys.stdout.flush()
    
    #ii_good = np.logical_and(np.log10(rscat_reion['mgrav']) >= 7.5,
    #                         np.log10(rscat_reion['mgrav'] <= 8.5))
    #ii_good = np.logical_and(np.log10(rscat_reion['mgrav']) >= 7.6,
    #                         np.log10(rscat_reion['mgrav'] <= 8.5))
    #ii_good = np.logical_and(np.log10(rscat_reion['mgrav']) >= 7.7,
    #                         np.log10(rscat_reion['mgrav'] <= 8.5))
    #ii_good = np.logical_and(np.log10(rscat_reion['mgrav']) >= 7.8,
    #                         np.log10(rscat_reion['mgrav'] <= 8.5))
    #ii_good = np.logical_and(np.log10(rscat_reion['mgrav']) >= 8,
    #                         np.log10(rscat_reion['mgrav'] <= 8.5))
    #ii_good = np.array(rscat_reion['vmax'] >= 9.48535156)

    relic_candidates = rscat_reion[ii_good]
    
    start = time.time()
    mtc = haloutils.load_mtc(hpath,indexbyrsid=True)
    print "Load mtc {:.2f}".format(time.time()-start)
    sys.stdout.flush()
    
    ## relic_indices matches rscat at z_r to mtc
    start = time.time()
    relic_indices = search_mtc_for_relics(relic_candidates, mtc, snap_r)
    print "Search mtc for relics {:.2f}".format(time.time()-start)
    sys.stdout.flush()
#    N_found = 0
#    for key in relic_indices:
#        _N = np.sum(relic_indices[key] != -1)
#        N_found += _N
#        #print key, _N
#    print "Lost {} candidates (not in MTC)".format(len(relic_candidates) -N_found)
    
    ## relic_rows is the rows of the MT that are relics (ignoring merging)
    relic_rows = get_relic_rows(relic_indices)
    N_found = 0
    for key in relic_rows:
        N_found += len(relic_rows[key])
    print "Lost {} candidates (not in MTC)".format(len(relic_candidates) -N_found)
    sys.stdout.flush()
    
    ## merged_relic_flags is the objects that have merged with > MX
    start = time.time()
    merged_relic_flags = tag_relics_as_merged(relic_rows, mtc, MX=MX)#,verbose=True)
    print "Time to find merged relics {:.2f}".format(time.time()-start)
    N_merged = 0
    for key in merged_relic_flags:
        N_merged += np.sum(merged_relic_flags[key])
    print "{}/{} are merged with MX ({:.2f})".format(N_merged,N_found,float(N_merged)/N_found)
    sys.stdout.flush()

    start = time.time()
    merged_relic_flags_any = tag_relics_as_merged(relic_rows, mtc, MX=1.0)#,verbose=True)
    print "Time to find merged relics {:.2f}".format(time.time()-start)
    N_merged = 0
    for key in merged_relic_flags_any:
        N_merged += np.sum(merged_relic_flags_any[key])
    print "{}/{} are merged with anything ({:.2f})".format(N_merged,N_found,float(N_merged)/N_found)
    sys.stdout.flush()

    #with open("UFDSEARCHTMP/H{}_7585.p".format(hid),'w') as f:
    #    pickle.dump([relic_indices,relic_rows,merged_relic_flags,merged_relic_flags_any],f)
    with open("UFDSEARCHTMP/H{}_7685.p".format(hid),'w') as f:
        pickle.dump([relic_indices,relic_rows,merged_relic_flags,merged_relic_flags_any],f)
    #with open("UFDSEARCHTMP/H{}_7785.p".format(hid),'w') as f:
    #    pickle.dump([relic_indices,relic_rows,merged_relic_flags,merged_relic_flags_any],f)
    #with open("UFDSEARCHTMP/H{}_7885.p".format(hid),'w') as f:
    #    pickle.dump([relic_indices,relic_rows,merged_relic_flags,merged_relic_flags_any],f)
    #with open("UFDSEARCHTMP/H{}_vmax12gd.p".format(hid),'w') as f:
    #    pickle.dump([relic_indices,relic_rows,merged_relic_flags,merged_relic_flags_any],f)

if __name__=="__main__":
#def tmp():
    #hid = 1387186
    lx = 14
    z_r = 8.0
    #z_r = 12.0
    MX = 1.e9

    good_cids = [1,2,3,4,5,6,8,9,10,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,29,31,33,36,37]
    if len(sys.argv)==2:
        num = int(sys.argv[1])
        assert num in [1,2,3,4,5,6]
        start = (num-1)*6 #ooooops should have been 5
        if num==6: # oops messed up hahaha
            good_cids = [good_cids[x] for x in [5+6*y for y in range(5)]]
        else:
            good_cids = good_cids[start:start+5]
    print "CIDs",good_cids

    for cid in good_cids:
        hid = haloutils.cid2hid[cid]
        try:
            ufdsearch(hid,lx,z_r,MX)
        except Exception as e:
            print e
            print "---H{} failed! Skipping---".format(hid)
        print
