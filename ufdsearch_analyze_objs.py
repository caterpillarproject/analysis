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

import ufdsearch

def find_time_and_mass_of_mergers(relic_rows, mb, mt, verbose=False):
    merge_times = np.zeros(len(relic_rows)) - 1.0
    merge_masses = np.zeros(len(relic_rows)) - 1.0
    merge_mpeaks = np.zeros(len(relic_rows)) - 1.0
    
    relic_cands = mt[relic_rows]
    id2row = ufdsearch.get_id2row_map(mt)
    mb_treeids = set(mb['id'])
    
    for i,row in enumerate(relic_rows):
        next_desc_row = id2row[mt[row]['desc_id']]
        mpeak = mt[row]['mvir']
        while next_desc_row != -1:
            row = next_desc_row
            if mt[row]['id'] in mb_treeids: break
            mpeak = max(mpeak, mt[row]['mvir']) #hasn't merged yet
            next_desc_row = id2row[mt[row]['desc_id']]

        if next_desc_row == -1:
            #got through tree without merging...
            if verbose: print "    --no mergers--: entry {} in relic_rows".format(i)
            merge_times[i] = np.nan
            merge_masses[i] = np.nan
            merge_mpeaks[i] = np.nan
        else:
            merge_times[i] = mt[row]['scale']
            merge_masses[i] = mt[row]['mvir']
            merge_mpeaks[i] = mpeak
    return merge_times, merge_masses, merge_mpeaks

if __name__=="__main__":

    am_model = greg_am.GarrisonKimmel()
    mpeak_upper = am_model.stellar_to_halo_mass(2e5)
    mpeak_lower = am_model.stellar_to_halo_mass(1e3)
    vmax_ach = 9.48535156
    vmax_filt= 23.54248047
    
    all_recorded_data = {}
    
    suffix = "_7585"
    for fnames in glob.glob("UFDSEARCHTMP/*"+suffix+".p"):
        if "MBs" in fnames: continue
        hid = os.path.basename(fnames)[:-7]
        
        try:
            hid = haloutils.hidint(hid)
        except Exception as e:
            print "ERROR?",fnames, hid
            print e
            continue
        hpath = haloutils.get_hpath_lx(hid,14)
        print hid, fnames
        sys.stdout.flush()
        
        rscat = haloutils.load_rscat(hpath,haloutils.get_lastsnap(hpath))
        zoomid = haloutils.load_zoomid(hpath)
        host = rscat.ix[zoomid]
        subs = rscat.get_all_subhalos_within_halo(zoomid)
        subids = set(subs['id'])
        
        ##print "Loading MTC"
        ##start = time.time()
        ##mtc = haloutils.load_mtc(hpath,indexbyrsid=True)
        ##print "{:.2f}".format(time.time()-start)
        ##sys.stdout.flush()
        print "Loading zoom MTC"
        start = time.time()
        mtc = haloutils.load_zoom_mtc(hpath,indexbyrsid=True)
        print "{:.2f}".format(time.time()-start)
        sys.stdout.flush()
        
        with open(fnames,"r") as f:
            data = pickle.load(f)
        relic_indices, relic_rows, merged_relic_flags, merged_relic_flags2 = data
        
        mtkeys = merged_relic_flags.keys()
        in_main_host = [mtkey in subids for mtkey in mtkeys]
        in_main_host = dict(zip(mtkeys, in_main_host))
        
        fieldids = set(mtkeys).difference([zoomid]).difference(subids)
        
        # Sort relics based on where they are at z=0
        surviving_relics_in_host  = {} # subhalos of host, no merger with > MX
        surviving_relics_in_field = {} # other objects without mergers
        relics_merged_into_host   = {} # found in mtc[zoomid]
        relics_merged_into_subs   = {} # found in mtc[subid]
        relics_merged_into_field  = {} # found in mtc[fieldid]
        
        num_host_relics = -1
        num_host_merged = -1
        record_host_mergetimes = []
        record_host_mergemasses = []
        record_host_mergempeaks = []
        
        record_sub_relics = []
        record_sub_merged = []
        record_sub_dr     = []
        record_sub_mpeak  = []

        # TODO: how to do this?
        sub_mergetimes = []
        sub_mergemasses = []
        sub_mergempeaks = []
        
        record_field_relics = []
        record_field_merged = []
        record_field_dr     = []
        
        for mtkey, flags in iteritems(merged_relic_flags):
            num_relic_cands = len(relic_rows[mtkey])
            num_merged = np.sum(flags)

            ### Get Greg's extant data (based on DwarfMethods.py get_extant_data())
            dataE = AE.read(hpath)
            data_reion = E.read(hpath)
            extant = pd.concat((dataE, data_reion), axis=1)
            extant.index = extant['rsid']
            ## Cut depth==0 which is only the actually extant objects
            extant = extant[extant['depth']==0]

            mask_pre = np.array(extant['vmax_12'] > vmax_ach)
            mask_post= np.array(extant['peak_vmax'] >= vmax_filt)
            in_range1 = mask_pre | mask_post
            in_range2 = mask_pre
            in_range3 = np.logical_and(np.array(extant['max_mass'] < mpeak_upper),
                                       np.array(extant['max_mass'] > mpeak_lower))
            num_greg_1 = np.sum(in_range1 & in_range3)
            num_greg_2 = np.sum(in_range2 & in_range3)

            if mtkey == zoomid: # the main host
                num_host_relics = num_relic_cands
                num_host_merged = num_merged
                try:
                    mt = mtc[mtkey]
                except:
                    print "  --ERROR: cannot find host in mtc???--"
                else:
                    mb = mt.getMainBranch(0)
                    merging_times, merging_masses, merging_mpeaks = find_time_and_mass_of_mergers(relic_rows[mtkey], mb, mt)
                    record_host_mergetimes = merging_times
                    record_host_mergemasses= merging_masses
                    record_host_mergempeaks= merging_mpeaks
                    
                    ## TODO find merging time and mass of FIRST merger

            elif mtkey in subids: # sub of the host
                record_sub_relics.append(num_relic_cands)
                record_sub_merged.append(num_merged)
                record_sub_dr.append(float(subs.ix[mtkey]['dr']))
                try:
                    mpeak = float(extant.ix[mtkey]['max_mass'])
                except:
                    record_sub_mpeak.append(np.nan)
                else:
                    record_sub_mpeak.append(mpeak)

                try:
                    mt = mtc[mtkey]
                except:
                    print "  --ERROR: cannot find {} in mtc--".format(mtkey)
                    sub_mergetimes.append([])
                    sub_mergemasses.append([])
                    sub_mergempeaks.append([])
                else:
                    mb = mt.getMainBranch(0)
                    merging_times, merging_masses, merging_mpeaks = find_time_and_mass_of_mergers(relic_rows[mtkey], mb, mt)
                    sub_mergetimes.append(merging_times)
                    sub_mergemasses.append(merging_masses)
                    sub_mergempeaks.append(merging_mpeaks)
                    
            elif mtkey in fieldids: # field halo
                record_field_relics.append(num_relic_cands)
                record_field_merged.append(num_merged)
                try:
                    dr = float(rscat.ix[mtkey]['dr'])
                except:
                    print "  --ERROR: found drcut halo in field {}--".format(mtkey)
                    record_field_dr.append(np.nan)
                else:
                    record_field_dr.append(dr)
            else:
                print "ERROR: MTKEY NOT RECOGNIZED!!!"
        
        record_sub_relics = np.array(record_sub_relics)
        record_sub_merged = np.array(record_sub_merged)
        record_sub_dr     = np.array(record_sub_dr)
        record_sub_mpeak  = np.array(record_sub_mpeak)
        record_field_relics = np.array(record_field_relics)
        record_field_merged = np.array(record_field_merged)
        record_field_dr     = np.array(record_field_dr)

        print "  The host (logM={:.2f}) has {} relics of which {} merge".format(np.log10(host['mgrav']),num_host_relics, num_host_merged)
        print "  There are {} subs with relics, total {} of which {} merge".format(len(record_sub_relics), np.sum(record_sub_relics), np.sum(record_sub_merged))
        print "    {} subs have relics but no merged relics".format(np.sum(record_sub_merged==0))
        num_mine = np.sum(np.logical_and(np.logical_and(record_sub_merged==0,record_sub_mpeak > mpeak_lower), record_sub_mpeak < mpeak_upper))
        print "    {} of those are in the UFD mass range".format(num_mine)
        print "    Greg gets {} (pre) and {} (pre+post)".format(num_greg_2, num_greg_1)
        print "    {} subs have bad mpeak".format(np.sum(np.isnan(record_sub_mpeak)))
        print "  There are {} field halos with relics, total {} of which {} merge".format(len(record_field_relics), np.sum(record_field_relics), np.sum(record_field_merged))
        
        all_recorded_data[hid] = [record_sub_relics,
                                  record_sub_merged,
                                  record_sub_dr    ,
                                  record_sub_mpeak ,
                                  record_field_relics,
                                  record_field_merged,
                                  record_field_dr    ,
                                  num_greg_1, num_greg_2, num_mine,
                                  num_host_relics, np.sum(record_sub_relics), np.sum(record_sub_merged),
                                  record_host_mergetimes, record_host_mergemasses, record_host_mergempeaks,
                                  sub_mergetimes, sub_mergemasses, sub_mergempeaks
                                  ]
    

    with open("UFDSEARCH_ANALYZE_OBJS.pkl",'w') as f:
        pickle.dump(all_recorded_data, f)

        # Pull out z=0 objects containing relics
        # Filter by those in the host
        #subz0_relics = rscat.ix[subids_of_relics_in_host]
        #num_relics_in_host = len(subz0_relics)
        
        #fieldz0_relics = rscat.ix[rsids_of_relics_in_field]
        #num_relics_in_field = len(fieldz0_relics)
        
        ## TODO: split of surviving, merged into main host, merged into subhalos

        ## TODO: mass/radial distribution of surviving objects
        
        ## TODO: merging time of host relics into main host
        ## TODO: merging time of host relics into any host
        
        ## TODO: merging time of sub relics into subhalos
        ## TODO: mass distribution of subhalos with sub relics
        
