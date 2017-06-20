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

if __name__=="__main__":

    am_model = greg_am.GarrisonKimmel()
    mpeak_upper = am_model.stellar_to_halo_mass(2e5)
    mpeak_lower = am_model.stellar_to_halo_mass(1e3)
    vmax_ach = 9.48535156
    vmax_filt= 23.54248047
    
    with open("UFDSEARCH_ANALYZE_OBJS.pkl",'r') as f:
        all_recorded_data = pickle.load(f)
        
    all_hids = all_recorded_data.keys()
    abins = np.linspace(0,1,21)
    Mbins = np.arange(5.0, 13.0, .2)
    
    Hs = []
    HMs= []
    num_host = []
    num_subs = []
    num_field= []
    num_total= []
    num_all  = []
    num_ufds = []

    HaMs = []
    for hid in all_hids:
        #print "H"+str(hid)
        data = all_recorded_data[hid]
        
        num_host_relics = data[10]
        num_sub_relics = np.sum(data[0])
        num_field_relics = np.sum(data[4])
        num_ufds.append(data[9])
        num_host.append(num_host_relics)
        num_subs.append(num_sub_relics)
        num_field.append(num_field_relics)
        num_total.append(num_host_relics+num_sub_relics)
        num_all.append(num_host_relics+num_sub_relics+num_field_relics)

        record_host_mergetimes = data[13]
        record_host_mergemasses= data[14]
        record_host_mergempeaks= data[15]

        h,x = np.histogram(record_host_mergetimes, bins=abins)
        Hs.append(h)
        h,x = np.histogram(np.log10(record_host_mergempeaks), bins=Mbins)
        HMs.append(h)
        
        h,x1,x2 = np.histogram2d(record_host_mergetimes, np.log10(record_host_mergempeaks), bins=[abins,Mbins])
        HaMs.append(h)
    Hs = np.array(Hs)
    Hmean = np.mean(Hs,axis=0)
    Hstd  = np.std(Hs,axis=0)
    Htot = np.sum(Hs,axis=0)
    
    HMs = np.array(HMs)
    HMtot = np.sum(HMs,axis=0)

    HaMs = np.array(HaMs)
    HaMtot = np.sum(HaMs, axis=0)

    num_host = np.array(num_host)*1.0
    num_subs = np.array(num_subs)*1.0
    num_ufds = np.array(num_ufds)*1.0
    num_field= np.array(num_field)*1.0
    num_total= np.array(num_total)*1.0
    num_all= np.array(num_all)*1.0
    
    frac_host = num_host/num_total
    frac_subs = num_subs/num_total
    frac_ufds = num_ufds/num_total
    
    print "{:.3f} +/- {:.3f} in host".format(np.mean(frac_host), np.std(frac_host))
    print "{:.3f} +/- {:.3f} in subs".format(np.mean(frac_subs), np.std(frac_subs))
    print "{:.3f} +/- {:.3f} in ufds".format(np.mean(frac_ufds), np.std(frac_ufds))

    
    fig,ax = plt.subplots()
    amid = (abins[1:]+abins[:-1])/2.
    ax.plot(amid, Htot, drawstyle='steps-mid')
    fig,ax = plt.subplots()
    Mmid = (Mbins[1:]+Mbins[:-1])/2.
    ax.plot(Mmid,HMtot, drawstyle='steps-mid')

    fig, ax = plt.subplots()
    im = ax.imshow(HaMtot.T,aspect='auto',interpolation='none',
              origin='lower', extent=[0,1,5,13], cmap='bone_r')
    fig.colorbar(im)
    plt.show()
    
