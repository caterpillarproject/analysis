from brendanlib.grifflib import getcandidatelist
import numpy as np
import sys,glob
import subprocess as sub

"""
Contact: Brendan Griffen <brendan.f.griffen@gmail.com>

This script will get the list of halos ids. It will also make
the folders for the halos if they don't already exist.

"""

MAKE_DIRECTORY_OPT = False   # Set if you want to make the directories - BE SURE!

base_path = "/bigbang/data/AnnaGroup/caterpillar/halos/"
candidates = getcandidatelist('./candidates.dat')

def get_candidates(lower_mass,upper_mass,number_of_halos):
    cmasses = candidates[:,1]

    # Don't get halos near the edge as constructing their lagr. volume is a tad more annoying.
    pos_mask = (candidates[:,4] <= 90.) & (candidates[:,4] >= 10.) \
               & (candidates[:,5] <= 90.) & (candidates[:,5] >= 10.) \
               & (candidates[:,6] <= 90.) & (candidates[:,6] >= 10.)
    
    middle_mass_mask = np.logical_and(cmasses >= lower_mass,cmasses <= upper_mass)
    mask_total = middle_mass_mask & pos_mask
    candidates_masked = candidates[mask_total]

    # Generate random index of size: number_of_halos
    idx = np.random.randint(len(candidates_masked),size=number_of_halos)

    for i in idx:
        print int(candidates_masked[i,0])

    return list(candidates_masked[:,0][idx])

print 
print "Middle bin halos"
middle_bin_halos = get_candidates(1e12,2e12,number_of_halos=48)
print
print "Higher bin halos"
right_bin_halos = get_candidates(2e12,3e12,number_of_halos=6)
print 
print "Lower bin halos"
left_bin_halos = get_candidates(0.7e12,1e12,number_of_halos=6)

print
if MAKE_DIRECTORY_OPT:
    for halo in middle_bin_halos:
        cmd_make_folder = "mkdir -p " + base_path + "H" + str(int(halo))
        #print cmd_make_folder
        mask_candidate = np.where(halo == candidates[:,0])
        mass_haloi = candidates[:,1][mask_candidate][0]
        print "Mass: %3.2e" % (float(mass_haloi))
        print cmd_make_folder
        sub.call([cmd_make_folder],shell=True)
