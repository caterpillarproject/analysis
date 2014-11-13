import readsnapshots.readsnapHDF5 as rsHD
import readsnapshots.readsnap as rs
import numpy as np
import random
import readhalos.RSDataReader as RSDataReader
import pylab as plt
import sys
from random import randint

"""
Contact: Brendan Griffen <brendan.f.griffen@gmail.com>

This script will select candidate caterpillar halos from the parent simulation.
The algorithm is as follows:

1. Select all halos within 0.7-3 x 10^12 Msol.
2. Iterate through each candidate in this mass range.
2a. Find minimum distance to 7 x 10^13 Msol halo.
2b. Find minimum distance to halo 0.5*Mhost[i]
3. Ensure minimum distance to (2a) is 7 Mpc.
4. Ensure minimum distance to (2b) is 2.8 Mpc.

This effectively follows the ELVIS suite for isolated candidates.

"""

# SET DATA DIRECTORIES
basepath = '/bigbang/data/AnnaGroup/caterpillar/parent/gL100X10/'
halopath = basepath + 'rockstar'
outputfile = './candidates.dat'

header = rsHD.snapshot_header(basepath + 'outputs/snapdir_127/snap_127.0.hdf5')
hubble = header.hubble

# SELECT CANDIDATE MASS RANGE
lowermasscut = 0.7e12
uppermasscut = 3e12

# LOAD ROCKSTAR HALOS
print "READING...",halopath
halodata = RSDataReader.RSDataReader(halopath,127,version=6)
allhalos = halodata.get_hosts()

mask_mass_range = np.logical_and(allhalos['mvir']/hubble>lowermasscut, allhalos['mvir']/hubble<uppermasscut)
mw_candidates = allhalos[mask_mass_range]

# LOAD HALOS LARGER THAN X MSOL FOR ISOLATION REQUIREMENT
halos_greater_than13 = allhalos[allhalos['mvir']/hubble>7e13] # HALOS LARGER THAN 7E13

number_mw_candidates = len(mw_candidates)
print "Number of MW candidates between 7e11 and 3e12:",number_mw_candidates
print "Number of halos larger than 7e13 Msol",len(halos_greater_than13)

xpos13 = np.array(np.float64(halos_greater_than13['posX']))
ypos13 = np.array(np.float64(halos_greater_than13['posY']))
zpos13 = np.array(np.float64(halos_greater_than13['posZ']))

# OUTPUT FILE
out = open(outputfile, 'w')
out.write('# id, msol, rvir, vmax, x, y, z, vx, vy, vz, distance to first 7e13 halo, distance to halo at least half as massive\n')

number_of_candidates = 0
for i in xrange(0,number_mw_candidates):

    # CYCLE THROUGH EACH CANDIDATE
    xposi = np.array(np.float64(mw_candidates['posX']))[i]
    yposi = np.array(np.float64(mw_candidates['posY']))[i]
    zposi = np.array(np.float64(mw_candidates['posZ']))[i]
    xveli = np.array(np.float64(mw_candidates['pecVX']))[i]
    yveli = np.array(np.float64(mw_candidates['pecVY']))[i]
    zveli = np.array(np.float64(mw_candidates['pecVZ']))[i]
    massi = np.array(np.float64(mw_candidates['mvir']))[i]/hubble
    rviri = np.array(np.float64(mw_candidates['rvir']))[i]
    vmaxi = np.array(np.float64(mw_candidates['vmax']))[i]
    idi = np.array(mw_candidates['id'])[i]

    # CALCULATE DISTANCE TO HALOS LARGER THAN 7e13
    R13 = np.sqrt((xposi-xpos13)**2.+(yposi-ypos13)**2.+(zposi-zpos13)**2.)/hubble

    # SELECT ALL HALOS LARGER THAN HALF THE SIZE OF THE CANDIDATE idi
    largerthanMW = allhalos[allhalos['mvir']/hubble >= 0.5*massi]
    
    #largerthanMW = allhalos[allhalos['mvir']/hubble >= 4e12]

    # SINCE 0.5*massi INCLUDES idi, NEED TO REMOVE IT FROM X,Y,Z POS CALC SO MIN(R) != 0.0
    idindex = np.where(largerthanMW['id'] != idi)
    xtmp = np.array(np.float64(largerthanMW['posX']))
    ytmp = np.array(np.float64(largerthanMW['posY']))
    ztmp = np.array(np.float64(largerthanMW['posZ']))

    xposMgtMW = np.float64(xtmp[idindex[0]])
    yposMgtMW = np.float64(ytmp[idindex[0]])
    zposMgtMW = np.float64(ztmp[idindex[0]])

    RMgtMW = np.sqrt((xposi-xposMgtMW)**2.+(yposi-yposMgtMW)**2.+(zposi-zposMgtMW)**2.)/hubble

    # NO HALO LARGER THAN 7e13 CLOSER THAN 4 MPC AND NO HALO HALF THE MASS OF THE HOST WITHIN 2.8 MPC
    if np.min(R13) >= 7. and np.min(RMgtMW) >= 2.8:
        number_of_candidates += 1
        out.write('%i %e %f %f %f %f %f %f %f %f %f %f \n' % (int(idi),massi,rviri,vmaxi,xposi,yposi,zposi,xveli,yveli,zveli,np.min(R13),np.min(RMgtMW))) 

out.close()

print "Number of candidates found:",number_of_candidates
print "Output written to:",outputfile

print "Now use 02_get_candidate_list.py to generate folders then 03_create_lagr.py to construct lagrangian region."
