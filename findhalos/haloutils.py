import numpy as np
import os
import sys
import subprocess
from optparse import OptionParser

import mergertrees.MTCatalogue as MTC

def get_numsnaps(outpath):
    return sum(1 for line in open(outpath+'/ExpansionList'))

def get_foldername(outpath):
    return os.path.basename(os.path.normpath(outpath))

def get_parent_hid(outpath):
    hidstr = get_foldername(outpath).split('_')[0]
    return int(hidstr[1:])

def get_zoom_params(outpath):
    """ return ictype, LX, NV """
    split = get_foldername(outpath).split('_')
    return split[1],int(split[5][2:]),int(split[7][2:])

def check_last_subfind_exists(outpath):
    numsnaps = get_numsnaps(outpath)
    lastsnap = numsnaps - 1; snapstr = str(lastsnap).zfill(3)
    group_tab = os.path.exists(outpath+'/outputs/groups_'+snapstr+'/group_tab_'+snapstr+'.0')
    subhalo_tab = os.path.exists(outpath+'/outputs/groups_'+snapstr+'/subhalo_tab_'+snapstr+'.0')
    return group_tab and subhalo_tab

def check_last_rockstar_exists(outpath,particles=False):
    numsnaps = get_numsnaps(outpath)
    lastsnap = numsnaps - 1; snapstr = str(lastsnap)
    halo_exists = os.path.exists(outpath+'/halos/halos_'+snapstr+'/halos_'+snapstr+'.0.bin')
    if not particles:
        return halo_exists
    part_exists = os.path.exists(outpath+'/halos/halos_'+snapstr+'/halos_'+snapstr+'.0.particles')
    return halo_exists and part_exists

def check_mergertree_exists(outpath,autoconvert=False):
    ascii_exists = os.path.exists(outpath+'/trees/tree_0_0_0.dat')
    binary_exists = os.path.exists(outpath+'/trees/tree.bin')
    if not binary_exists and autoconvert:
        print "---check_mergertree_exists: Automatically converting ascii to binary"
        MTC.convertmt(outpath+'/trees',version=3)
        binary_exists = os.path.exists(outpath+'/trees/tree.bin')
    return ascii_exists and binary_exists

def find_halo_paths(basepath="/bigbang/data/AnnaGroup/caterpillar/halos",
                    nrvirlist=[3,4,5,6],levellist=[11,12,13,14],ictype="BB",
                    require_rockstar=False,
                    require_mergertree=False,autoconvert_mergertree=False,
                    onlychecklastsnap=False,verbose=False,hdf5=True):
    """ Returns a list of paths to halos that have gadget completed/rsynced
        with the specified nrvirlist/levellist/ictype """
    if verbose:
        print "basepath:",basepath
        print "nrvirlist:",nrvirlist
        print "levellist:",levellist
        print "ictype:",ictype
    def gadget_finished(outpath):
        numsnaps = get_numsnaps(outpath)
        gadgetpath = outpath+'/outputs'
        if (not os.path.exists(gadgetpath)):
            if verbose: print "  Gadget folder not present in "+get_foldername(outpath)
            return False
        if onlychecklastsnap: #only check last snap
            snapstr = str(numsnaps-1).zfill(3)
            snappath = gadgetpath+"/snapdir_"+snapstr+"/snap_"+snapstr+".0"
            if hdf5: snappath += ".hdf5"
            if (not os.path.exists(snappath)):
                if verbose: print "  Snap "+snapstr+" not in "+get_foldername(outpath)
                return False
            else:
                return True
        for snap in xrange(numsnaps): # check that all snaps are there
            snapstr = str(snap).zfill(3)
            snappath = gadgetpath+"/snapdir_"+snapstr+"/snap_"+snapstr+".0"
            if hdf5: snappath += ".hdf5"
            if (not os.path.exists(snappath)):
                if verbose: print "  Snap "+snapstr+" not in "+get_foldername(outpath)
                return False
        return True

    halopathlist = []
    haloidlist = []
    for filename in os.listdir(basepath):
        if filename[0] == "H":
            haloidlist.append(filename)
    for haloid in haloidlist:
        subdirnames = basepath + "/" + haloid
        halosubdirlist = []
        try:
            for filename in os.listdir(subdirnames):
                halosubdirlist.append(filename)
                fileparts =  filename.split("_")
                levelmax = float(fileparts[5][2:4])
                nrvir = fileparts[7][-1]
                haloid = fileparts[0]
                if (int(levelmax) in levellist and int(nrvir) in nrvirlist and fileparts[1]==ictype):
                    outpath = basepath+"/"+haloid+"/"+filename
                    if gadget_finished(outpath):
                        halopathlist.append(outpath)
        except:
            continue

    if require_rockstar:
        newhalopathlist = []
        for outpath in halopathlist:
            if check_last_rockstar_exists(outpath):
                newhalopathlist.append(outpath) 
        halopathlist = newhalopathlist
    if require_mergertree:
        newhalopathlist = []
        for outpath in halopathlist:
            if check_mergertree_exists(outpath,autoconvert=autoconvert_mergertree):
                newhalopathlist.append(outpath) 
        halopathlist = newhalopathlist
    return halopathlist
