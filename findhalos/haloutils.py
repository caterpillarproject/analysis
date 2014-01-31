import numpy as np
import os
import sys
import subprocess
from optparse import OptionParser

def get_numsnaps(outpath):
    return sum(1 for line in open(outpath+'/ExpansionList'))

def get_foldername(outpath):
    return os.path.basename(os.path.normpath(outpath))

def check_last_subfind_exists(outpath):
    numsnaps = get_numsnaps(outpath)
    lastsnap = numsnaps - 1; snapstr = str(lastsnap).zfill(3)
    group_tab = os.path.exists(outpath+'/outputs/groups_'+snapstr+'/group_tab_'+snapstr+'.0')
    subhalo_tab = os.path.exists(outpath+'/outputs/groups_'+snapstr+'/subhalo_tab_'+snapstr+'.0')
    return group_tab and subhalo_tab

def check_last_rockstar_exists(outpath):
    numsnaps = get_numsnaps(outpath)
    lastsnap = numsnaps - 1; snapstr = str(lastsnap)
    halo_exists = os.path.exists(outpath+'/halos/halos_'+snapstr+'/halos_'+snapstr+'.0.bin')
    part_exists = os.path.exists(outpath+'/halos/halos_'+snapstr+'/halos_'+snapstr+'.0.particles')
    return halo_exists and part_exists

def find_halo_paths(basepath="/bigbang/data/AnnaGroup/caterpillar/halos",
                    nrvirlist=[3,4,5,6],levellist=[11,12,13,14],ictype="BB",
                    onlychecklastsnap=False,verbose=False):
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
            if (not os.path.exists(gadgetpath+"/snapdir_"+snapstr+"/snap_"+snapstr+".0")):
                if verbose: print "  Snap "+snapstr+" not in "+get_foldername(outpath)
                return False
            else:
                return True
        for snap in xrange(numsnaps): # check that all snaps are there
            snapstr = str(snap).zfill(3)
            if (not os.path.exists(gadgetpath+"/snapdir_"+snapstr+"/snap_"+snapstr+".0")):
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
    return halopathlist
