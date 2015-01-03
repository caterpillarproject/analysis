
import numpy as np
import pylab as plt
import sys
from random import randint
from colorsys import hsv_to_rgb
from random import randint, uniform
from matplotlib import *
import glob as glob

import readsnapshots.readsnapHDF5 as rsHD
import readsnapshots.readsnap as rs
import readhalos.readsubf
import readsnapshots.readids
import mergertrees.MTCatalogue as MT

from brendanlib.grifflib import convert_pid_zid,create_mt_image
import alexlib.haloutils as htils


"""
Contact: Brendan Griffen <brendan.f.griffen@gmail.com>

This script will print out two images relating to the merger 
tree of the caterpillar halos *from the parent*.

IMAGE 1:
H{haloid}_all_mt_properties.png
This will be a panel of 15 images of every quantity 
rockstar provides as a function of scale factor. 

IMAGE 2:
H{haloid}_stacked_mvir_vmax.png
This will plot a stacked image of all the mvir and vmax
evolution of all of the host halos in one plot. This helps
observe the scatter in their evolutionary histories.

"""

WANT_ZOOM = True

# THIS NEEDS TO BE MODIFIED TO RECEIVE THE NEW CANDIDATE LIST
base_halo_path = "/bigbang/data/AnnaGroup/caterpillar/halos/"

# PLOTTING SETTINGS
ticksize = 11
#ncols = 0

# SET PATHS
image_path = '/bigbang/data/bgriffen/caterpillar/zoom_mts/'
ext = "/outputs/snapdir_000/snap_000"
extics = 'ics/ics'

haloidlist_full = glob.glob(base_halo_path + "H*/H*LX14*")

if WANT_ZOOM:
    for hpath in haloidlist_full:
        if os.path.isdir(hpath+"/halos/"):
            pid = htils.get_parent_hid(hpath.split("/")[-1])
            zoomid = convert_pid_zid(pid,14)
            header = rsHD.snapshot_header(hpath + ext + ".0.hdf5")
            create_mt_image(hpath+"/halos/",header,zoomid,pid,image_path)

else:
    hpath = "/bigbang/data/AnnaGroup/caterpillar/parent/gL100X10/"
    candidatelist = []

    for haloid in haloidlist_full:
        candidatelist.append(int(haloid.split("/H")[-2]))

    for haloid in candidatelist:
        header = rsHD.snapshot_header(basepath + ext + ".0.hdf5")
        create_mt_image(hpath,header,haloid,haloid,image_path)
