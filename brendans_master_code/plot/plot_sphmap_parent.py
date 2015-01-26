import matplotlib.pyplot as plt
import numpy as np
import pylab, glob, os
import subprocess as sub

import haloutils as htils
import ctils
import asciitable

base_path = "/bigbang/data/AnnaGroup/caterpillar/parent/gL100X9"

ndim=256
snapshot = 63

def run_algo(hpath):
    img_path = "/bigbang/data/AnnaGroup/Dropbox/caterpillar_plots/zoom/density_projections/static/png/"
    snapbase = "/outputs/snapdir_063/"
    
    print "Generating hsml values..."
    hsmlpath = hpath+"/analysis/hsml_ngb4.npy"
    if not os.path.isfile(hsmlpath) and os.path.isdir(hpath+snapbase): ctils.generate_hsml(hpath,snapnum=snapshot,ndim=256)

    print "Generating map..."
    map_path = hpath + "/analysis/projected_xy_density_field.npy"
    if not os.path.isfile(map_path) and os.path.isfile(hsmlpath): ctils.construct_map_parent(hpath,ndim=256)

    print "Generating image..."
    image_filename = img_path+"PARENT_DPROJ_DIM"+str(ndim)+"_NGB4_SNAP"+str(127).zfill(3)+".png"
    if not os.path.isfile(image_filename) and os.path.isfile(map_path): ctils.generate_img_parent(hpath,image_filename)  

run_algo(base_path)
