import matplotlib.pyplot as plt
import numpy as np
import pylab, glob, os
import subprocess as sub

import haloutils as htils
import ctils
import asciitable

base_path = "/bigbang/data/AnnaGroup/caterpillar/halos/"
halo_dirs = glob.glob(base_path+"H*/H*LX14*")

ndim=2048

def run_algo(hpath,ndim):
    img_path = "/bigbang/data/AnnaGroup/Dropbox/caterpillar_plots/zoom/density_projections/static/png/"
    
    ictype, LX, NV = htils.get_zoom_params(hpath)
    parentid = htils.get_parent_hid(hpath)
    #rvir = htils.get_quant_zoom(hpath,"rvmax")

    print
    print 'Processing:', hpath.split("/")[-1]
    
    cmd_make = "mkdir -p " + hpath + "/analysis/"
    sub.call([cmd_make],shell=True)

    snapbase = "/outputs/snapdir_255/"
    try: 

        hsmlpath = hpath+"/analysis/hsml_ngb4.npy"
        if not os.path.isfile(hsmlpath) and os.path.isdir(hpath+snapbase): ctils.generate_hsml(hpath)

        map_path = hpath + "/analysis/projected_xy_density_field.npy"
        if not os.path.isfile(map_path) and os.path.isfile(hsmlpath): ctils.construct_map(hpath)

        image_filename = img_path+"H"+str(parentid)+"_DPROJ_LX"+str(LX)+"_DIM"+str(ndim)+"_NGB4_SNAP"+str(255).zfill(3)+".png"
        if not os.path.isfile(image_filename) and os.path.isfile(map_path): ctils.generate_img(hpath,image_filename)  
    
    except ValueError:
        print "H"+str(parentid)+" not in index!"

#import multiprocessing as mp

#pool = mp.Pool(processes=4)
#results = [pool.apply_async(run_algo, args=(hpath,ndim)) for hpath in halo_dirs]

#print(results)
#[1, 8, 27, 64, 125, 216]
#pool = mp.Pool(processes=4)
#
#print(results)

for hpath in halo_dirs:
    run_algo(hpath,ndim)
