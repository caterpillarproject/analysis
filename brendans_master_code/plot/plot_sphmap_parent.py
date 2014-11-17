import matplotlib.pyplot as plt
import numpy as np
import pylab,sys,glob,os

import calc.CalcHsml as ca
import viz.SphMap as SphMap
import readsnapshots.readhsml as readhsml
import readhalos.readgroup as readgroup
import readhalos.readsubf as readsubf
import readsnapshots.readsnapHDF5_greg as rs

from brendanlib.grifflib import makecolormap
#import asciitable

base_path = "/bigbang/data/AnnaGroup/caterpillar/parent/gL100X10"

#list of snapshots
snapnums=np.array(np.arange(127))

#map resolution
nx=512
ny=512
numngb=4

#x and y directions (projected along third)
dim0=0
dim1=1

mode = 3
periodic_bool = 1

#mode
# 1 "COLUMN MASS MAP
# 2 "QUANTITY MASS-WEIGHTED MAP
# 3 "COLUMN DENSITY MAP

boxsize = 100.
box=np.array([boxsize,boxsize,boxsize], dtype="float32")

totlen = 0
hubble = 0.6711
center=np.array([boxsize/2.,boxsize/2.,boxsize/2.], dtype="float32")
cmap = makecolormap()
print 'Processing ', base_path
for num in snapnums:
    if not os.path.isfile("./data/parent_proj_dim"+str(nx)+"_snap"+str(num).zfill(3)+".dat"):
        print "Doing: snapdir_" + str(num).zfill(3) + "/snap_" +str(num).zfill(3)

        snapbase = "/outputs/snapdir_" + str(num).zfill(3) + "/snap_" +str(num).zfill(3)
        
        #
        pos = rs.read_block(base_path + snapbase, "POS ",parttype=1).astype("float32")
        mass = rs.read_block(base_path + snapbase, "MASS",parttype=1).astype("float32")
    
        hsml = readhsml.hsml_file(base_path+"/outputs/",num)
        hsmlvals = hsml.Hsml
    
        boxmass=mass[(abs(pos[:,0]-center[0]) < box[0]/2.) & (abs(pos[:,1]-center[1]) < box[1]/2.) & (abs(pos[:,2]-center[2]) < box[2]/2.)].sum()
        map = SphMap.CalcDensProjection(pos, hsmlvals, mass, mass, nx, ny, box[0], box[1], box[2], center[0], center[1], center[2], dim0, dim1, mode, periodic_bool)
        
        filename_out="./data/parent_proj_dim"+str(nx)+"_snap"+str(num).zfill(3)+".dat"
        #f = open(filename_out,'wb')
        #np.array([nx],dtype=np.uint32).tofile(f)
        #np.array([ny],dtype=np.uint32).tofile(f)
        #map.astype("float32").tofile(f)
        np.savetxt(filename_out,np.array(map))
        #f.close()
        #make image
 #       map = map * 10.0**10.0
    
 #       fig = plt.figure( figsize = (10.,10.) )
        #pylab.spectral()
 #       ma=map.max()/2.
 #       mi=ma/10000.0
        
 #       ma=np.log10(ma)
 #       mi=np.log10(mi) 
#        print 'mean density of map = ', map.mean()
#        map=np.log10(map)       
    
#        print map.min(), map.max()
#        print mi,ma
            
 #       map[map<mi]=mi
#        map[map>ma]=ma
            
#        print 'min/max map=', map.min(), map.max()  
#        cmap = makecolormap()
#        plt.imshow(np.array(map).T,origin='lower',cmap= cmap,interpolation='nearest',extent=[center[0]-box[0]/2., center[0]+box[0]/2., center[1]-box[1]/2., center[1]+box[1]/2.])#,vmin=np.log10(min_dens), vmax=np.log10(max_dens))
#    
#        plt.axis([center[0]-box[0]/2., center[0]+box[0]/2., center[1]-box[1]/2., center[1]+box[1]/2.])
#        plt.xlabel('Mpc/$h$', fontsize=20)
#        plt.ylabel('Mpc/$h$', fontsize=20)  
#    
#        plt.savefig("./imgs/PARENT_L"+str(int(boxsize))+"MPC_SNAP"+str(num).zfill(3)+".eps",bbox_inches="tight")
#        plt.savefig("./imgs/PARENT_L"+str(int(boxsize))+"MPC_SNAP"+str(num).zfill(3)+".png",bbox_inches="tight")
#        plt.close(fig)
    #plt.show()
