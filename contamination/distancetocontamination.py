#IMPORT CORE MODULES
import numpy as np
import os, platform, asciitable

#IMPORT CATERPILLAR MODULES
import readhalos.RSDataReader as rstar
import readsnapshots.readsnapHDF5 as rhdf5

#IMPORT PERSONAL LIBRARIS
from brendanlib.grifflib import determinebasepath

basepath = determinebasepath(platform.node())

# SET PATHS
parentgadpath = basepath + 'caterpillar/parent/gL100X10/outputs'
parenthalopath = basepath + 'caterpillar/parent/gL100X10/rockstar/'
halopath = basepath + 'caterpillar/halos/'

parentid_index = []
ictype_index = []
LX_index = []
NV_index = []
zoomid_index = []
x_index = []
y_index = []
z_indez = []

halosindex = asciitable.read(halopath+"parent_zoom_index.txt", Reader=asciitable.FixedWidth)
haloparent = rstar.RSDataReader(parenthalopath,127,version=6)
allparenthosts = haloparent.get_hosts()

haloidlist = []
for filename in os.listdir(halopath):
    if filename[0] == "H":
        haloidlist.append(filename)

hubble = 0.6711
print "#folder name, distance to closest contam. particle [kpc/h], parent mass, zoom mass, ratio"
for halo in haloidlist:
    for folder in os.listdir(halopath+halo):
        halofullpath = halopath + halo + "/" + folder
        if os.path.isdir(halofullpath + "/outputs/snapdir_255/") and os.path.isdir(halofullpath + "/halos/halos_255/"):
            foldername = halofullpath.split("/")[-1]
            LX = int(foldername.split("_")[5][-2:])
            NV = int(foldername.split("_")[-1][-1])
            htype = foldername.split("_")[1]
            mask = ((int(halo[1:]) == halosindex['parentid']) & (htype == halosindex['ictype']) & (NV == halosindex['NV']) & (LX == halosindex['LX']))
            if mask.any():
                pos = rhdf5.read_block(halofullpath + "/outputs/snapdir_255/snap_255","POS ",parttype=2)
                host_xpos = halosindex['x'][mask]
                host_ypos = halosindex['y'][mask]
                host_zpos = halosindex['z'][mask]
                    
                dx = host_xpos - pos[:,0]
                dy = host_ypos - pos[:,1]
                dz = host_zpos - pos[:,2]
                R = np.sqrt(dx**2+dy**2+dz**2)*1000./hubble
                
                halodata = rstar.RSDataReader(halofullpath+"/halos",255,version=6)

                allhalos = halodata.get_hosts()
                zoomid = halosindex['zoomid'][mask]
                
                parentmass = float(allparenthosts.ix[halosindex['parentid'][mask]]['mvir']/hubble)
                zoommass = float(allhalos.ix[zoomid]['mvir']/hubble)
                ratiomass = zoommass/parentmass

                if np.abs(ratiomass-1) > 0.1:
                    print "%40s | %7.2f | %3.2e | %3.2e | %3.2f *" % (folder,R.min(),parentmass,zoommass,ratiomass)
                else:
                    print "%40s | %7.2f | %3.2e | %3.2e | %3.2f" % (folder,R.min(),parentmass,zoommass,ratiomass)