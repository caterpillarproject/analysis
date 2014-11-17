import matplotlib.pyplot as plt
import numpy as np
import pylab,sys,glob,os
from matplotlib.ticker import NullFormatter

import calc.CalcHsml as ca
import viz.SphMap as SphMap
import readsnapshots.readhsml as readhsml
import readhalos.readgroup as readgroup
import readhalos.readsubf as readsubf
import readsnapshots.readsnapHDF5_greg as rs
import mergertrees.MTCatalogue as MT

from brendanlib.grifflib import makecolormap
import brendanlib.constants as const
import asciitable

WANT_ALL_PARTICLE_TYPES = True
base_path = "/bigbang/data/AnnaGroup/caterpillar/halos/"

#halosindex = asciitable.read("/bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index_nofix.txt", Reader=asciitable.FixedWidth)
halosindex = asciitable.read("/bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index.txt", Reader=asciitable.FixedWidth)

halo_dirs = glob.glob(base_path+"H*/H*LX13*")

#list of snapshots
snapnums=np.array([255])

#map resolution
nx=512
ny=512
numngb=4

#x and y directions (projected along third)
dim0=0
dim1=1

mode = 3
periodic = False

if periodic:
    periodic_bool = 1
else:
    periodic_bool = 0

#mode
# 1 "COLUMN MASS MAP
# 2 "QUANTITY MASS-WEIGHTED MAP
# 3 "COLUMN DENSITY MAP

boxsize = 1.
box=np.array([boxsize,boxsize,boxsize], dtype="float32")


def makeimage(map,basepath,filename,verbose=True):
    #make image
    map = map * 10.0**10.0
    #pylab.spectral()
    ma=map.max()/2.
    mi=ma/10000.0
    
    ma=np.log10(ma)
    mi=np.log10(mi) 
    
    map=np.log10(map)
    
    if verbose:
        print 'mean density of map = ', map.mean()
        print map.min(), map.max()
        print mi,ma
        
    map[map<mi]=mi
    map[map>ma]=ma
        
    print 'min/max map=', map.min(), map.max()  
    cmap = makecolormap()

    fig = plt.figure( figsize = (10.,10.) )
    ax = fig.add_subplot(111)
    ax.imshow(np.array(map).T,origin='lower',cmap=cmap, interpolation='nearest',extent=[center0[0]-box[0]/2., center0[0]+box[0]/2., center0[1]-box[1]/2., center0[1]+box[1]/2.])#,vmin=np.log10(min_dens), vmax=np.log10(max_dens))
    ax.axis([center0[0]-box[0]/2., center0[0]+box[0]/2., center0[1]-box[1]/2., center0[1]+box[1]/2.])
    ax.set_xlabel('Mpc', fontsize=20)
    ax.set_ylabel('Mpc', fontsize=20)  
    plt.savefig(basepath+"/eps/"+filename+".eps",bbox_inches="tight")
    plt.savefig(basepath+"/png/"+filename+".png",bbox_inches="tight")
    plt.close(fig)

    fig = plt.figure( figsize = (10.,10.) )
    ax = fig.add_subplot(111)
    ax.imshow(np.array(map).T,origin='lower',cmap=cmap, interpolation='nearest',extent=[center0[0]-box[0]/2., center0[0]+box[0]/2., center0[1]-box[1]/2., center0[1]+box[1]/2.])#,vmin=np.log10(min_dens), vmax=np.log10(max_dens))
    ax.axis([center0[0]-box[0]/2., center0[0]+box[0]/2., center0[1]-box[1]/2., center0[1]+box[1]/2.])
    nullfmt = NullFormatter()
    ax.xaxis.set_major_formatter(nullfmt)
    ax.yaxis.set_major_formatter(nullfmt)
    plt.savefig(basepath+"/eps/"+filename+"_NOTICK.eps",bbox_inches="tight")
    plt.savefig(basepath+"/png/"+filename+"_NOTICK.png",bbox_inches="tight")
    plt.close(fig)

def writemap(filename,nx,ny):
    f = open(filename_out,'wb')
    np.array([nx],dtype=np.uint32).tofile(f)
    np.array([ny],dtype=np.uint32).tofile(f)
    map.astype("float32").tofile(f)
    f.close()

for halofullpath in halo_dirs:
    if os.path.isdir(halofullpath + "/halos/halos_255/") and os.path.isfile(halofullpath+"/halos/trees/tree.bin") and os.path.isdir(halofullpath + "/outputs/hsmldir_255/"):
        foldername = halofullpath.split("/")[-1]
        LX = int(foldername.split("_")[5][-2:])
        NV = int(foldername.split("_")[-1][-1])
        htype = foldername.split("_")[1]
        halo = foldername.split("_")[0]
        
        mask = ((int(halo[1:]) == halosindex['parentid']) & (htype == halosindex['ictype']) & (NV == halosindex['NV']) & (LX == halosindex['LX']))
        
        if mask.any():
            parentid = halosindex['parentid'][mask]
            zoomid = int(halosindex['zoomid'][mask])
            cat = MT.MTCatalogue(halofullpath + '/halos/trees',indexbyrsid=False,haloids=[zoomid])
            center_x = cat[0].data['posX'][0]
            center_y = cat[0].data['posY'][0]
            center_z = cat[0].data['posZ'][0]
            center=np.array([center_x,center_y,center_z], dtype="float32")

            for num in snapnums:
                snapbase = "/outputs/snapdir_" + str(num).zfill(3) + "/snap_" +str(num).zfill(3)
                print
                print 'Processing ', halofullpath.split("/")[-1]
                pos = rs.read_block(halofullpath + snapbase, "POS ",parttype=1).astype("float32")
                mass = rs.read_block(halofullpath + snapbase, "MASS",parttype=1).astype("float32")
            
                hsml = readhsml.hsml_file(halofullpath+"/outputs/",255)
                hsmlvals = hsml.Hsml

                if WANT_ALL_PARTICLE_TYPES:
                    for typei in range(2,6):
                        print
                        print "Adding particle type:",typei
                        pos_tmp = rs.read_block(halofullpath + snapbase, "POS ",parttype=typei).astype("float32")
                        mass_tmp = rs.read_block(halofullpath + snapbase, "MASS",parttype=typei).astype("float32")
                        hsml_tmp=ca.CalcHsml(pos_tmp, mass_tmp, float(numngb), float(0), float(0.0), float(1.0), float(1.0), float(0.0))
                        hsmlvals = np.hstack([hsmlvals,hsml_tmp])
                        
                    pos = rs.read_block(halofullpath + snapbase, "POS ").astype("float32")
                    mass = rs.read_block(halofullpath + snapbase, "MASS").astype("float32")

                pos = (pos - center)/const.h
                center0 = np.array([0.,0.,0.], dtype="float32")
                boxmass=mass[(abs(pos[:,0]-center0[0]) < box[0]/2.) & (abs(pos[:,1]-center0[1]) < box[1]/2.) & (abs(pos[:,2]-center0[2]) < box[2]/2.)].sum()
                map = SphMap.CalcDensProjection(pos, hsmlvals, mass, mass, nx, ny, box[0], box[1], box[2], center0[0], center0[1], center0[2], dim0, dim1, mode, periodic_bool)

                base_out = '/bigbang/data/bgriffen/caterpillar/densityproj/'
                filename_out="H"+str(parentid[0])+"_proj_density_field.dat"
                writemap(base_out+filename_out,nx,ny)
                filename_img = "H"+str(parentid[0])+"_DPROJ_LX"+str(LX)+"_L"+str(int(boxsize))+"MPC_SNAP"+str(num).zfill(3)
                makeimage(map,base_out,filename_img)



                    #hsmlvals = hsml.Hsml[0:totlen]
                    #hsml = np.array([1.25]*len(pos))
                    #print 'HSML values:'
                    #print hsmlvals
                    #sys.exit()
                    

                    #pos=np.vstack([pos,rs.read_block(filename, "POS ", parttype=1,doubleprec=False)])
                    #mass=np.hstack([mass,rs.read_block(filename, "MASS", parttype=1,doubleprec=False)])
                
                    #pos=rs.read_block(filename, "POS ", parttype=1,doubleprec=False)

                                        #print map[map == 0]
                    #save data
                    #plt.clf()
                    #pos -= cat.GroupCM[0,:]
                    #    totlen += len(pos)
                    #    print totlen
                    #mass=rs.read_block(filename, "MASS", parttype=1,doubleprec=False)
                    #u=np.float64(rs.read_block(filename, "U   ", parttype=1,doubleprec=False))
                    #ne=np.float64(rs.read_block(filename, "NE  ", parttype=1,doubleprec=False))
                    #temp=np.float32(co.GetTemp(u, ne, 5./3.))
                    #hsml=1.25*rs.read_block(filename, "HSML", parttype=1,doubleprec=False)
                    
                    #hsml=ca.CalcHsml(pos, mass, float(numngb), float(0), float(0.0), float(1.0), float(1.0), float(0.0))
                    
                    #sys.exit()
                
