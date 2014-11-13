import readhalos.RSDataReader as RSDataReader
import readsnapshots.readsnapHDF5_greg as rsHD
from brendanlib.grifflib import CorrectPos,COM
import readsnapshots.readsnap as rs
import glob
import numpy as np
import os


"""

Contact: Brendan Griffen <brendan.f.griffen@gmail.com>

This script will construct the lagrangian files for a given 
set of halo IDs set by what exists inside halos/.

"""

writelagrfile = True

base_halo_path = "/bigbang/data/AnnaGroup/caterpillar/halos/"
basepath = "/bigbang/data/AnnaGroup/caterpillar/parent/gL100X10/"
ext = "outputs/snapdir_127/snap_127"
extics = 'ics/ics'
halopath = basepath + 'rockstar'
lagroutputdir = "/bigbang/data/AnnaGroup/caterpillar/ics/lagr/"

# THIS NEEDS TO BE CONVERTED TO A LIST OF IDS NOT OF DIRECTORIES ALEX!
haloidlist_full = glob.glob(base_halo_path + "H*")

haloidlist = []
for haloid in haloidlist_full:
    haloidlist.append(int(haloid.split("/H")[1]))

nrvir_list = [4,5,6]

print "Reading halo catalogue from parent..."
halodata = RSDataReader.RSDataReader(halopath,127,version=6)

print "Getting host halos only..."
allhalos = halodata.get_hosts()

print "Reading z = 0 snapshot header..."
header=rsHD.snapshot_header(basepath+ext)

print "Reading z = 0 particle ids..."
snapIDS=rsHD.read_block(basepath+ext,"ID  ")

print "Reading z = 0 particle positions..."
snapPOS = rsHD.read_block(basepath+ext,"POS ")

print "Reading lagrangian particle ids..."
snapIDSlagr = rs.read_block(basepath+extics,"ID  ")

print "Reading lagrangian particle positions..."
snapPOSlagr = rs.read_block(basepath+extics,"POS ",doubleprec=False)

for haloid in haloidlist:
       
    print "Creating:",haloid
    rvircand = allhalos.ix[haloid]['rvir']
    mvircand = allhalos.ix[haloid]['mvir']
    posXcand = allhalos.ix[haloid]['posX']
    posYcand = allhalos.ix[haloid]['posY']
    posZcand = allhalos.ix[haloid]['posZ']
        
    print "------------------------------------------------"
    print "Rockstar ID inside parent simulation: ",haloid
    print "------------------------------------------------"
    print "            x-pos:",'{:.2f}'.format(float(allhalos.ix[haloid]['posX'])), "   \ [Mpc/h]"
    print "            y-pos:",'{:.2f}'.format(float(allhalos.ix[haloid]['posY'])), "   \ [Mpc/h]"
    print "            z-pos:",'{:.2f}'.format(float(allhalos.ix[haloid]['posZ'])), "   \ [Mpc/h]"
    print "      virial mass:",'{0:.2e}'.format(float(allhalos.ix[haloid]['mvir'])/header.hubble),"\ [Msol]"
    print "    virial radius:",'{:.2f}'.format(float(allhalos.ix[haloid]['rvir'])),"  \ [kpc]"
    print "------------------------------------------------"
    
    for nrviri in nrvir_list:
        headerfilename = lagroutputdir + 'H' + str(haloid) + 'NRVIR' + str(int(nrviri)) + '.head'
            
        if not os.path.isfile(headerfilename):
            print headerfilename 
            print 'Constructing: H' + str(haloid) + 'NRVIR' + str(int(nrviri))
            nrviri = float(nrviri)

            dx = np.array(allhalos.ix[haloid]['posX']) - snapPOS[:,0]
            dy = np.array(allhalos.ix[haloid]['posY']) - snapPOS[:,1]
            dz = np.array(allhalos.ix[haloid]['posZ']) - snapPOS[:,2]
            R = np.sqrt(dx**2. + dy**2. + dz**2.)
            Rindex = np.where(R <= nrviri*rvircand/1000.)
            currentpos = snapPOS[Rindex[0]]
        
            regionIDS = snapIDS[Rindex[0]]
        
            mask = np.in1d(snapIDSlagr,regionIDS,assume_unique=True)
        
            lagrPos=snapPOSlagr[mask]
            del mask
        
            CorrectPos(lagrPos[:,0],header.boxsize)
            CorrectPos(lagrPos[:,1],header.boxsize)
            CorrectPos(lagrPos[:,2],header.boxsize)

            comV=COM(lagrPos[:,0],lagrPos[:,1],lagrPos[:,2])

            xmax = np.max(lagrPos[:,0]); xmin = np.min(lagrPos[:,0])
            dx = xmax-xmin
            centx = (xmax+xmin)/2.0
            ymax = np.max(lagrPos[:,1]); ymin = np.min(lagrPos[:,1])
            dy = ymax-ymin
            centy = (ymax+ymin)/2.0
            zmax = np.max(lagrPos[:,2]); zmin = np.min(lagrPos[:,2])
            dz = zmax-zmin
            centz = (zmax+zmin)/2.0


            centx = centx/header.boxsize
            centy = centy/header.boxsize
            centz = centz/header.boxsize

            extx = dx/header.boxsize
            exty = dy/header.boxsize
            extz = dz/header.boxsize
            
            #centx -= 
            #extx=2.0*dx/header.boxsize
            #exty=2.0*dy/header.boxsize
            #extz=2.0*dz/header.boxsize
        
            #extx=dx
            #exty=dy
            #extz=dz

            print "Writing lagrangian file: H" + str(haloid) + "NRVIR" + str(int(nrviri)) + ".head"
            if writelagrfile == True:
                
                # WRITE OUT BOX IC CENTERS AND EXTENTS
                headerfilename = lagroutputdir + '/H' + str(haloid) + 'NRVIR' + str(int(nrviri)) + '.head'
                f1=open(headerfilename,'w')
                f1.write('#' + str(centx) + '\n')
                f1.write('#' + str(centy) + '\n')
                f1.write('#' + str(centz) + '\n')
                f1.write('#' + str(extx) + '\n')
                f1.write('#' + str(exty) + '\n')
                f1.write('#' + str(extz) + '\n')
                f1.close()
                
                print centx,centy,centz,extx,exty,extz
                # WRITE OUT LAGRANGIAN REGION FILE
                filename = lagroutputdir + '/H' + str(haloid) + 'NRVIR' + str(int(nrviri))
                f2=open(filename,'w')
                for iv in xrange(0,len(lagrPos[:,0])):
                    f2.write(str(lagrPos[iv,0]/header.boxsize)+' '+str(lagrPos[iv,1]/header.boxsize)+' '+ str(lagrPos[iv,2]/header.boxsize)+'\n')                
                f2.close()
    
            print "> COMPLETED <"
            #sys.exit()
