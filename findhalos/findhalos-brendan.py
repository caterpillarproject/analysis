#IMPORT CORE MODULES
import numpy as np
import pylab as plt
import sys, os, platform,glob
import matplotlib
#IMPORT CATERPILLAR MODULES
import readhalos.readsubf as readsubf
import readhalos.readgroup as readgroup
import readhalos.RSDataReaderv2 as RSDataReader
import readsnapshots.readsnap as rs
import readsnapshots.readsnapHDF5 as rhdf5

#IMPORT PERSONAL LIBRARIS
from brendanlib.grifflib import *

#IMPORTANT STUFF
hubble = 0.6711
snapnum = 255

basepath = determinebasepath(platform.node())
listdir = glob.glob(basepath + "caterpillar/halos/H*/*_BB_*_NV3/outputs/groups_255")

fig = plt.figure(figsize=(9,9))
axa = fig.add_subplot(221)
axb = fig.add_subplot(222)
axc = fig.add_subplot(223)
axd = fig.add_subplot(224)

listdir = glob.glob(basepath + "caterpillar/halos/H*/")

nvlist = ["NV3"]
halotypelist = ["BB"]
levelmaxlist = ["LX11","LX12","LX13","LX14"]

halolist = ["H121869","H140666","H147781","H173411","H1930","H21047",
"H233776","H241932","H242183","H260266","H268422","H69786","H96624"]

halolist = ["H268422"]

for halo in listdir:
    masslist = []
    vmaxlist = []
    posxlist = []
    posylist = []
    npfoflist = []
    npsublist = []


    haloname = halo.split("/")[6]
    print "----------------------------------------------"
    for haloin in glob.glob(halo+"/H*_BB_*"):
   # for haloin in 
        
        haloparts = haloin.split("/")[7].split("_")
        levelmax = haloparts[5]
        halotype = haloparts[1]
        nv = haloparts[7]
        if nv in nvlist and halotype in halotypelist and levelmax in levelmaxlist and haloname in halolist:
            outputsdir = haloin+"/outputs"

            if os.path.isdir(outputsdir+"/groups_255/"):
#               header = rs.snapshot_header(outputsdir+"/snapdir_255/snap_255.0")
                s = readsubf.subfind_catalog(outputsdir, snapnum)
                id = readsubf.subf_ids(outputsdir, snapnum, 0, 0, read_all=1)
                
                submass = s.sub_mass*10**10/hubble
                subposx = s.sub_pos[:,0]
                subposy = s.sub_pos[:,1]
                subposz = s.sub_pos[:,2]
                
                groupx = s.group_pos[0,0]
                groupy = s.group_pos[0,1]
                groupz = s.group_pos[0,2]
                
                R = np.sqrt((groupx-subposx)**2 + (groupy-subposy)**2 + (groupz-subposz)**2)
                
                mask = (R < 0.5) & (submass <= 5e12) & (submass >= 7e11)

                if not mask.any():
                    mask = (R < 0.5) & (submass <= 9e12) & (submass >= 7e11)

                level = int(levelmax[2:])
        
                if level == 11:
                    style = 'ro'
                elif level == 12:
                    style = 'go'
                elif level == 13:
                    style = 'bo'
                elif level == 14:
                    style = 'ko'

                if mask.any():
                    masslist.append(submass[mask])
                    vmaxlist.append(s.sub_vmax[mask])
                    posxlist.append(subposx[mask])
                    posylist.append(subposy[mask])
                    npfoflist.append(s.group_len[0])
                    npsublist.append(np.log10(s.sub_len[0]))
                    axa.plot(posxlist,posylist,style,markeredgewidth=0)
                    axc.plot(npfoflist,npsublist,style,markeredgewidth=0)
                    axb.plot(masslist,vmaxlist,style,markeredgewidth=0)

                    #print ""
                    print "++",haloin.replace("/bigbang/data/AnnaGroup/caterpillar/halos/","")
                    print "+",levelmax,nv
                    #haloin.replace(basepath+"caterpillar/halos/","")
                    print "MASS: %3.2e" % (submass[mask])
                    print "VMAX: %3.2f" % (s.sub_vmax[mask])
                    print "SUB-NPART: %i6" % (s.sub_len[mask])
                    print "FOF-NPART:%i6" % (s.group_len[0])
                    print "X-Y-Z: %3.2f,%3.2f,%3.2f" % (subposx[mask],subposy[mask],subposz[mask])

#                    print "mp: %3.2e" % header.massarr[1]*10**10

    masslist = np.array(masslist)
    vmaxlist = np.array(vmaxlist)
    posxlist = np.array(posxlist)
    posylist = np.array(posylist)
    npfoflist = np.array(npfoflist)
    npsublist = np.array(npsublist)

    axa.plot(posxlist,posylist,'b-')
    axc.plot(npfoflist,npsublist,'b-')
    axb.plot(masslist,vmaxlist,'b-')

axa.set_xlabel(r'x-pos [Mpc/h]')
axa.set_ylabel(r'y-pos [Mpc/h]')
axb.set_ylabel(r'vmax [km/s]')
axb.set_xlabel(r'mass [M$_\odot$]')
axc.set_ylabel(r'log$_{10}$ [# particles in subhalo]')
axc.set_xlabel(r'log$_{10}$ [# particles in FOF group]')
plt.show()
