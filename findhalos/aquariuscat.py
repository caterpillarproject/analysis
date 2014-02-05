import numpy as np
import matplotlib.pyplot as plt
import pynbody as pnb
from brendanlib.grifflib import makecolormap
import readhalos.readgroup as RG
import readhalos.readsubf as readsubf
import readsnapshots.readsnap as rs
from operator import itemgetter
import glob as glob
import os

hubble = 0.6711
WANTHEATMAP =  False
WANTRADIAL = True
FULLSTDOUT = True
fofrcut = 100.5 # Mpc/h

basepath = "/bigbang/data/AnnaGroup/caterpillar/halos/"

def readblock(path,parttype):
    pos = rs.read_block(path,"POS ",parttype=parttype,doubleprec=True,verbose=False)
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]
    #print "-- number of particles:",len(x)
    return x,y,z

def getmaxsnapnum(halo):
    if halo == "Aq-A":
        snap = 1023
        swapval = False
        long_ids = True
    elif halo == "Aq-F":
        snap = 111
        swapval = False
        long_ids = False
    elif halo == "Aq-E":
        snap = 127
        swapval = True
        long_ids = False
    else:
        snap = 127
        swapval = False
        long_ids = True

    return snap,swapval,long_ids

halolist = ["Aq-A","Aq-B","Aq-C","Aq-D","Aq-E","Aq-F"]

fig = plt.figure()
ax = fig.add_subplot(111)
colorwheel = ['r','b','c','m','g','k']
for haloi,colori in zip(halolist,colorwheel):
    maxsnap,swapval,long_ids = getmaxsnapnum(haloi)
    maxsnap = str(maxsnap)
    path = "/bigbang/data/AnnaGroup/aquarius/data/" + haloi + "/2/"
    numblocks = []
    #path = haloi
    if os.path.isdir(path+'groups_'+maxsnap):
         #and "LX14" in haloi
        
        s = readsubf.subfind_catalog(path, maxsnap,long_ids = long_ids, swap = swapval)
        contamNR = s.group_contamination_count
        contamMR = s.group_contamination_mass*10**10/hubble
        minindex = min(enumerate(s.group_contamination_count[0:3]), key=itemgetter(1))[0]
        pcmass_contam = (s.group_contamination_mass/s.group_mass)*100.
        allgroupmass = s.group_mass*10**10/hubble
        # GROUPS
        candgrouppos = s.group_pos[minindex]
        candgroupmass = allgroupmass[minindex]
        # SUBS
        maxsubindex = np.argsort(s.sub_len)
        submass = s.sub_mass[maxsubindex][-1]*10**10/hubble
        subpos = s.sub_pos[maxsubindex][-1]
        
        if FULLSTDOUT:
            print "--------------------------------------------------------"
            print haloi.replace(basepath,"")
            print "FOF With Smallest Contam. Position:",candgrouppos
            print "FOF With Smallest Contam. Mass: %0.2e" % (candgroupmass)
            print "Subhalo Mass With Most Particles: %0.2e" % (submass)
            print "Subhalo Position With Most Particles:",subpos
            print "Contamination Num.:",contamNR[minindex]
            print "Contamination Mass: %0.2e [%3.3f pc]" % (contamMR[minindex]/0.6711,(contamMR[minindex]/candgroupmass)*100)
        
        if contamMR[minindex]/0.6711 > 0 and not FULLSTDOUT:
            print "--------------------------------------------------------"
            print haloi.replace(basepath,"")
            print "FOF With Smallest Contam. Mass: %0.2e" % (candgroupmass)
            print "Contamination Mass: %0.2e [%3.3f pc]" % (contamMR[minindex]/0.6711,(contamMR[minindex]/candgroupmass)*100)
        
        rfof = np.sqrt((candgrouppos[0]-s.group_pos[:,0])**2 + (candgrouppos[1]-s.group_pos[:,1])**2 + (candgrouppos[2]-s.group_pos[:,2])**2)
            
        if WANTRADIAL:
            idxrfof = np.argsort(rfof)
            contam_n_radial = np.cumsum(contamNR[idxrfof])
            contam_m_radial = np.cumsum(contamMR[idxrfof])

            ax.plot(rfof[idxrfof],np.log10(contam_n_radial),linestyle='-',color=colori,linewidth=2,label=haloi)
            ax.set_xlim([0,7])
            ax.set_ylim([0,5])
            ax.set_ylabel(r'$\mathrm{\Sigma\ log_{10}\ N_{CP}(<R)}$', color=colori,fontsize=14)

            # set dual side axes for both mass and number contamination count
            for tl in ax.get_yticklabels():
                tl.set_color('b')
            axb = ax.twinx()
            axb.tick_params(axis='both', which='major', labelsize=12)
            axb.set_ylim([8,14])
            axb.plot(rfof[idxrfof],np.log10(contam_m_radial),linestyle='--',color=colori,linewidth=2)
            axb.set_ylim([8,14])
            for tl in axb.get_yticklabels():
                tl.set_color('r')
            axb.set_ylabel(r'$\mathrm{log_{10}\ \Sigma\ M_{CP}(<R)\ [M_\odot/h]}$', color=colori,fontsize=14)
            ax.set_xlabel(r'$\mathrm{R_{FOF}\ [Mpc/h]}$',fontsize=14)
            ax.set_xlim([0,5])
            #ax.set_title(haloi)
            
        # create heatmap for % of FOF mass contaminated versus FOF group mass.
        if WANTHEATMAP:
            mask = np.where(rfof<=fofrcut)
            xvar = np.log10(s.group_mass[mask]*10**10/hubble)
            yvar = pcmass_contam[mask]
            heatmap, xedges, yedges = np.histogram2d(xvar,yvar, bins=128)
            fig = plt.figure()
            ax = fig.add_subplot(111)
            if len(mask) != len(allgroupmass):
                sc1 = ax.scatter(xvar,yvar,c=rfof[mask],marker='o', edgecolors='none')
                cbar1 = fig.colorbar(sc1)
                cbar1.set_label('FOF distance')
            else:
                extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
                heatmap = np.flipud(np.rot90(heatmap))
                sc1 = ax.imshow(np.log10(heatmap),extent = extent,cmap = 'jet', origin='lower')
                cbar1 = fig.colorbar(sc1)
                cbar1.set_label('# FOF Groups')
                ax.set_aspect('auto')
            ax.set_title(haloi)
            ax.set_ylabel('% mass contaminated')
            ax.set_xlabel('log [FOF mass]')
            
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels,loc='lower right')
        #if os.path.isdir(path+'snapdir_'+maxsnap):
        #    for parttype in xrange(2,6):
        #            x,y,z = readblock(path+'snapdir_'+maxsnap+'/snap_'+maxsnap,parttype=parttype)
        #            dr = np.sqrt((subpos[0]-x)**2 + (subpos[1]-y)**2 + (subpos[2]-z)**2)
        #            if FULLSTDOUT:
        #                print "-- Particle Type %i is %3.1f kpc away from subhalo" % (parttype,min(dr)*1000.)

plt.show()