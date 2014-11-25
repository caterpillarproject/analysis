import asciitable
import glob
import numpy as np
import pylab as plt
import sys
import matplotlib
import matplotlib.patches as patches
from numpy import random

from brendanlib.grifflib import placetext_direct,makecolormap,placetext,plotxyzprojr,getlagrxyz,getcentext
import alexlib.haloutils as htils
import readhalos.RSDataReader as RSDataReader

"""
This code will plot 5 rows and 3 columns.
Each column is a different projections (x-y,x-z,y-z)

Each row is as follows:
1: Lagrangian volume from parent simulation.
2: Parent simulation halo location and surroundings.
3-5: Increasing resolutions for zoomed halos and surroundings (if available).

You need to change the mass cut offs and boxwidths below.
Also check the file paths as things may have moved

"""

datadir = '/bigbang/data/AnnaGroup/caterpillar/halos/'
icbasepath = '/bigbang/data/AnnaGroup/caterpillar/ics/lagr/'
parentpath = "/bigbang/data/AnnaGroup/caterpillar/parent/gL100X10/"
parenthalopath = parentpath + 'rockstar/'
img_path = '/bigbang/data/bgriffen/caterpillar/halos/'

mcut = 5e10
hubble = 0.6711
boxwidth = 4. #2.*4.*hubble

htable = asciitable.read("/bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index.txt",Reader=asciitable.FixedWidth)
key = [str(pid)+'_'+str(lx) for pid,lx in zip(htable['parentid'],htable['LX'])]
hindex = dict(zip(key,htable['zoomid']))

circut1 = 1.
circut2 = 3.
circut3 = 4.

halo_dirs = glob.glob(datadir+"H*")
#halo_dirs = halo_dirs + glob.glob(datadir+"/largehalos_forlater/H*")
#halo_list = [1194083,1292049,1327707,231858,1725139,230667,1269360,649524,706754]

# get host halos
parenthalos = RSDataReader.RSDataReader(parenthalopath,127,digits=3,version=6)
hosts = parenthalos.get_hosts()
subs = parenthalos.get_subs()

# remove hubble component
xposhosts_all = np.array(hosts['posX'])/hubble
yposhosts_all = np.array(hosts['posY'])/hubble
zposhosts_all = np.array(hosts['posZ'])/hubble
rvirhosts_all = np.array(hosts['rvir'])
mvirhosts_all = np.array(hosts['mvir'])/hubble

xpossubs_all = np.array(subs['posX'])/hubble
ypossubs_all = np.array(subs['posY'])/hubble
zpossubs_all = np.array(subs['posZ'])/hubble
rvirsubs_all = np.array(subs['rvir'])
mvirsubs_all = np.array(subs['mvir'])/hubble

mvirall = []

def get_mask(x,y,z,center,delta):
    cond1h = (np.array(x) >= center[0] - delta) & (np.array(x) <= center[0] + delta)
    cond2h = (np.array(y) >= center[1] - delta) & (np.array(y) <= center[1] + delta)
    cond3h = (np.array(z) >= center[2] - delta) & (np.array(z) <= center[2] + delta)
    return cond1h & cond2h & cond3h

def set_axes(ax1,ax2,ax3,width):
    ax1.set_xlim([-width,width])
    ax1.set_ylim([-width,width])
    ax2.set_xlim([-width,width])
    ax2.set_ylim([-width,width])
    ax3.set_xlim([-width,width])
    ax3.set_ylim([-width,width])

for haloi in halo_dirs:
    #haloid_check = haloi.split("NRVIR")[0].split("lagr/")[1][1:]
    haloid = haloi.split("/")[-1][1:]
    #if int(haloid) in halo_list:
    #print
    sub_dirs = glob.glob(haloi+"/H*")
    f, axs = plt.subplots(5, 3, figsize=(12,14))
    f.subplots_adjust(hspace=0,wspace=0)
    
    # ///////////////////////////////////////////////////////////
    #           plot lagrangian region of the parent
    # ///////////////////////////////////////////////////////////
    print "> CREATING LAGRANGIAN VOLUME PARTICLES PLOT! [%s]" % (haloid)
    pointfile = icbasepath + 'H'+str(haloid) + 'NRVIR4'
    #print "reading...",pointfile
    x_ic,y_ic,z_ic = getlagrxyz(pointfile)
    centx_ic,centy_ic,centz_ic,extx_ic,exty_ic,extz_ic = getcentext(pointfile + '.head')
    axi = 0
    ax1 = axs[axi,0]
    ax2 = axs[axi,1]
    ax3 = axs[axi,2]
    ax1.plot(x_ic,y_ic,'b.')
    ax2.plot(x_ic,z_ic,'b.')
    ax3.plot(y_ic,z_ic,'b.')
    xbox_ic = centx_ic - 0.5*extx_ic
    ybox_ic = centy_ic - 0.5*exty_ic
    zbox_ic = centz_ic - 0.5*extz_ic
    ax1.add_patch(patches.Rectangle((xbox_ic,ybox_ic),extx_ic,exty_ic,facecolor='none',
                                                         edgecolor='r',
                                                         linestyle='dashed'))
    ax2.add_patch(patches.Rectangle((xbox_ic,zbox_ic),extx_ic,extz_ic,facecolor='none',
                                                         edgecolor='r',
                                                         linestyle='dashed'))
    ax3.add_patch(patches.Rectangle((ybox_ic,zbox_ic),exty_ic,extz_ic,facecolor='none',
                                                         edgecolor='r',
                                                         linestyle='dashed'))

    ax1.set_xlim(x_ic.min()*0.9,x_ic.max()*1.1)
    ax1.set_ylim(y_ic.min()*0.9,y_ic.max()*1.1)
    ax2.set_xlim(x_ic.min()*0.9,x_ic.max()*1.1)
    ax2.set_ylim(z_ic.min()*0.9,z_ic.max()*1.1)
    ax3.set_xlim(y_ic.min()*0.9,y_ic.max()*1.1)
    ax3.set_ylim(z_ic.min()*0.9,z_ic.max()*1.1)
    ax1.set_ylabel('(y,y,z) [boxwidth]')
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax3.set_xticklabels([])
    ax3.set_yticklabels([])
    placetext(ax1,0.05,0.9,'LAGR','bold',12)
    # ///////////////////////////////////////////////////////////
    #                   plot parent host halos
    # ///////////////////////////////////////////////////////////
    print "> CREATING HOST HALO PLOT!"
    axi += 1
    ax1 = axs[axi,0]
    ax2 = axs[axi,1]
    ax3 = axs[axi,2]
    
    # extract host information
    #mask = (parenthalos['id']==haloid)
    haloid_int = int(haloid)
    xcen = parenthalos.ix[haloid_int]['posX']/hubble
    ycen = parenthalos.ix[haloid_int]['posY']/hubble
    zcen = parenthalos.ix[haloid_int]['posZ']/hubble
    rvirch = parenthalos.ix[haloid_int]['rvir']/1000.
    mvirch = parenthalos.ix[haloid_int]['mvir']/hubble
  
    # construct position and mass mask for region centered on host
    mask_pos_hosts = get_mask(xposhosts_all,yposhosts_all,zposhosts_all,[xcen,ycen,zcen],boxwidth/2.)
    mask_mass_hosts = (mvirhosts_all >= mcut/10.)
    condh = mask_pos_hosts & mask_mass_hosts
    
    xposhosts = np.array(xposhosts_all[condh])-xcen
    yposhosts = np.array(yposhosts_all[condh])-ycen
    zposhosts = np.array(zposhosts_all[condh])-zcen
    rvirhosts = np.array(rvirhosts_all[condh])/1000.
    mvirhosts = mvirhosts_all[condh]
    pos_hosts = np.array([xposhosts,yposhosts,zposhosts]).T
    for xposi,yposi,zposi,mviri,rviri in zip(xposhosts,yposhosts,zposhosts,mvirhosts,rvirhosts):
        if mviri < mvirch and mviri > mcut:
            mstring = '{0:.2e}'.format(mviri)
            placetext_direct(ax1,xposi,yposi,mstring,'normal',10)
            placetext_direct(ax2,xposi,zposi,mstring,'normal',10)
            placetext_direct(ax3,yposi,zposi,mstring,'normal',10)
    
    # construct position and mass mask for subhalos
    mask_pos_sub = get_mask(xpossubs_all,ypossubs_all,zpossubs_all,[xcen,ycen,zcen],boxwidth/2.)
    mask_mass_sub = (mvirsubs_all >= mcut/10.)

    conds = mask_pos_sub & mask_mass_sub

    xpossubs = np.array(xpossubs_all[conds])-xcen
    ypossubs = np.array(ypossubs_all[conds])-ycen
    zpossubs = np.array(zpossubs_all[conds])-zcen
    rvirsubs = np.array(rvirsubs_all[conds])/1000.
    pos_subs = np.array([xpossubs,ypossubs,zpossubs]).T
    # draw subhalos which mate the cut
    plotxyzprojr(ax1,ax2,ax3,pos_subs,rvirsubs,axislabels=False,color='r',linestyle='-',linewidth=1)
    # draw hosts which made the cut
    plotxyzprojr(ax1,ax2,ax3,pos_hosts,rvirhosts,axislabels=False,color='b',linestyle='-',linewidth=1)
    # draw central host
    plotxyzprojr(ax1,ax2,ax3,np.array([0,0,0]).T,rvirch,axislabels=False,color='b',linestyle='-',linewidth=2)
    # exclusion zones
    plotxyzprojr(ax1,ax2,ax3,np.array([0,0,0]).T,circut1,axislabels=False,color='c',linestyle='--',linewidth=2)
    plotxyzprojr(ax1,ax2,ax3,np.array([0,0,0]).T,circut2,axislabels=False,color='c',linestyle='-.',linewidth=2)
    plotxyzprojr(ax1,ax2,ax3,np.array([0,0,0]).T,circut3,axislabels=False,color='c',linestyle='-',linewidth=2)
    # add label for host
    placetext_direct(ax1,0.,0.,'{0:.2e}'.format(mvirch),'bold',12)
    placetext_direct(ax2,0.,0.,'{0:.2e}'.format(mvirch),'bold',12)
    placetext_direct(ax3,0.,0.,'{0:.2e}'.format(mvirch),'bold',12)
    placetext(ax1,0.05,0.9,'PARENT','bold',12)
    print '  -- # hosts:',len(xposhosts)
    print '  -- # subs:  ',len(xpossubs)

    # plot formatting
    ax1.set_ylabel('(y,y,z) [Mpc]')
    set_axes(ax1,ax2,ax3,boxwidth/2.)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    ax3.set_xticklabels([])
    ax3.set_yticklabels([])
    # ///////////////////////////////////////////////////////////
    
    print "> PLOTTING ZOOM IN HALOS!"

    for sub_diri in sub_dirs:
        LX = sub_diri.split("LX")[1][:2]
        if haloid+'_'+LX in hindex.keys():
            print "> DOING HALO:",haloid,", LX:",LX
            axi += 1
            ax1 = axs[axi,0]
            ax2 = axs[axi,1]
            ax3 = axs[axi,2]
            
            # load catalogues
            halodata = RSDataReader.RSDataReader(sub_diri+"/halos",255,version=7)
            allhalos = halodata.data

            hosts = halodata.get_hosts()
            condhosts = (np.array(hosts['mvir'])/hubble >= mcut/10.)

            # central host halo
            zoomid = hindex[haloid+'_'+LX]
            candhalos = allhalos.ix[zoomid]
            xhost = float(allhalos.ix[zoomid]['posX'])/hubble
            yhost = float(allhalos.ix[zoomid]['posY'])/hubble
            zhost = float(allhalos.ix[zoomid]['posZ'])/hubble
            mvirhost = float(allhalos.ix[zoomid]['mvir'])/hubble
            rvirhost = float(allhalos.ix[zoomid]['rvir'])
            plotxyzprojr(ax1,ax2,ax3,np.array([0,0,0]).T,rvirhost,axislabels=False,color='b',linestyle='-',linewidth=2)
            # draw hosts
            xposhosts = np.array(hosts['posX'][condhosts])/hubble-xhost
            yposhosts = np.array(hosts['posY'][condhosts])/hubble-yhost
            zposhosts = np.array(hosts['posZ'][condhosts])/hubble-zhost
            rvirhosts = hosts['rvir'][condhosts]/1000.
            mvirhosts = hosts['mvir'][condhosts]/hubble
            hostids = hosts['id'][condhosts]

            pos_zoom_hosts = np.array([xposhosts,yposhosts,zposhosts]).T
            plotxyzprojr(ax1,ax2,ax3,pos_zoom_hosts,rvirhosts,axislabels=False,color='b',linestyle='-',linewidth=1)
            
            mstring = '{0:.2e}'.format(mvirhost)
            placetext_direct(ax1,0.,0.,mstring,'bold',12)
            placetext_direct(ax2,0.,0.,mstring,'bold',12)
            placetext_direct(ax3,0.,0.,mstring,'bold',12)
            placetext(ax1,0.05,0.9,"LX: "+LX,'bold',12)

            for idi,x,y,z,rvir,mvir in zip(hostids,xposhosts,yposhosts,zposhosts,rvirhosts,mvirhosts):
                if mvir > mcut and idi != haloid and mvir < mvirhost and \
                    x > -boxwidth/2. and x < boxwidth/2. and \
                    y > -boxwidth/2. and y < boxwidth/2. and \
                    z > -boxwidth/2. and z < boxwidth/2.:
                        mstring = '{0:.2e}'.format(mvir)
                        placetext_direct(ax1,x,y,mstring,'normal',10)
                        placetext_direct(ax2,x,z,mstring,'normal',10)
                        placetext_direct(ax3,y,z,mstring,'normal',10)

            # draw subhalos
            subs = halodata.get_subs()
            condsubs = (np.array(subs['mvir'])/hubble >= mcut/.10)
            
            xpossubs = np.array(subs['posX'][condsubs]/hubble-xhost)
            ypossubs = np.array(subs['posY'][condsubs]/hubble-yhost)
            zpossubs = np.array(subs['posZ'][condsubs]/hubble-zhost)
            rvirsubs = subs['rvir'][condsubs]/1000.
            mvirsubs = subs['mvir'][condsubs]/0.6711

            pos_zoom_subs = np.array([xpossubs,ypossubs,zpossubs]).T
            plotxyzprojr(ax1,ax2,ax3,pos_zoom_subs,rvirsubs,axislabels=False,color='r',linestyle='-',linewidth=2)
            
            plotxyzprojr(ax1,ax2,ax3,np.array([0,0,0]).T,circut1,axislabels=False,color='c',linestyle='--',linewidth=2)
            plotxyzprojr(ax1,ax2,ax3,np.array([0,0,0]).T,circut2,axislabels=False,color='c',linestyle='-.',linewidth=2)
            plotxyzprojr(ax1,ax2,ax3,np.array([0,0,0]).T,circut3,axislabels=False,color='c',linestyle='-',linewidth=2)
            #ax1.set_xlabel('x-pos [Mpc]')
            ax1.set_ylabel('(y,z,z) [Mpc]')
            if LX != '13':
                ax2.set_xticklabels([])
                ax3.set_xticklabels([])
                
            ax2.set_yticklabels([])
            ax3.set_yticklabels([])
            if LX == '13':
                ax1.set_xlabel('x [Mpc]')
                ax2.set_xlabel('x [Mpc]')
                ax3.set_xlabel('y [Mpc]')
            set_axes(ax1,ax2,ax3,boxwidth/2.)
            
    f.savefig(img_path + 'H'+str(haloid)+'_XYZ_HALOS.png',bbox_inches='tight')
    print "Writing to file:",img_path + 'H'+str(haloid)+'_xyz.png'
    
    #f.savefig('/bigbang/data/bgriffen/caterpillar/spatial/H'+str(haloid)+'-xyz-zoom_nofix.png',bbox_inches='tight')
    
    plt.close(f)
        #sys.exit()
        