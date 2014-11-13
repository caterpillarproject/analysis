import readsnapshots.readsnapHDF5 as rsHD
import readsnapshots.readsnap as rs
import readhalos.readsubf
import readsnapshots.readids
import numpy as np
import pylab as plt
import sys
from random import randint
from matplotlib import *
from brendanlib.grifflib import *
import mergertrees.MTCatalogue as MT
from colorsys import hsv_to_rgb
from random import randint, uniform
from brendanlib.grifflib import getcaterpillarcandidates

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


# THIS NEEDS TO BE MODIFIED TO RECEIVE THE NEW CANDIDATE LIST
base_halo_path = "/bigbang/data/AnnaGroup/caterpillar/halos/"
haloidlist_full = glob.glob(base_halo_path + "H*")
candidatelist = []
for haloid in haloidlist_full:
    candidatelist.append(int(haloid.split("/H")[1]))

# PLOTTING SETTINGS
hubble = 0.6711
ticksize = 11
ncols = 0
icand = 0
legendfontsize = 10
axislabelfontsize = 14

# SET PATHS
img_path = './images/'
basepath = "/bigbang/data/AnnaGroup/caterpillar/parent/gL100X10/"
ext = "outputs/snapdir_127/snap_127"
extics = 'ics/ics'
halopath = basepath + 'rockstar'
header = rsHD.snapshot_header(basepath + ext + ".0.hdf5")


fig2 = plt.figure(figsize=(15.0,6.0))
ax1a = fig2.add_subplot(1,2,1)
ax2a = fig2.add_subplot(1,2,2)

for haloid in candidatelist:

    fig1 = plt.figure(figsize=(22.0,12.0))
    ax1 = fig1.add_subplot(3,5,1)
    ax2 = fig1.add_subplot(3,5,2)
    ax3 = fig1.add_subplot(3,5,3)
    ax4 = fig1.add_subplot(3,5,4)
    ax5 = fig1.add_subplot(3,5,5)
    ax6 = fig1.add_subplot(3,5,6)
    ax7 = fig1.add_subplot(3,5,7)
    ax8 = fig1.add_subplot(3,5,8)
    ax9 = fig1.add_subplot(3,5,9)
    ax10 = fig1.add_subplot(3,5,10)
    ax11 = fig1.add_subplot(3,5,11)
    ax12 = fig1.add_subplot(3,5,12)
    ax13 = fig1.add_subplot(3,5,13)
    ax14 = fig1.add_subplot(3,5,14)
    ax15 = fig1.add_subplot(3,5,15)

    plt.subplots_adjust(hspace=0.1,wspace=0.3)

    haloid = int(haloid)
    
    cat = MT.MTCatalogue(halopath + '/trees',indexbyrsid=False,haloids=[haloid])
    tree = cat[0]

    mainbranch = tree.getMainBranch()
    
    scale = mainbranch['scale']
    rvir = mainbranch['rvir']
    posX = mainbranch['posX']/hubble
    posY = mainbranch['posY']/hubble
    posZ = mainbranch['posZ']/hubble
    mvir = mainbranch['mvir']/header.hubble
    mall = mainbranch['m200c_all']/header.hubble
    m200b = mainbranch['m200b']/header.hubble
    vmax = mainbranch['vmax']
    vrms = mainbranch['vrms']
    rs = mainbranch['rs']
    xoff = mainbranch['xoff']
    voff = mainbranch['voff']
    virrat = mainbranch['T/|U|']
    bta = mainbranch['b_to_a']
    cta = mainbranch['c_to_a']
    spin = mainbranch['spin']
    spinb = mainbranch['spin_bullock']
    pecvx = mainbranch['pecVX']
    pecvy = mainbranch['pecVY']
    pecvz = mainbranch['pecVZ']
    jx = mainbranch['Jx']
    jy = mainbranch['Jy']
    jz = mainbranch['Jz']
    Ax = mainbranch['A[x]']
    Ay = mainbranch['A[y]']
    Az = mainbranch['A[z]']
    ctb = bta*(1/cta)
    Tch = calcT(ctb)

    npart = mvir/(header.massarr[1]*10**10/header.hubble)

    normmvir = mvir/mvir[0]
    normvmax = vmax**2/vmax[0]**2

    icand += 1

    print "---------------------------------------"
    print "         Candidate:",icand
    print "---------------------------------------"
    print "     rockstar id: %i" % (haloid)
    print "  merger tree id: %i" % (mainbranch['id'][0])
    print "           x-pos:",'{:.2f}'.format(posX[0]), "   \ [Mpc]"
    print "           y-pos:",'{:.2f}'.format(posY[0]), "   \ [Mpc]"
    print "           z-pos:",'{:.2f}'.format(posZ[0]), "   \ [Mpc]"
    print "            vmax:",'{:.2f}'.format(vmax[0]),"  \ [km/s]"
    print "     virial mass:",'{0:.2e}'.format(mvir[0]),"\ [Msol]"
    print "   virial radius:",'{:.2f}'.format(rvir[0]),"  \ [kpc]"
    print "---------------------------------------"

    plotquantl(ax1,scale,np.log10(mvir),'virial')
    plotquantl(ax1,scale,np.log10(mall),'+unbound')
    plotquantl(ax1,scale,np.log10(m200b),'m200')
    plotquant(ax2,scale,normmvir)
    plotquantl(ax3,scale,vmax,'max')
    plotquantl(ax3,scale,vrms,'rms')
    plotquant(ax4,scale,normvmax)
    plotquantl(ax5,scale,rvir,'virial')  
    plotquantl(ax5,scale,rs,'scale')
    plotquantl(ax6,scale,spin,'normal')
    plotquantl(ax6,scale,spinb,'bullock')
    plotquantl(ax7,scale,xoff,'position')
    plotquantl(ax7,scale,voff,'velocity')
    plotquant(ax8,scale,virrat)
    plotquantl(ax9,scale,np.log10(npart),haloid)
    plotquantl(ax10,scale,bta,'b/a')
    plotquantl(ax10,scale,cta,'c/a')
    plotquant(ax11,scale,Tch)
    plotquantl(ax12,scale,pecvx,'Vx')
    plotquantl(ax12,scale,pecvy,'Vy')
    plotquantl(ax12,scale,pecvz,'Vz')
    plotquantl(ax13,scale,jx,'Jx')
    plotquantl(ax13,scale,jy,'Jy')
    plotquantl(ax13,scale,jz,'Jz')
    plotquantl(ax14,scale,Ax,'Ax')
    plotquantl(ax14,scale,Ay,'Ay')
    plotquantl(ax14,scale,Az,'Az')

    plotquant(ax2a,scale,normvmax)
    plotquantl(ax1a,scale,normmvir,haloid)
    
    if icand % 3 != 0:
        ncols += 1
    
  
    new_tick_locations = np.array([.2, .5, .9])
    
    ax11.text(0.5, 0.05,'oblate', horizontalalignment='center', verticalalignment='center', color='black', fontsize=11, transform = ax11.transAxes)
    ax11.text(0.5, 0.95,'prolate', horizontalalignment='center', verticalalignment='center', color='black', fontsize=11, transform = ax11.transAxes)
    
    ax15.text(0.5, 0.5,'spare', horizontalalignment='center', verticalalignment='center', color='black', fontsize=11, transform = ax15.transAxes)
    
    ax1.set_xlim([0,1])
    ax2.set_xlim([0,1])
    ax3.set_xlim([0,1])
    ax4.set_xlim([0,1])
    ax5.set_xlim([0,1])
    ax6.set_xlim([0,1])
    ax7.set_xlim([0,1])
    ax8.set_xlim([0,1])
    ax9.set_xlim([0,1])
    ax10.set_xlim([0,1])
    ax11.set_xlim([0,1])
    ax12.set_xlim([0,1])
    ax13.set_xlim([0,1])
    ax14.set_xlim([0,1])
    ax15.set_xlim([0,1])
    ax11.set_ylim([0,1])
    
    plotlegend(ax1,legendfontsize,location=4)
    plotlegend(ax3,legendfontsize)
    plotlegend(ax5,legendfontsize)
    plotlegend(ax6,legendfontsize)
    plotlegend(ax7,legendfontsize)
    plotlegend(ax10,legendfontsize,location=4)
    plotlegend(ax12,legendfontsize)
    plotlegend(ax13,legendfontsize)
    plotlegend(ax14,legendfontsize)
    
    new_tick_locations = np.array([0.,.2, .4, 0.6, 0.8, 1.])
    
    ax1top = ax1.twiny()
    ax1top.set_xticklabels(tick_function(new_tick_locations))
    
    ax2top = ax2.twiny()
    ax2top.set_xticklabels(tick_function(new_tick_locations))
    
    ax3top = ax3.twiny()
    ax3top.set_xticklabels(tick_function(new_tick_locations))
    
    ax4top = ax4.twiny()
    ax4top.set_xticklabels(tick_function(new_tick_locations))
    
    ax5top = ax5.twiny()
    ax5top.set_xticklabels(tick_function(new_tick_locations))
    
    ax1top.set_xlabel(r'$\mathrm{redshift}$',size=14)
    ax2top.set_xlabel(r'$\mathrm{redshift}$',size=14)
    ax3top.set_xlabel(r'$\mathrm{redshift}$',size=14)
    ax4top.set_xlabel(r'$\mathrm{redshift}$',size=14)
    ax5top.set_xlabel(r'$\mathrm{redshift}$',size=14)
    
    ax1.set_xticks([])
    ax2.set_xticks([])
    ax3.set_xticks([])
    ax4.set_xticks([])
    ax5.set_xticks([])
    ax6.set_xticks([])
    ax7.set_xticks([])
    ax8.set_xticks([])
    ax9.set_xticks([])
    ax10.set_xticks([])
    
    ax1.set_ylabel(r'$\mathrm{log_{10}\ M(z)\ [M_\odot]}$',size=axislabelfontsize)
    ax2.set_ylabel(r'$\mathrm{M_v(z)/M_v(z=0)}$',size=axislabelfontsize)
    ax3.set_ylabel(r'$\mathrm{V(z)\ [km/s]}$',size=axislabelfontsize)
    ax4.set_ylabel(r'$\mathrm{V_{max}(z)^2/V_{max}(z=0)^2}$',size=axislabelfontsize)
    ax5.set_ylabel(r'$\mathrm{radius\ [kpc]}$',size=axislabelfontsize)
    ax6.set_ylabel(r'$\mathrm{spin}$',size=axislabelfontsize)
    ax7.set_ylabel(r'$\mathrm{offset}$',size=axislabelfontsize)
    ax8.set_ylabel(r'$\mathrm{virial\ ratio}$',size=axislabelfontsize)
    ax9.set_ylabel(r'$\mathrm{log_{10}\ number\ of\ particles}$',size=axislabelfontsize)
    ax10.set_ylabel(r'$\mathrm{axis\ ratios}$',size=axislabelfontsize)
    ax11.set_ylabel(r'$\mathrm{triaxiality\ parameter}$',size=axislabelfontsize)
    ax12.set_ylabel(r'$\mathrm{peculiar\ V_x,\ V_y,\ V_z\ [km/s]}$',size=axislabelfontsize)
    ax13.set_ylabel(r'$\mathrm{J_x, J_y, J_z}$',size=axislabelfontsize)
    ax14.set_ylabel(r'$\mathrm{ellipticity\ axis}$',size=axislabelfontsize)
    
    ax11.set_xlabel(r'$\mathrm{scale\ factor}$',size=axislabelfontsize)
    ax12.set_xlabel(r'$\mathrm{scale\ factor}$',size=axislabelfontsize)
    ax13.set_xlabel(r'$\mathrm{scale\ factor}$',size=axislabelfontsize)
    ax14.set_xlabel(r'$\mathrm{scale\ factor}$',size=axislabelfontsize)
    ax15.set_xlabel(r'$\mathrm{scale\ factor}$',size=axislabelfontsize)
    
    ax1a.set_xlim([0,1])
    ax2a.set_xlim([0,1])
    ax1a.set_ylim([0,1.1])
    ax2a.set_ylim([0,1.1])

    ax1a.set_ylabel(r'$\mathrm{M_v(z)/M_v(z=0)}$',size=axislabelfontsize)
    #ax1a.set_ylabel(r'$\mathrm{log_{10}\ M\ [M_\odot]}$',size=axislabelfontsize)
    ax2a.set_ylabel(r'$\mathrm{V_{max}(z)^2/V_{max}(z=0)^2}$',size=axislabelfontsize)
    #ax2a.set_ylabel(r'$\mathrm{V_{circ}\ [km/s]}$',size=axislabelfontsize)
    
    ax1atop = ax1a.twiny()
    ax1atop.set_xticklabels(tick_function(new_tick_locations))
    
    ax2atop = ax2a.twiny()
    ax2atop.set_xticklabels(tick_function(new_tick_locations))
    
    ax1atop.set_xlabel(r'$\mathrm{redshift}$',size=14)
    ax2atop.set_xlabel(r'$\mathrm{redshift}$',size=14)
    ax1a.set_xlabel(r'$\mathrm{scale\ factor}$',size=axislabelfontsize)
    ax2a.set_xlabel(r'$\mathrm{scale\ factor}$',size=axislabelfontsize)

    fig1.savefig(img_path + 'H'+str(haloid)+'_all_mt_properties.png',bbox_inches='tight')
    plt.close(fig1)

fig2.savefig(img_path + 'all_stacked_mvir_vmax.png',bbox_inches='tight')

