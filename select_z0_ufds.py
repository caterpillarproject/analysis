import numpy as np
import haloutils
import time,sys
import pandas as pd
import matplotlib.pyplot as plt

import MTanalysis3 as mta
import MTaddition as mtadd
from FindMiniHalos import mcrit

import sys
sys.path.append("./greg_dwarfs")
import DwarfMethods as dm
import abundance_matching

if __name__=="__main__":
    hpaths = dm.get_hpaths(field=False, lx=14)
    hpath = hpaths[0]  # or whatever hpath you want
    #hpath = hpaths[4]
    #hpath = hpaths[7]

#def compute_histograms(hpath):
    hid = haloutils.get_parent_hid(hpath)
    h0 = .6711

    ## Abundance Matching for Mpeak limits
    ## UFD defined as within 1000 to 2e5 Mstar
    ## More massive is classical dSph, less massive is HFD
    am = abundance_matching.GarrisonKimmel()
    LMmin, LMmax = am.stellar_to_halo_mass([1000., 2.e5])

    #### From MTaddition.py:
    #### zsnaps = [90,78,67,59,53,47,42,38,34]  # corresponds to z = 6.33, 7.26, 8.346, 9.33, 10.22, 11.28, 12.33, 13.31, 14.44
    ## Get snap_r
    z_r = 8.346
    snap_r = 67
    #print np.log10(mcrit(1e4, z_r))
    
    print "starting rs z={} load".format(z_r); start = time.time()
    reion_rscat = haloutils.load_rscat(hpath, snap_r, rmaxcut=False)
    z0_rscat = haloutils.load_rscat(hpath, 319, rmaxcut=False)
    logm200_zr_all = np.log10(reion_rscat["altm2"]) #m200c
    logm200_zr_all = logm200_zr_all[np.isfinite(logm200_zr_all)]
    vmax_zr_all = reion_rscat["vmax"]
    print time.time()-start
    
    ## load data
    data=dm.get_extant_data(hpath,field=False)

    logm200_zr = np.log10(data["m200_8"]/h0)
    logm200_zr = logm200_zr[np.isfinite(logm200_zr)]
    vmax_zr = data["vmax_8"]
    vmax_z0 = np.array(z0_rscat.ix[data["rsid"]]["vmax"])
    
    ## Mmax criteria from abundance matching
    Mmaxm  = data["max_mass"]/h0
    ii_maxm = np.logical_and(Mmaxm > LMmin, Mmaxm < LMmax)
    
    ## Hargis et al. 2014 criteria
    ## Note Vpeak in ELVIS is Vmax at Mmax!
    Vpeak = data["max_mass_vmax"]
    zpeak = haloutils.get_z_snap(hpath, data["max_mass_snap"])
    ii_h14_cdSph = Vpeak <= 25.0
    # Massive in the past: Vpeak > 12 km/s, ~200 subhalos
    ii_h14m = Vpeak > 15.0
    # Formed before reionization: z > 8, Npart > 32, ~170 subhalos
    #ii_h14r = np.logical_and(Vpeak < 20.0, data["m200_8"]/h0 > mcrit(1e4, z_r))
    ii_h14r = np.logical_and(Vpeak < 20.0, data["m200_8"]/h0 > 1e7)
    #ii_h14r = np.logical_and(Vpeak < 20.0, data["m200_8"]/h0 > 6080000.)
    # Earliest infall: zpeak > 3, Vmax@z=0 > 8, ~56 subhalos
    #ii_h14i = np.logical_and(zpeak > 3.0, vmax_zr > 8.0)
    ii_h14i = np.logical_and(zpeak > 3.0, vmax_z0 > 8.0)
    
    ii_maxm = np.logical_and(ii_h14_cdSph, ii_maxm)
    ii_h14m = np.logical_and(ii_h14_cdSph, ii_h14m)
    ii_h14r = np.logical_and(ii_h14_cdSph, ii_h14r)
    ii_h14i = np.logical_and(ii_h14_cdSph, ii_h14i)

    #print "starting mtc load"; start = time.time()
    #mtc = haloutils.load_zoom_mtc(hpath, indexbyrsid=True)
    #print time.time()-start
    
    logMbins = np.arange(4,10,.1)
    vmaxbins = np.arange(0,50,1)

    fig, axes = plt.subplots(2,2,figsize=(10,10))
    ax = axes[0,0]
    all_h_m = []
    h,x,p = ax.hist(logm200_zr_all,      label='all',  histtype='step', bins=logMbins, color='k'); all_h_m.append(h)
    h,x,p = ax.hist(logm200_zr,          label='surv', histtype='step', bins=logMbins); all_h_m.append(h)
    h,x,p = ax.hist(logm200_zr[ii_maxm], label='maxm', histtype='step', bins=logMbins); all_h_m.append(h)
    h,x,p = ax.hist(logm200_zr[ii_h14m], label='h14m', histtype='step', bins=logMbins); all_h_m.append(h)
    h,x,p = ax.hist(logm200_zr[ii_h14r], label='h14r', histtype='step', bins=logMbins); all_h_m.append(h)
    h,x,p = ax.hist(logm200_zr[ii_h14i], label='h14i', histtype='step', bins=logMbins); all_h_m.append(h)
    ax.set_xlabel("M200c at z=8")
    ax.set_yscale('log')
    ax.text(.95,.93,haloutils.hidstr(hid), transform=ax.transAxes, ha='right')
    ax.text(.95,.87,"{:.2e}".format(float(haloutils.load_haloprops(hpath)[0])), transform=ax.transAxes, ha='right')
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.add_patch(plt.Rectangle((xlim[0],ylim[0]),6-xlim[0],ylim[1]-ylim[0], alpha=.2, color='k'))

    ax = axes[0,1]
    all_h_v = []
    h,x,p = ax.hist(vmax_zr_all,      label='all',  histtype='step', bins=vmaxbins, color='k'); all_h_v.append(h)
    h,x,p = ax.hist(vmax_zr,          label='surv', histtype='step', bins=vmaxbins); all_h_v.append(h)
    h,x,p = ax.hist(vmax_zr[ii_maxm], label='maxm', histtype='step', bins=vmaxbins); all_h_v.append(h)
    h,x,p = ax.hist(vmax_zr[ii_h14m], label='h14m', histtype='step', bins=vmaxbins); all_h_v.append(h)
    h,x,p = ax.hist(vmax_zr[ii_h14r], label='h14r', histtype='step', bins=vmaxbins); all_h_v.append(h)
    h,x,p = ax.hist(vmax_zr[ii_h14i], label='h14i', histtype='step', bins=vmaxbins); all_h_v.append(h)
    ax.set_xlabel("vmax at z=8")
    ax.legend(loc='upper right', fontsize=10)
    ax.set_yscale('log')
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.add_patch(plt.Rectangle((0,ylim[0]),4,ylim[1]-ylim[0], alpha=.2, color='k'))
    
    ax = axes[1,0]
    h_all, h_surv, h_maxm, h_h14m, h_h14r, h_h14i = all_h_m
    h_surv = h_surv.astype(float)
    logMbinsmid = (logMbins[1:]+logMbins[:-1])/2.
    ax.plot(logMbinsmid, h_surv/h_surv, drawstyle='steps-mid')
    ax.plot(logMbinsmid, h_maxm/h_surv, drawstyle='steps-mid')
    ax.plot(logMbinsmid, h_h14m/h_surv, drawstyle='steps-mid')
    ax.plot(logMbinsmid, h_h14r/h_surv, drawstyle='steps-mid')
    ax.plot(logMbinsmid, h_h14i/h_surv, drawstyle='steps-mid')
    ax.set_xlabel("M200c at z=8")
    ax.set_ylabel("frac of surviving")

    ax.text(.05,.95,"Nsurv={}".format(len(Mmaxm)), transform=ax.transAxes)
    ax.text(.05,.90,"Nmaxm={}".format(np.sum(ii_maxm)), transform=ax.transAxes)
    ax.text(.05,.85,"Nh14m={}".format(np.sum(ii_h14m)), transform=ax.transAxes)
    ax.text(.05,.80,"Nh14r={}".format(np.sum(ii_h14r)), transform=ax.transAxes)
    ax.text(.05,.75,"Nh14i={}".format(np.sum(ii_h14i)), transform=ax.transAxes)

    ax = axes[1,1]
    ax.hist(np.log10(Mmaxm),          label='surv', histtype='step', bins=logMbins)
    ax.hist(np.log10(Mmaxm[ii_maxm]), label='maxm', histtype='step', bins=logMbins)
    ax.hist(np.log10(Mmaxm[ii_h14m]), label='h14m', histtype='step', bins=logMbins)
    ax.hist(np.log10(Mmaxm[ii_h14r]), label='h14r', histtype='step', bins=logMbins)
    ax.hist(np.log10(Mmaxm[ii_h14i]), label='h14i', histtype='step', bins=logMbins)
    ax.text(.05,.94,"Sharp cutoff at 7.5 imposed by Greg", transform=ax.transAxes)
    ax.set_xlabel("Mvir at Mmax")
    ax.set_xlim(7,10)
    plt.show()

