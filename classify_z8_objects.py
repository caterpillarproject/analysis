import numpy as np
import time,sys
import pandas as pd
import cPickle as pickle
import matplotlib.pyplot as plt

import haloutils
from caterpillaranalysis import MassAccrPlugin

import sys
sys.path.append("./greg_dwarfs")
import DwarfMethods as dm
import abundance_matching

from collections import Counter

from select_z8_objects import zin_to_zr_snapr

import seaborn as sns
sns.set(context='poster',style='ticks',font='serif',palette='colorblind', font_scale=1.5)
sns.set_style({"xtick.direction":"in","ytick.direction":"in"})

## Global logMbins
logMbins = np.arange(4,10,.1)
logMbinsmid = (logMbins[1:]+logMbins[:-1])/2.
allufdtypes = ["maxm","h14m","h14r","h14i"]
ufdtypescolors = ['b','g','m','c']
ufdlinestyles = ['-','--',':','-.']

## Global Vmaxbins
logVmaxbins = np.arange(0,2.6,.05)
logVmaxbinsmid = (logVmaxbins[1:]+logVmaxbins[:-1])/2.

## Global Concentration bins
concbins = np.arange(0,20,.5)
concbinsmid = (concbins[1:]+concbins[:-1])/2.

## Global distance bins (in kpc/h)
logDbins = np.arange(-1,4,.1)
logDbinsmid = (logDbins[1:]+logDbins[:-1])/2.

## Global spin bins
spinbins = np.arange(0,.5,.01)
spinbinsmid = (spinbins[1:]+spinbins[:-1])/2.

## T/|U| bins
TUbins = np.arange(0,2.5,.1)
TUbinsmid = (TUbins[1:]+TUbins[:-1])/2.
logetabins = np.arange (0,1.01,.05)
logetabinsmid = np.arange (0,1.01,.05)

## NOTE: replacing Vmax with T/|U|
#all_bins = [logMbins, logVmaxbins, spinbins, concbins, logDbins, TUbins]
#all_bins_mid = [logMbinsmid, logVmaxbinsmid, spinbinsmid, concbinsmid, logDbinsmid, TUbinsmid]
all_bins = [logMbins, TUbins, spinbins, concbins, logDbins]
all_bins_mid = [logMbinsmid, TUbinsmid, spinbinsmid, concbinsmid, logDbinsmid]
#all_bins = [logMbins, logetabins, spinbins, concbins, logDbins]
#all_bins_mid = [logMbinsmid, logetabinsmid, spinbinsmid, concbinsmid, logDbinsmid]

## Property labels
#prop_labels = ["logM (Msun)", "logVmax (km/s)", "spin", "conc", "logdist (kpc/h)"]
#prop_xlims = [(6.5,10), (0.5,2.5), (0,0.5), (0,20), (0,4)]
prop_labels = ["logM (Msun)", "T/|U|", "spin", "conc", "logdist (kpc/h)"]
#prop_labels = ["logM (Msun)", "log 2T/|U|", "spin", "conc", "logdist (kpc/h)"]
prop_xlims = [(6.5,10), (0.0,2.5), (0,0.5), (0,20), (0,4)]
#prop_xlims = [(6.5,10), (0.0,1.0), (0,0.5), (0,20), (0,4)]

def load_one_halo_data(zin,hpath,use_vmaxconc=False):
    """
    Load quantities at zin for halo (hpath) that could correlate.
    """
    h0 = 0.6711
    z_r, snap_r = zin_to_zr_snapr(zin, verbose=False)
    
    hid = haloutils.get_parent_hid(hpath)


    ## Merger tree quantities
    # Recall that I have cut logM > 10**6.5
    with open("UFDSEARCH_Z0/{}_ufdids.pkl".format(haloutils.hidstr(hid)),"r") as fp:
        allufdids = pickle.load(fp)
    zrobjs = np.load("UFDSEARCH_Z0/{}_z{}halos.npy".format(haloutils.hidstr(hid),zin))
    
    zrlogmass = np.log10(zrobjs['mvir']/h0)
    zrlogvmax = np.log10(zrobjs['vmax'])
    
    # Spin: 
    zrspin = zrobjs['spin']

    # Distance to host
    plug = MassAccrPlugin()
    mb = plug.read(hpath)
    ii = mb['snap'] == snap_r
    assert np.sum(ii) == 1, "{} has {}".format(hid,np.sum(ii))
    halopos = zrobjs[['posX','posY','posZ']].view(np.float).reshape(-1,3)
    hostpos = np.array(mb[ii][['x','y','z']]).view(np.float).reshape(-1,3)
    zrlogD = np.log10(np.sqrt(np.sum((halopos - hostpos)**2,axis=1))) + 3
    
    # Concentration: have to load rscat, will lose phantoms

    # Concentration with rvmax and vmax
    # Dooley et al. 2014 Eqn 9
    if use_vmaxconc:
        rscat = haloutils.load_rscat(hpath, snap_r, rmaxcut=False)
        rscatobjs = rscat.ix[zrobjs['origid']]
        conc_const = 0.21639 * (10/h0)**2 # scaling constant * (kpc/(km/s)/H0)**2
        conc = conc_const * np.array(rscatobjs['vmax']/rscatobjs['rvmax'])**2
    else:
        #conc = rscatobjs['rvir']/rscatobjs['rs']
        rscatobjs = None
        conc = zrobjs['rvir']/zrobjs['rs']

    return zrobjs, rscatobjs, zrlogmass, zrlogvmax, zrspin, zrlogD, conc

def histogram_one_halo(zin,hpath,full_output=False, use_vmaxconc=False):
    zrobjs, rscatobjs, zrlogmass, zrlogvmax, zrspin, zrlogD, conc = load_one_halo_data(zin, hpath,
                                                                                       use_vmaxconc=use_vmaxconc)
    zrTU = zrobjs["T/|U|"]

    hid = haloutils.get_parent_hid(hpath)
    with open("UFDSEARCH_Z0/{}_ufdids.pkl".format(haloutils.hidstr(hid)),"r") as fp:
        allufdids = pickle.load(fp)
    # Mass
    h1, _ = np.histogram(zrlogmass, bins=logMbins)
    ### Vmax
    ##h2, _ = np.histogram(zrlogvmax, bins=logVmaxbins)
    # T/|U|
    h2, _ = np.histogram(zrTU, bins=TUbins)
    # Spin
    h3, _ = np.histogram(zrspin, bins=spinbins)
    # Conc
    _conc = conc[np.isfinite(conc)]
    h4, _ = np.histogram(_conc, bins=concbins)
    # Dist
    h5, _ = np.histogram(zrlogD, bins=logDbins)

    h1ufds = []
    h2ufds = []
    h3ufds = []
    h4ufds = []
    h5ufds = []
    for ufdids, ufdtype in zip(allufdids,allufdtypes):
        ufdids = np.array(ufdids).astype(int)
        ii_ufd = np.array(map(lambda x: x in ufdids, zrobjs['mtkey']))
        
        hufd, _ = np.histogram(zrlogmass[ii_ufd], bins=logMbins)
        h1ufds.append(hufd)
        #hufd, _ = np.histogram(zrlogvmax[ii_ufd], bins=logVmaxbins)
        hufd, _ = np.histogram(zrTU[ii_ufd], bins=TUbins)
        h2ufds.append(hufd)
        hufd, _ = np.histogram(zrspin[ii_ufd], bins=spinbins)
        h3ufds.append(hufd)
        _conc = conc[ii_ufd]
        _conc = _conc[np.isfinite(_conc)]
        hufd, _ = np.histogram(_conc, bins=concbins)
        h4ufds.append(hufd)
        hufd, _ = np.histogram(zrlogD[ii_ufd], bins=logDbins)
        h5ufds.append(hufd)

    if full_output:
        return [h1,h2,h3,h4,h5],[h1ufds,h2ufds,h3ufds,h4ufds,h5ufds],[zrobjs, rscatobjs, zrlogmass, zrlogvmax, zrspin, zrlogD, conc]
    return [h1,h2,h3,h4,h5],[h1ufds,h2ufds,h3ufds,h4ufds,h5ufds]

def histogram_all_halos(zin):
    ## match this to histogram_one_halo
    num_hists = len(all_bins_mid)
    all_data = [[] for i in range(num_hists)]
    all_data2= [[] for i in range(num_hists)]

    hpaths = dm.get_hpaths(field=False, lx=14)
    for hpath in hpaths:
        start = time.time()
        try:
            h_out,hufd_out = histogram_one_halo(zin, hpath)
            for x1,x2,h1,h2s in zip(all_data, all_data2, h_out, hufd_out):
                x1.append(h1)  # 32 x Nbin
                x2.append(h2s) # 32 x 4 x Nbin
        except Exception as e:
            print e
        print "{}: {:.2f}".format(haloutils.hidstr(haloutils.get_parent_hid(hpath)), time.time()-start)
            
    to_save = []
    for m,h1,h2s in zip(all_bins_mid,all_data,all_data2):
        to_save.append([m,np.array(h1),np.array(h2s)])
    with open("UFDSEARCH_Z0/haloprops_at_z{}.pkl".format(zin),"w") as fp:
        pickle.dump(to_save, fp)

def plot_histograms(zin, plot_type="count"):
    assert plot_type in ["count","ratio"], plot_type

    with open("UFDSEARCH_Z0/haloprops_at_z{}.pkl".format(zin),"r") as fp:
        all_data = pickle.load(fp)
    fig, axes = plt.subplots(3,2,figsize=(8*2,8*3))
    labels = prop_labels
    xlims = prop_xlims
    for data,ax,label,xlim in zip(all_data,axes.flat[0:5],labels,xlims):
        x, h1, h2s = data
        y1,y2,y3 = medianscatter(h1, percentiles=[50-95/2.,50.,50+95/2.], axis=0)
        if plot_type == "count":
            ax.plot(x, y2, 'k', drawstyle='steps-mid', label="All halos")
            ax.fill_between(x,y1,y3,color='k',facecolor='k',alpha=.3,step='mid')
        for j,(ufdtype, color, ls) in enumerate(zip(allufdtypes, ufdtypescolors, ufdlinestyles)):
            h2 = h2s[:,j,:]
            if plot_type == "count":
                y1,y2,y3 = medianscatter(h2, percentiles=[50-95/2.,50.,50+95/2.], axis=0)
            elif plot_type == "ratio":
                y1,y2,y3 = medianscatter(h2/h1.astype(float), percentiles=[50-95/2.,50.,50+95/2.], axis=0)
            ax.plot(x, y2, color, linestyle=ls, drawstyle='steps-mid', label=ufdtype)
            ax.fill_between(x,y1,y3,color=color,facecolor=color,alpha=.3,step='mid')

        ax.set_xlabel(label)
        ax.set_xlim(xlim)
        if plot_type == "count":
            ax.set_yscale("log")
        else:
            ax.set_ylim(0,1)
    ## Make legend in 6th axis
    ax = axes[0,0]
    handles,labels = ax.get_legend_handles_labels()
    lax = axes[2,1]
    lax.legend(handles,labels,loc='center')

    if plot_type == "count":
        fig.savefig("haloprops_at_z{}.pdf".format(zin))
    elif plot_type == "ratio":
        fig.savefig("halopropsratio_at_z{}.pdf".format(zin))
    return fig

def medianscatter(x, percentiles=None, **kwargs):
    if percentiles is None:
        percentiles = [0, 50-95/2., 50-68/2., 50, 50+68/2., 50+95/2., 100]
    out = []
    for p in percentiles:
        out.append(np.nanpercentile(x, p, **kwargs))
    return out
def calculate_total_list(alldata):
    all_total = []
    for data in alldata:
        host, ufds, subs = data
        total = np.sum([host, ufds[0], subs[0]], axis=0)
        #for i in [1,2,3]:
        #    assert (total == np.sum([host, ufds[i], subs[i]], axis=0)).all()
        all_total.append(total)
    return all_total


if __name__=="__main__":
    #try:
    #    zin = int(sys.argv[1])
    #    out = zin_to_zr_snapr(zin)
    #except:
    #    zin = 6
    #    print "Default zin: {}".format(zin)
    ##histogram_all_halos(zin)
    #plot_histograms(zin, plot_type="count")
    #plot_histograms(zin, plot_type="ratio")

    for zin in [4,6,8,10,12]:
        #histogram_all_halos(zin)
        fig = plot_histograms(zin, plot_type="count")
        plt.close(fig)
        fig = plot_histograms(zin, plot_type="ratio")
        plt.close(fig)

def run_count(zin):
    h0 = 0.6711
    z_r, snap_r = zin_to_zr_snapr(zin)

    #z_r = 8.346
    #snap_r = 67
    with open("UFDSEARCH_Z0/summary_data_z{}.pkl".format(zin),"r") as fp:
        all_hids, all_hists_m, all_hists_v, all_hists_maxm, all_num_missing = \
            pickle.load(fp)
        
    ## OUTPUT DATA:
    # A list of tuples, containing histograms wrt logMbins
    # objs in host, 4xobjs in ufds, 4xobjs in other subs
    all_data = []

    hpaths = dm.get_hpaths(field=False, lx=14)
    for hpath in hpaths:
        hid = haloutils.get_parent_hid(hpath)
        with open("UFDSEARCH_Z0/{}_ufdids.pkl".format(haloutils.hidstr(hid)),"r") as fp:
            allufdids = pickle.load(fp)
        z8objs = np.load("UFDSEARCH_Z0/{}_z{}halos.npy".format(haloutils.hidstr(hid),zin))
        #z8objs = z8objs[z8objs['mvir']/h0 > 10**logMmin]
        z8mass = np.log10(z8objs['mvir']/h0)

        ## Histogram Host
        zoomid = haloutils.load_zoomid(hpath)
        h_host, _ = np.histogram(z8mass[z8objs['mtkey'] == zoomid], logMbins)
        
        ## Histogram UFDs
        hs_ufds = []
        for (ufdids, ufdtype) in zip(allufdids,allufdtypes):
            ufdids = np.array(ufdids).astype(int)
            ii_ufd = np.array(map(lambda x: x in ufdids, z8objs['mtkey']))
            h, _ = np.histogram(z8mass[ii_ufd], logMbins)
            hs_ufds.append(h)

        ## Histogram Non-UFD Subs
        hs_subs = []
        for (ufdids, ufdtype) in zip(allufdids,allufdtypes):
            ufdids = np.array(ufdids).astype(int)
            ii_sub = np.array(map(lambda x: (x != zoomid) and (x not in ufdids), z8objs['mtkey']))
            h, _ = np.histogram(z8mass[ii_sub], logMbins)
            hs_subs.append(h)

        this_data = h_host, hs_ufds, hs_subs
        all_data.append(this_data)
    return all_data

def make_host_figure(alldata,plot_frac=True):
    all_total = calculate_total_list(alldata)
    all_hs = []
    all_fracs = []
    for data,total in zip(alldata, all_total):
        h = data[0]
        frac = h.astype(float)/total
        all_hs.append(h)
        all_fracs.append(frac)
    fig, ax = plt.subplots(figsize=(8,8))
    if plot_frac:
        toplot = medianscatter(all_fracs, axis=0)
    else:
        toplot = medianscatter(all_hs, axis=0)
    ax.plot(logMbinsmid, toplot[3], 'k')
    ax.fill_between(logMbinsmid, toplot[1], toplot[5], color='k', facecolor='k', alpha=.3)
    ax.fill_between(logMbinsmid, toplot[2], toplot[4], color='k', facecolor='k', alpha=.3)
    ax.set_xlabel(r'$\log M/M_\odot$ at $z_r={}$'.format(zin))
    ax.set_ylabel('fraction merged into host')
    return fig
    
def make_ufds_figure(alldata):
    all_total = calculate_total_list(alldata)
    ufdtypes_all_hs = []
    ufdtypes_all_fracs = []

    for ufdtype in allufdtypes:
        ufdtypes_all_hs.append([])
        ufdtypes_all_fracs.append([])

    for data,total in zip(alldata, all_total):
        for ufdtype, h, all_hs, all_fracs in \
        zip(allufdtypes, data[1], ufdtypes_all_hs, ufdtypes_all_fracs):
            frac = h.astype(float)/total
            all_hs.append(h)
            all_fracs.append(frac)
    fig, ax = plt.subplots(figsize=(8,8))
    for j in range(len(allufdtypes)):
        toplot = medianscatter(ufdtypes_all_fracs[j], axis=0)
        color = ufdtypescolors[j]
        ls = ufdlinestyles[j]
        ax.plot(logMbinsmid, toplot[3], color=color, label=allufdtypes[j], ls=ls)
        ax.fill_between(logMbinsmid, toplot[1], toplot[5], color=color, facecolor=color, alpha=.3)
        ax.fill_between(logMbinsmid, toplot[2], toplot[4], color=color, facecolor=color, alpha=.3)
    ax.legend(loc='upper right')
    ax.set_xlabel(r'$\log M/M_\odot$ at $z_r={}$'.format(zin))
    ax.set_ylabel('fraction survived as z=0 UFD')
    return fig

def make_subs_figure(alldata):
    all_total = calculate_total_list(alldata)
    ufdtypes_all_hs = []
    ufdtypes_all_fracs = []

    for ufdtype in allufdtypes:
        ufdtypes_all_hs.append([])
        ufdtypes_all_fracs.append([])

    for data,total in zip(alldata, all_total):
        for ufdtype, h, all_hs, all_fracs in \
        zip(allufdtypes, data[2], ufdtypes_all_hs, ufdtypes_all_fracs):
            frac = h.astype(float)/total
            all_hs.append(h)
            all_fracs.append(frac)
    fig, ax = plt.subplots(figsize=(8,8))
    for j in range(len(allufdtypes)):
        toplot = medianscatter(ufdtypes_all_fracs[j], axis=0)
        color = ufdtypescolors[j]
        ls = ufdlinestyles[j]
        ax.plot(logMbinsmid, toplot[3], color=color, label=allufdtypes[j], ls=ls)
        ax.fill_between(logMbinsmid, toplot[1], toplot[5], color=color, facecolor=color, alpha=.3)
        ax.fill_between(logMbinsmid, toplot[2], toplot[4], color=color, facecolor=color, alpha=.3)
    ax.legend(loc='lower left')
    ax.set_xlabel(r'$\log M/M_\odot$ at $z_r={}$'.format(zin))
    ax.set_ylabel('fraction merged into z=0 subhalo')
    return fig

def tmp():
#    out = run_count(zin)
#    with open("UFDSEARCH_Z0/histograms_of_z{}.pkl".format(zin),"w") as fp:
#        pickle.dump(out,fp)
    with open("UFDSEARCH_Z0/histograms_of_z{}.pkl".format(zin),"r") as fp:
        alldata = pickle.load(fp)
    fig = make_host_figure(alldata)
    fig.savefig("frac_of_z{}_objs_in_host.pdf".format(zin),bbox_inches='tight')
    fig.savefig("frac_of_z{}_objs_in_host.png".format(zin),bbox_inches='tight')
    fig = make_host_figure(alldata, plot_frac=False)
    fig.axes[0].set_xlim(6.5,9.0)
    fig.savefig("num_of_z{}_objs_in_host.pdf".format(zin),bbox_inches='tight')
    fig.savefig("num_of_z{}_objs_in_host.png".format(zin),bbox_inches='tight')
    fig = make_ufds_figure(alldata)
    fig.axes[0].set_xlim(6.5,9.0)
    fig.savefig("frac_of_z{}_objs_in_ufds.pdf".format(zin),bbox_inches='tight')
    fig.savefig("frac_of_z{}_objs_in_ufds.png".format(zin),bbox_inches='tight')
    fig = make_subs_figure(alldata)
    fig.axes[0].set_xlim(6.5,9.0)
    fig.savefig("frac_of_z{}_objs_in_subs.pdf".format(zin),bbox_inches='tight')
    fig.savefig("frac_of_z{}_objs_in_subs.png".format(zin),bbox_inches='tight')

    plt.close('all')
#    plt.show()

#    run_count(7.0)
#    run_count(7.5)
#    plt.boxplot(np.array(out))
#    plt.gca().set_xticklabels(["maxm","h14m","h14r","h14i"])
#    plt.savefig("fig_for_mohammad.pdf")
#    run_count(8.0)
#    run_count(8.5)
#    run_count(9.0)
