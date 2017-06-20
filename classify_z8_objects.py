import numpy as np
import haloutils
import time,sys
import pandas as pd
import cPickle as pickle
import matplotlib.pyplot as plt

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

def medianscatter(x, percentiles=None, **kwargs):
    if percentiles is None:
        percentiles = [0, 50-95/2., 50-68/2., 50, 50+68/2., 50+95/2., 100]
    out = []
    for p in percentiles:
        out.append(np.percentile(x, p, **kwargs))
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
        ax.fill_between(logMbinsmid, toplot[1], toplot[5], color=color, facecolor=color, alpha=.2)
        ax.fill_between(logMbinsmid, toplot[2], toplot[4], color=color, facecolor=color, alpha=.2)
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
        ax.fill_between(logMbinsmid, toplot[1], toplot[5], color=color, facecolor=color, alpha=.2)
        ax.fill_between(logMbinsmid, toplot[2], toplot[4], color=color, facecolor=color, alpha=.2)
    ax.legend(loc='lower left')
    ax.set_xlabel(r'$\log M/M_\odot$ at $z_r={}$'.format(zin))
    ax.set_ylabel('fraction merged into z=0 subhalo')
    return fig

if __name__=="__main__":
    zin = 12
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
