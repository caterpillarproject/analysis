import numpy as np
import time,sys
import pandas as pd
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MultipleLocator

import haloutils
from caterpillaranalysis import MassAccrPlugin

sys.path.append("./greg_dwarfs")
import DwarfMethods as dm

from classify_z8_objects import plot_histograms as plot_zin_histograms
from plot_props_by_massbin import plot_extant_massbin_histograms

from classify_z8_objects import medianscatter
from classify_z8_objects import allufdtypes, ufdtypescolors, ufdlinestyles
from classify_z8_objects import logMbins, logVmaxbins, concbins, logDbins, spinbins, TUbins
from classify_z8_objects import logMbinsmid, logVmaxbinsmid, concbinsmid, logDbinsmid, spinbinsmid, TUbins
from trace_z0_ufds_to_zr import AlexExtantDataPlugin

prop_cols = ["logmvir", "T/|U|", "spin", "conc", "logD"]
prop_short_labels = ["logM", "T/|U|", "spin", "conc", "logD"]
prop_labels = ["logM (Msun)", "T/|U|", "spin", "conc", "logdist (kpc/h)"]
prop_xlims = [(6.5,9), (0.5,2.5), (0,0.2), (0,20), (1,4)]
from classify_z8_objects import all_bins, all_bins_mid
all_bins = [np.arange(6,9,.2), np.arange(.5,2.5,.05),
            np.arange(0,.5,.02), np.arange(0,20,.2),
            np.arange(1,4,.1)]
all_bins_mid = [(bins[1:]+bins[:-1])/2. for bins in all_bins]

logmass_bigbins = np.arange(6.5,8.6,0.5)
logmass_bigbinsmid = (logmass_bigbins[1:]+logmass_bigbins[:-1])/2.

global_percentiles = [50-95/2.,50.,50+95/2.]
survtypes = ["surv"]+allufdtypes

import seaborn as sns
sns.set(context='poster',style='ticks',font='serif',palette='muted')
sns.set_style({"xtick.direction":"in","ytick.direction":"in"})



def load_all_hids():
    hpaths = dm.get_hpaths(field=False, lx=14)
    hids = [haloutils.get_parent_hid(hpath) for hpath in hpaths]
    return hids

def load_all_data(zin_to_use = None, hid_to_use = None):
    if zin_to_use is None: zin_to_use = [4,6,8,10,12]

    h0 = 0.6711
    hpaths = dm.get_hpaths(field=False, lx=14)
    all_dfs = []
    for zin in zin_to_use:
        for hpath in hpaths:
            hid = haloutils.get_parent_hid(hpath)
            if hid_to_use is not None:
                if hid not in hid_to_use:
                    continue
            data = pd.DataFrame(np.load("UFDSEARCH_Z0/{}_z{}haloprops.npy".format(
                        haloutils.hidstr(hid), zin)))
            data["hid"] = hid
            data["zin"] = zin
            all_dfs.append(data)
    alldata = pd.concat(all_dfs, ignore_index=True)
    alldata["logmvir"] = np.log10(alldata["mvir"]/h0)
    return alldata

def histogram_by_halo(alldata, zin, propcol, bins, survtype=None):
    """
    Create a Nhalo x Nbin array histogramming propcol
    """
    assert propcol in alldata.columns, propcol
    data = alldata[alldata["zin"]==zin]
    if survtype is not None:
        data = data[data[survtype]]
    groups = data.groupby(["hid"], sort=False)
    hids = load_all_hids()
    
    output = np.zeros((len(hids), len(bins)-1))
    for i,hid in enumerate(hids):
        df = groups.get_group(hid)
        h,x = np.histogram(df[propcol], bins=bins)
        output[i,:] = h
    return output

def histogram_2d_by_halo(alldata, zin, propcol1, propcol2, bins1, bins2):
    """
    Create a Nhalo x Nbin1 x Nbin2 array histogramming propcol1/2
    """
    assert propcol1 in alldata.columns, propcol1
    assert propcol2 in alldata.columns, propcol2
    data = alldata[alldata["zin"]==zin]
    groups = data.groupby(["hid"], sort=False)
    hids = load_all_hids()
    
    output = np.zeros((len(hids), len(bins1)-1, len(bins2)-1))
    for i,hid in enumerate(hids):
        df = groups.groupby(hid)
        H,x,y = np.histogram2d(df[propcol1], df[propcol2], bins=[bins1,bins2])
        output[i,:,:] = H
    return output

def plot_medscat(ax, x, yarr, label, color, ls):
    y1,y2,y3 = medianscatter(yarr, percentiles=global_percentiles, axis=0)
    ax.plot(x, y2, drawstyle='steps-mid', color=color, ls=ls, label=label)
    ax.fill_between(x, y1, y3, color=color, facecolor=color, alpha=.3, step='mid')

def make_fig_1(alldata=None):
    zin = 8
    if alldata is None:
        alldata = load_all_data(zin_to_use=[zin])
    fig, axes = plt.subplots(2,3, figsize=(8*3, 8*2))
    for propcol,ax,bins,x,xlabel,slabel,xlim in \
            zip(prop_cols, axes.flat[:-1], all_bins, all_bins_mid, 
                prop_labels, prop_short_labels, prop_xlims):
        harr = histogram_by_halo(alldata, zin, propcol, bins)
        plot_medscat(ax, x, harr, "All halos", 'k', '-')
        harr = histogram_by_halo(alldata, zin, propcol, bins, survtype="surv")
        plot_medscat(ax, x, harr, "surv", 'gray', ':')
        for ufdtype,color,ls in zip(allufdtypes,ufdtypescolors,ufdlinestyles):
            harr = histogram_by_halo(alldata, zin, propcol, bins, survtype=ufdtype)
            plot_medscat(ax, x, harr, ufdtype, color, ls)
        ax.set_xlabel(xlabel)
        ax.set_xlim(xlim)
        ax.set_yscale('log')
        ax.set_ylabel("N")
    ax = axes[0,0]
    handles,labels = ax.get_legend_handles_labels()
    lax = axes[1,2]
    lax.legend(handles,labels,loc='center')
    lax.text(.5,.8,"z={}".format(zin),color='k',transform=lax.transAxes,ha='center')
    lax.set_xticks([]); lax.set_yticks([])
    fig.savefig("fig1.pdf")
    return fig

def make_fig_2(alldata=None):
    zin = 8
    if alldata is None:
        alldata = load_all_data(zin_to_use=[zin])
    fig, axes = plt.subplots(2,3, figsize=(8*3, 8*2))
    for propcol,ax,bins,x,xlabel,slabel,xlim in \
            zip(prop_cols, axes.flat[:-1], all_bins, all_bins_mid, 
                prop_labels, prop_short_labels, prop_xlims):
        allharr = histogram_by_halo(alldata, zin, propcol, bins)
        allharr = allharr.astype(float)
        #plot_medscat(ax, x, harr, "All halos", 'k', '-')
        harr = histogram_by_halo(alldata, zin, propcol, bins, survtype="surv")
        plot_medscat(ax, x, harr/allharr, "surv", 'gray', ':')
        for ufdtype,color,ls in zip(allufdtypes,ufdtypescolors,ufdlinestyles):
            harr = histogram_by_halo(alldata, zin, propcol, bins, survtype=ufdtype)
            plot_medscat(ax, x, harr/allharr, ufdtype, color, ls)
        ax.set_xlabel(xlabel)
        ax.set_xlim(xlim)
        ax.set_yscale('log')
        ax.set_ylabel("P(surv | {})".format(slabel))
    ax = axes[0,0]
    handles,labels = ax.get_legend_handles_labels()
    lax = axes[1,2]
    lax.legend(handles,labels,loc='center')
    lax.text(.5,.8,"z={}".format(zin),color='k',transform=lax.transAxes,ha='center')
    lax.set_xticks([]); lax.set_yticks([])
    fig.savefig("fig2.pdf")
    return fig

def make_fig_3(alldata=None):
    zin = 8
    if alldata is None:
        alldata = load_all_data(zin_to_use=[zin])
    
    Nbigbins = len(logmass_bigbinsmid)
    iprops = [1,2,3,4]
    Nprops = len(iprops)
    
    fig, axes = plt.subplots(Nbigbins, Nprops, figsize=(8*Nprops, 8*Nbigbins))

    def _N_by_halo(df):
        hids = load_all_hids()
        groups = df.groupby("hid")
        Narr = []
        for hid in hids:
            Narr.append(len(groups.get_group(hid)))
        return np.array(Narr)

    bigbin_indices = np.digitize(alldata["logmvir"], bins=logmass_bigbins)
    for bigbin_ix in range(Nbigbins):
        bigbin_ii = bigbin_indices == (bigbin_ix + 1)
        bigbin_data = alldata[bigbin_ii]
        for j,iprop in enumerate(iprops):
            ax = axes[bigbin_ix, j]
            bins = all_bins[iprop]
            binsmid = all_bins_mid[iprop]
            propcol = prop_cols[iprop]
            dx = bins[1]-bins[0]
            
            # Compute normalized PDF of all halos in bigbin
            harr = histogram_by_halo(bigbin_data, zin, propcol, bins)
            Narr = _N_by_halo(bigbin_data).astype(float)
            harr = (harr.T/Narr).T
            harr /= dx # form the PDF
            plot_medscat(ax, binsmid, harr, "All halos", 'k', '-')
            
            for k, (ufdtype, color, ls) in enumerate(zip(allufdtypes, ufdtypescolors, ufdlinestyles)):
                ufd_data = bigbin_data[bigbin_data[ufdtype]]
                harr = histogram_by_halo(ufd_data, zin, propcol, bins)
                Narr = _N_by_halo(ufd_data).astype(float)
                harr = (harr.T/Narr).T
                harr /= dx # form the PDF
                plot_medscat(ax, binsmid, harr, ufdtype, color, ls)
                
            ax.set_ylabel("p({} | logM)".format(prop_short_labels[iprop]))
            ax.set_xlabel(prop_labels[iprop])
            ax.set_title("logM={:.1f}-{:.1f}".format(logmass_bigbins[bigbin_ix], logmass_bigbins[bigbin_ix+1]))
            ax.set_yscale('log')

    fig.savefig("fig3.pdf")
    return fig
    #fig = plot_extant_massbin_histograms(8, plot_type="count", log=False)
    #fig.savefig("fig3.pdf",bbox_inches="tight")
    #plt.close(fig)

def make_fig_4():
    fig = plot_extant_massbin_histograms(8, plot_type="ratio", log=True)
    fig.savefig("fig4.pdf",bbox_inches="tight")
    plt.close(fig)

def make_fig_5():
    pass

def make_fig_6():
    pass

def make_fig_7():
    pass

def make_fig_8():
    pass

def make_fig_9():
    pass

if __name__=="__main__":
    alldataz8 = load_all_data(zin_to_use=[8])
    #make_fig_1(alldataz8)
    #make_fig_2(alldataz8)
    make_fig_3(alldataz8)

