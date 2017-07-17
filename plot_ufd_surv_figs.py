import numpy as np
import time,sys
import pandas as pd
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MultipleLocator

global_cmap="viridis"
global_cmap = "tab20c_r"

h0 = 0.6711

import haloutils
from caterpillaranalysis import MassAccrPlugin

sys.path.append("./greg_dwarfs")
import DwarfMethods as dm

from classify_z8_objects import plot_histograms as plot_zin_histograms
from plot_props_by_massbin import plot_extant_massbin_histograms
from plot_joint_properties import plot_2d_hist, plot_2d_psurv

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
            np.arange(0,.5,.02), np.arange(0,20,.5),
            np.arange(1,4,.1)]
all_bins_mid = [(bins[1:]+bins[:-1])/2. for bins in all_bins]

logmass_bigbins = np.arange(6.5,8.6,0.5)
logmass_bigbins = np.arange(6.5,9.1,0.5)
logmass_bigbinsmid = (logmass_bigbins[1:]+logmass_bigbins[:-1])/2.

global_percentiles = [50-95/2.,50.,50+95/2.]
survtypes = ["surv"]+allufdtypes
fiducial_ufdtype = "maxm"

import seaborn as sns
redshift_colors = sns.color_palette('colorblind')[0:5]
redshift_linestyles = ['-','--','-.',':','-']

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
        try:
            df = groups.get_group(hid)
        except KeyError:
            pass # That host has no halos
        else:
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

def make_fig_1(zin=8, alldata=None):
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
    fig.savefig("fig1_z{}.pdf".format(zin))
    return fig

def make_fig_2(zin=8, alldata=None):
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
    fig.savefig("fig2_z{}.pdf".format(zin))
    return fig

def make_fig_3(zin=8, alldata=None):
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
            try:
                Narr.append(len(groups.get_group(hid)))
            except KeyError:
                Narr.append(0)
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

    fig.savefig("fig3_z{}.pdf".format(zin))
    return fig

def make_fig_4(zin=8, alldata=None):
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
            
            # Compute histogram of all halos
            harr_all = histogram_by_halo(bigbin_data, zin, propcol, bins)
            harr_all = harr_all.astype(float)
            
            for k, (ufdtype, color, ls) in enumerate(zip(allufdtypes, ufdtypescolors, ufdlinestyles)):
                ufd_data = bigbin_data[bigbin_data[ufdtype]]
                harr = histogram_by_halo(ufd_data, zin, propcol, bins)
                to_plot = harr/harr_all
                plot_medscat(ax, binsmid, to_plot, ufdtype, color, ls)
                
            ax.set_ylabel("p(surv ufd | logM, {})".format(prop_short_labels[iprop]))
            ax.set_xlabel(prop_labels[iprop])
            ax.set_title("logM={:.1f}-{:.1f}".format(logmass_bigbins[bigbin_ix], logmass_bigbins[bigbin_ix+1]))
            ax.set_yscale('log')
            ax.set_ylim(1e-3, 1)

    fig.savefig("fig4_z{}.pdf".format(zin))
    return fig

def make_fig_5(zin=8):
    fig = plot_2d_hist(zin, fiducial_ufdtype)
    fig.savefig("fig5_z{}.pdf".format(zin))

def make_fig_6(zin=8, alldata=None):
    #fig = plot_2d_psurv(8, fiducial_ufdtype)
    #fig.savefig("fig6.pdf")
    if alldata is None:
        alldata = load_all_data(zin_to_use=[zin])

    iprops = [1,2,3,4]
    Nprops = len(iprops)
    Nsurvtypes = len(survtypes)

    iyval = 0
    yvals = alldata[prop_cols[iyval]]
    ybins = all_bins[iyval]
    ybinsmid = all_bins_mid[iyval]
    ylabel = prop_labels[iyval]
    ylim = prop_xlims[iyval]

    fig, axes = plt.subplots(Nsurvtypes, Nprops, figsize=(8*Nprops, 8*Nsurvtypes))
    for i,survtype in enumerate(survtypes):
        ii = alldata[survtype]
        for j, iprop in enumerate(iprops):
            ax = axes[i,j]
            xvals = alldata[prop_cols[iprop]]
            xbins = all_bins[iprop]
            xbinsmid = all_bins_mid[iprop]
            xlabel = prop_labels[iprop]
            xlim = prop_xlims[iprop]
            
            Hall, xe, ye = np.histogram2d(xvals, yvals, bins=[xbins,ybins])
            H, xe, ye = np.histogram2d(xvals[ii], yvals[ii], bins=[xbins,ybins])
            Hall = Hall.astype(float)
            
            im = ax.imshow(np.log10(H/Hall).T, origin="lower", vmin=-2.5, vmax=0,
                           cmap=global_cmap,
                           extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]], aspect="auto")
            
            ax.set_title(survtype)
            ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
            ax.set_xlim(xlim); ax.set_ylim(ylim)

    fig.subplots_adjust(right=.85)
    cax = fig.add_axes([.87,.15,.03,.7])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("log P(survive)")
    fig.savefig("fig6_z{}.pdf".format(zin))

def make_fig_7():
    zin_to_use = [4,6,8,10,12]
    alldata = load_all_data(zin_to_use=zin_to_use)
    #fig, axes = plt.subplots(len(zin_to_use), len(prop_cols),
    #                         figsize=(8*len(prop_cols), 8*len(zin_to_use))
    fig, axes = plt.subplots(1, len(prop_cols),
                             figsize=(8*len(prop_cols), 8))
    for i,(zin, color, ls) in enumerate(
        zip(zin_to_use, redshift_colors, redshift_linestyles)):
        for j, (propcol,bins,x,xlabel,slabel,xlim) in enumerate(
            zip(prop_cols, all_bins, all_bins_mid, prop_labels, prop_short_labels, prop_xlims)):
            #ax = axes[i,j]
            ax = axes[j]
            
            allharr = histogram_by_halo(alldata, zin, propcol, bins).astype(float)
            ufdharr = histogram_by_halo(alldata, zin, propcol, bins, survtype=fiducial_ufdtype)
            plot_medscat(ax, x, ufdharr/allharr, "z={}".format(zin), color, ls)

            ax.set_xlabel(xlabel)
            ax.set_xlim(xlim)
            ax.set_yscale('log')
            ax.set_ylabel("P(surv as {} ufd | {})".format(fiducial_ufdtype, slabel))
            ax.set_ylim(1e-3,1)
    axes[2].legend(loc='upper left', ncol=2)
    fig.savefig("fig7.pdf")

def make_fig_8():
    zin_to_use = [4,6,8,10,12]
    allz_alldata = load_all_data(zin_to_use=zin_to_use)

    iprops = [1,2,3,4]
    Nprops = len(iprops)
    Nsurvtypes = len(survtypes)

    iyval = 0
    ybins = all_bins[iyval]
    ybinsmid = all_bins_mid[iyval]
    ylabel = prop_labels[iyval]
    ylim = prop_xlims[iyval]

    fig, axes = plt.subplots(Nsurvtypes, Nprops, figsize=(8*Nprops, 8*Nsurvtypes))
    for i,zin in enumerate(zin_to_use):
        alldata = allz_alldata[allz_alldata["zin"]==zin]
        yvals = alldata[prop_cols[iyval]]
        ii = alldata[fiducial_ufdtype]
        for j, iprop in enumerate(iprops):
            ax = axes[i,j]
            xvals = alldata[prop_cols[iprop]]
            xbins = all_bins[iprop]
            xbinsmid = all_bins_mid[iprop]
            xlabel = prop_labels[iprop]
            xlim = prop_xlims[iprop]
            
            Hall, xe, ye = np.histogram2d(xvals, yvals, bins=[xbins,ybins])
            H, xe, ye = np.histogram2d(xvals[ii], yvals[ii], bins=[xbins,ybins])
            Hall = Hall.astype(float)
            
            im = ax.imshow(np.log10(H/Hall).T, origin="lower", vmin=-2.5, vmax=0,
                           cmap=global_cmap,
                           extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]], aspect="auto")
            
            ax.set_title("z={}".format(zin))
            ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
            ax.set_xlim(xlim); ax.set_ylim(ylim)

    fig.subplots_adjust(right=.85)
    cax = fig.add_axes([.87,.15,.03,.7])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("log P(survive as {})".format(fiducial_ufdtype))
    fig.savefig("fig8.pdf")

def make_fig_9():
    pass

if __name__=="__main__":
    for zin in [4,6,8,10,12]:
        alldatazin = load_all_data(zin_to_use=[zin])
        make_fig_1(zin,alldatazin)
        make_fig_2(zin,alldatazin)
        make_fig_3(zin,alldatazin)
        make_fig_4(zin,alldatazin)
        make_fig_5(zin)
        make_fig_6(zin,alldatazin)
    make_fig_7()
    make_fig_8()
