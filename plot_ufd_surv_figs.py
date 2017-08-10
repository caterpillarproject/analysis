import numpy as np
import time,sys,os
import pandas as pd
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MultipleLocator

global_cmap="viridis"
global_cmap = "tab20c_r"

h0 = 0.6711

import haloutils

sys.path.append("./greg_dwarfs")
import DwarfMethods as dm
import abundance_matching

from classify_z8_objects import medianscatter
from classify_z8_objects import allufdtypes, ufdtypescolors, ufdlinestyles
from trace_z0_ufds_to_zr import AlexExtantDataPlugin

## Default
#global_use_phantoms=True
#prop_cols = ["logmvir", "T/|U|", "spin", "conc", "logD"]
#prop_short_labels = ["logM", "T/|U|", "spin", "conc", "logD"]
#prop_labels = ["logM (Msun)", "T/|U|", "spin", "conc", "logdist (kpc/h)"]

## Vmaxconc
global_use_phantoms=False
#prop_cols = ["logmvir", "T/|U|", "spin", "vmaxconc", "logD"]
#prop_short_labels = ["logM", "T/|U|", "spin", "vmaxconc", "logD"]
#prop_labels = ["logM (Msun)", "T/|U|", "spin", "vmaxconc", "logdist (kpc/h)"]
prop_cols = ["logmvir", "logeta", "spin", "vmaxconc", "logD"]
prop_short_labels = ["logM", "logeta", "spin", "vmaxconc", "logD"]
prop_labels = ["logM (Msun)", "log 2T/|U|", "spin", "vmaxconc", "logdist (kpc/h)"]

#from classify_z8_objects import all_bins, all_bins_mid
#prop_xlims = [(6.5,9), (0.5,2.5), (0,0.2), (0,20), (1,4)]
#all_bins = [np.arange(6,9,.2), np.arange(.5,2.5,.05),
#            np.arange(0,.5,.02), np.arange(0,20,.5),
#            np.arange(1,4,.1)]
#all_bins_mid = [(bins[1:]+bins[:-1])/2. for bins in all_bins]
prop_xlims = [(6.5,9), (-0.1,1.5), (0,0.2), (0,20), (1,4)]
all_bins = [np.arange(6,9.1,.2), np.arange(-0.1,1.51,.05),
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



def hidstr(hid):
    """ converts halo ID to str Hxxxxxx (copied from haloutils.py) """
    if type(hid)==int or type(hid)==np.int64: return 'H'+str(hid)
    if type(hid)==str:
        if hid[0]=='H': return hid
        return 'H'+hid
    raise ValueError("hid must be int or str, is "+str(type(hid)))

def load_all_hids():
    #hpaths = dm.get_hpaths(field=False, lx=14)
    #hids = [haloutils.get_parent_hid(hpath) for hpath in hpaths]
    hids = \
        [1631506,  264569, 1725139,  447649,    5320,  581141, 1130025, 1387186,  581180,
         1354437, 1725272, 1195448, 1292085,  796175,  388476, 1079897,   94638,   95289,
         1232164, 1422331,  196589, 1268839, 1599988, 1195075, 1631582, 1422429,   65777,
         1232423,  196078, 1599902,  795802, 1104787]
    return hids

def load_all_zoomids():
    #hids = load_all_hids()
    #zoomids = []
    #for hid in hids:
    #    hpath = haloutils.get_hpath_lx(hid, 14)
    #    zoomid = haloutils.load_zoomid(hpath)
    #    zoomids.append(zoomid)
    zoomids = \
        [162271, 222203, 186440, 256463, 161947, 155220, 178802, 150793, 182707, 218773,
         219939,  31921, 217812, 188364, 291771,  53628, 151604, 187355, 145878, 348676,
         167590, 230245, 189407, 216985, 188532, 141967,  68893, 182140, 310128, 200806,
         324544, 168305]        
    return zoomids

def load_all_ufdids():
    hids = load_all_hids()
    ufdids_by_halo = []
    for hid in hids:
        with open("UFDSEARCH_Z0/{}_ufdids.pkl".format(hidstr(hid)),"r") as fp:
            allufdids = pickle.load(fp)
        ufdids_by_halo.append(allufdids)
    return ufdids_by_halo

def load_all_data(zin_to_use = None, hid_to_use = None, use_phantoms=True):
    if zin_to_use is None: zin_to_use = [4,6,8,10,12]

    h0 = 0.6711
    hids = load_all_hids()
    all_dfs = []
    for zin in zin_to_use:
        for hid in hids:
            if hid_to_use is not None:
                if hid not in hid_to_use:
                    continue
            data = pd.DataFrame(np.load("UFDSEARCH_Z0/{}_z{}haloprops.npy".format(
                        hidstr(hid), zin)))
            data["hid"] = hid
            data["zin"] = zin
            if not use_phantoms:
                data = data[data["phantom"]==0]
            all_dfs.append(data)
    alldata = pd.concat(all_dfs, ignore_index=True)
    alldata["logmvir"] = np.log10(alldata["mvir"]/h0)
    alldata["logeta"] = np.log10(2*alldata["T/|U|"])
    return alldata

def classify_halos_at_z0(zin, recalc=False):
    fname = "UFDSEARCH_Z0/classification_z{}.npy".format(zin)
    if os.path.exists(fname) and not recalc:
        # load the data
        classification = np.load(fname)
        return classification
    
    
    alldata = load_all_data(zin_to_use=[zin], use_phantoms=global_use_phantoms)
    
    am = abundance_matching.GarrisonKimmel()
    LMmin, LMmax = am.stellar_to_halo_mass([1000., 2.e5])
    
    hids = load_all_hids()
    zoomids = load_all_zoomids()
    ufdids_by_halo = load_all_ufdids()
    
    # Classificationify each halo into what object it ends up in
    # host, classical dSph, four ufd types
    classification = np.full((len(alldata), 6), False, dtype=bool)

    start = time.time()
    ii_big = np.array(alldata["extant_maxmass"] > LMmax)
    for hid, zoomid, allufdids in zip(hids, zoomids, ufdids_by_halo):
        ii_this = np.array(alldata["hid"] == hid)
        ii_host = np.array(alldata["mtkey"] == zoomid) & ii_this
        classification[ii_host, 0] = True

        ii_cdsph = ii_big & (~ii_host) & ii_this
        classification[:,1][ii_cdsph] = True
        
        print hid,0,np.sum(ii_host)
        print hid,1,np.sum(ii_cdsph)

        for j, ufdids in enumerate(allufdids):
            ufdids = np.array(ufdids).astype(int)
            col = j+2
            _ii_ufd = np.array(map(lambda x: x in ufdids, np.array(alldata["mtkey"])))
            ii_ufd = np.logical_and(ii_this, _ii_ufd)
            print hid,col,np.sum(ii_ufd)
            classification[ii_ufd, col] = True

    print time.time()-start
    print "Saving to {}".format(fname)
    np.save(fname, classification)
    return classification

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
            vals_to_hist = np.array(df[propcol])
            vals_to_hist = vals_to_hist[np.isfinite(vals_to_hist)]
            if len(vals_to_hist) != len(df):
                print "halo {} prop {} has {} bad values".format(hid, propcol, len(df)-len(vals_to_hist))
            h,x = np.histogram(vals_to_hist, bins=bins)
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
        vals_to_hist1 = np.array(df[propcol1]); vals_to_hist1 = vals_to_hist1[np.isfinite(vals_to_hist1)]
        vals_to_hist2 = np.array(df[propcol2]); vals_to_hist2 = vals_to_hist2[np.isfinite(vals_to_hist2)]
        H,x,y = np.histogram2d(vals_to_hist1, vals_to_hist2, bins=[bins1,bins2])
        output[i,:,:] = H
    return output

def plot_medscat(ax, x, yarr, label, color, ls):
    y1,y2,y3 = medianscatter(yarr, percentiles=global_percentiles, axis=0)
    ax.plot(x, y2, drawstyle='steps-mid', color=color, ls=ls, label=label)
    ax.fill_between(x, y1, y3, color=color, facecolor=color, alpha=.3, step='mid')

def make_fig_1(zin=8, alldata=None):
    if alldata is None:
        alldata = load_all_data(zin_to_use=[zin], use_phantoms=global_use_phantoms)
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
        alldata = load_all_data(zin_to_use=[zin], use_phantoms=global_use_phantoms)
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
        alldata = load_all_data(zin_to_use=[zin], use_phantoms=global_use_phantoms)
    
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
        alldata = load_all_data(zin_to_use=[zin], use_phantoms=global_use_phantoms)
    
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

def make_fig_5(zin=8, alldata=None, subs_only=False, hosts_only=False):
    #fig = plot_2d_hist(zin, fiducial_ufdtype, use_phantoms=global_use_phantoms)
    if alldata is None:
        alldata = load_all_data(zin_to_use=[zin], use_phantoms=global_use_phantoms)

    if subs_only and hosts_only:
        raise ValueError("Must specify only one of subs or hosts only")
    if subs_only:
        alldata = alldata[alldata["pid"] != -1]
        figname_extra = "_subs"
    elif hosts_only:
        alldata = alldata[alldata["pid"] == -1]
        figname_extra = "_hosts"
    else:
        figname_extra = ""
    
    # remove nans
    ii_to_remove = np.zeros(len(alldata), dtype=bool)
    for col in prop_cols:
        ii_to_remove = ii_to_remove | np.isnan(alldata[col])
    if np.sum(ii_to_remove) > 0:
        print "Removing {} bad points".format(np.sum(ii_to_remove))
        alldata = alldata[~ii_to_remove]

    Nfigs = len(all_bins)
    fig, axes = plt.subplots(Nfigs, Nfigs, figsize=(8*Nfigs, 8*Nfigs))
    for i, (ybins, ybins_mid, ylabel, ylim) in enumerate(zip(all_bins, all_bins_mid, prop_labels, prop_xlims)):
        # First do i == j
        ax = axes[i,i]
        h,x,p = ax.hist(alldata[prop_cols[i]], bins=ybins)
        ax.set_xlabel(ylabel)
        ax.set_xlim(ylim)
        for j, (xbins, xbins_mid, xlabel, xlim) in enumerate(zip(all_bins, all_bins_mid, prop_labels, prop_xlims)):
            if j==i: continue
            ax = axes[i,j]
            H, xe, ye = np.histogram2d(alldata[prop_cols[j]], alldata[prop_cols[i]],
                                       bins=[xbins,ybins])
            H = np.log10(H)
            im = ax.imshow(H.T, origin="lower", vmin=0, vmax=3.5, cmap=global_cmap,
                           extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]], aspect="auto")

            ax.set_xlim(xlim)
            ax.set_xlabel(xlabel)
            ax.set_ylim(ylim)
            ax.set_ylabel(ylabel)
        
    fig.subplots_adjust(right=.85)
    cax = fig.add_axes([.87,.15,.03,.7])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("log N")

    fig.savefig("fig5_z{}{}.pdf".format(zin,figname_extra))

def make_fig_6(zin=8, alldata=None, subs_only=False, hosts_only=False):
    #fig = plot_2d_psurv(8, fiducial_ufdtype)
    #fig.savefig("fig6.pdf")
    if alldata is None:
        alldata = load_all_data(zin_to_use=[zin], use_phantoms=global_use_phantoms)

    if subs_only and hosts_only:
        raise ValueError("Must specify only one of subs or hosts only")
    if subs_only:
        alldata = alldata[alldata["pid"] != -1]
        figname_extra = "_subs"
    elif hosts_only:
        alldata = alldata[alldata["pid"] == -1]
        figname_extra = "_hosts"
    else:
        figname_extra = ""
           

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
    fig.savefig("fig6_z{}{}.pdf".format(zin,figname_extra))

def make_fig_7():
    zin_to_use = [4,6,8,10,12]
    alldata = load_all_data(zin_to_use=zin_to_use, use_phantoms=global_use_phantoms)
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
    allz_alldata = load_all_data(zin_to_use=zin_to_use, use_phantoms=global_use_phantoms)

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

def make_fig_9(zin=8, alldata=None):
    """ This is actually Fig 2/3 """
    if alldata is None:
        alldata = load_all_data(zin_to_use=[zin], use_phantoms=global_use_phantoms)
    # 0: host
    # 1: cdsph
    # 2-5: maxm, h14m, h14r, h14i
    propcol = prop_cols[0] #logmvir
    x = all_bins_mid[0]
    xlabel = prop_labels[0]
    xlim = prop_xlims[0]
    xbins = all_bins[0]

    classification = classify_halos_at_z0(zin)
    
    fig, axes = plt.subplots(1,4, figsize=(8*4, 8*1))
    allharr = histogram_by_halo(alldata, zin, propcol, xbins)
    
    ax = axes[0]
    harr = histogram_by_halo(alldata[classification[:,0]], zin, propcol, xbins)
    plot_medscat(ax, x, harr/allharr, "host", 'k', '-')
    
    ax = axes[1]
    harr = histogram_by_halo(alldata[classification[:,1]], zin, propcol, xbins)
    plot_medscat(ax, x, harr/allharr, "classical dSph", 'k', '-')
    
    for j, (ufdtype, color, ls) in enumerate(zip(allufdtypes, ufdtypescolors, ufdlinestyles)):
        classcol = j+2
        
        # UFD
        iiufdsurv = np.array(alldata[ufdtype])
        #harr = histogram_by_halo(alldata[classification[:,classcol]], zin, propcol, xbins)
        harr = histogram_by_halo(alldata[iiufdsurv], zin, propcol, xbins)
        plot_medscat(axes[2], x, harr/allharr, ufdtype, color, ls)
        
        # Not UFD, host, or cdSph
        #ii = (np.sum(classification[:,[0,1,classcol]], axis=1) == 0)
        #harr = histogram_by_halo(alldata[ii], zin, propcol, xbins)
        #plot_medscat(axes[3], x, harr/allharr, "not "+ufdtype, color, ls)
        
        # In UFD but not main branch progenitor
        iiufddead = np.logical_and(classification[:,classcol], ~iiufdsurv)
        harr = histogram_by_halo(alldata[iiufddead], zin, propcol, xbins)
        plot_medscat(axes[3], x, harr/allharr, "merged into "+ufdtype, color, ls)
        
    for ax in axes:
        ax.set_xlim(xlim)
        ax.set_ylim(0,1)
        ax.set_xlabel(xlabel)
        ax.legend(loc="upper right", fontsize=12)
        ax.text(.05,.9,"z={}".format(zin), ha='left', transform=ax.transAxes)
    axes[0].set_ylabel("p(merge into host | M)")
    axes[1].set_ylabel("p(merge into classical dSph | M)")
    axes[2].set_ylabel("p(main progenitor of UFD | M)")
    axes[3].set_ylabel("p(other progenitor of UFD | M)")
    fig.savefig("fig9_z{}.pdf".format(zin))

def make_fig_10(zin=8, alldata=None):
    """ This is actually Fig 3 """
    zin_to_use = [4,6,8,10,12]
    #allz_alldata = load_all_data(zin_to_use=zin_to_use, use_phantoms=global_use_phantoms)

    # 0: host
    # 1: cdsph
    # 2-5: maxm, h14m, h14r, h14i
    propcol = prop_cols[0] #logmvir
    x = all_bins_mid[0]
    xlabel = prop_labels[0]
    xlim = prop_xlims[0]
    xbins = all_bins[0]

    fig, axes = plt.subplots(1,4, figsize=(8*4, 8*1))
    for i, (zin, color, ls) in enumerate(zip(zin_to_use, redshift_colors, redshift_linestyles)):
        alldata = load_all_data(zin_to_use=[zin], use_phantoms=global_use_phantoms)
        classification = classify_halos_at_z0(zin)
        allharr = histogram_by_halo(alldata, zin, propcol, xbins)
        
        ax = axes[0]
        harr = histogram_by_halo(alldata[classification[:,0]], zin, propcol, xbins)
        plot_medscat(ax, x, harr/allharr, "z={}".format(zin), color, ls)
        
        ax = axes[1]
        harr = histogram_by_halo(alldata[classification[:,1]], zin, propcol, xbins)
        plot_medscat(ax, x, harr/allharr, "z={}".format(zin), color, ls)
        
        ufdtype = "maxm"
        classcol = 2
        
        # UFD
        iiufdsurv = np.array(alldata[ufdtype])
        #harr = histogram_by_halo(alldata[classification[:,classcol]], zin, propcol, xbins)
        harr = histogram_by_halo(alldata[iiufdsurv], zin, propcol, xbins)
        plot_medscat(axes[2], x, harr/allharr, "{} z={}".format(ufdtype,zin), color, ls)
        
        # In UFD but not main branch progenitor
        iiufddead = np.logical_and(classification[:,classcol], ~iiufdsurv)
        harr = histogram_by_halo(alldata[iiufddead], zin, propcol, xbins)
        plot_medscat(axes[3], x, harr/allharr, "{} z={}".format(ufdtype,zin), color, ls)
        
    for ax in axes:
        ax.set_xlim(xlim)
        ax.set_ylim(0,1)
        ax.set_xlabel(xlabel)
        ax.legend(loc="upper left", fontsize=12)
        #ax.text(.05,.9,"z={}".format(zin), ha='left', transform=ax.transAxes)
    axes[0].set_ylabel("p(merge into host | M)")
    axes[1].set_ylabel("p(merge into classical dSph | M)")
    axes[2].set_ylabel("p(main progenitor of UFD | M)")
    axes[3].set_ylabel("p(other progenitor of UFD | M)")
    fig.savefig("fig10.pdf")

if __name__=="__main__":
    make_fig_10()
    #for zin in [4,6,8,10,12]:
    #    make_fig_9(zin)
    #    alldatazin = load_all_data(zin_to_use=[zin], use_phantoms=global_use_phantoms)
    #    make_fig_1(zin,alldatazin)
    #    make_fig_2(zin,alldatazin)
    #    make_fig_3(zin,alldatazin)
    #    make_fig_4(zin,alldatazin)
    #    make_fig_5(zin,alldatazin)
    #    make_fig_5(zin,alldatazin,subs_only=True)
    #    make_fig_5(zin,alldatazin,hosts_only=True)
    #    make_fig_6(zin,alldatazin)
    #    make_fig_6(zin,alldatazin,subs_only=True)
    #    make_fig_6(zin,alldatazin,hosts_only=True)
    #    plt.close("all")
    #make_fig_7()
    #make_fig_8()
    plt.close("all")
