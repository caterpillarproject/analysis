import numpy as np
import time,sys
import pandas as pd
import cPickle as pickle
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
#global_cmap = "viridis"
global_cmap = "tab20c_r"

import haloutils
from caterpillaranalysis import MassAccrPlugin

sys.path.append("./greg_dwarfs")
import DwarfMethods as dm

import seaborn as sns
sns.set(context='poster',style='ticks',font='serif',palette='muted')
sns.set_style({"xtick.direction":"in","ytick.direction":"in"})

h0 = .6711
from classify_z8_objects import load_one_halo_data
from classify_z8_objects import allufdtypes, ufdtypescolors, ufdlinestyles
#from classify_z8_objects import logMbins, logVmaxbins, concbins, logDbins, spinbins, TUbins
#from classify_z8_objects import logMbinsmid, logVmaxbinsmid, concbinsmid, logDbinsmid, spinbinsmid, TUbins
from classify_z8_objects import all_bins, all_bins_mid, prop_labels, prop_xlims

#from plot_ufd_surv_figs import allufdtypes, ufdtypescolors, ufdlinestyles
#from plot_ufd_surv_figs import all_bins, all_bins_mid, prop_labels, prop_xlims
from trace_z0_ufds_to_zr import AlexExtantDataPlugin

def plot_2d_hist(zin, which_to_plot, use_phantoms=True):
    ## HACK TO AVOID CIRCULAR REFERENCE
    from plot_ufd_surv_figs import all_bins, all_bins_mid, prop_labels, prop_xlims, prop_cols
    hpaths = dm.get_hpaths(field=False, lx=14)
    all_dfs = []
    for hpath in hpaths:
        hid = haloutils.get_parent_hid(hpath)
        data = pd.DataFrame(np.load("UFDSEARCH_Z0/{}_z{}haloprops.npy".format(
                    haloutils.hidstr(hid), zin)))
        data["hid"] = hid
        if not use_phantoms:
            data = data[data["phantom"] == 0]
        all_dfs.append(data)
    df = pd.concat(all_dfs, ignore_index=True)
    df["logmvir"] = np.log10(df["mvir"]/h0)
    if which_to_plot != "all":
        df = df[df[which_to_plot]]
    
    #cols_to_plot = ["logmvir","T/|U|","spin","conc","logD"]
    cols_to_plot = prop_cols

    Nfigs = len(all_bins)
    fig, axes = plt.subplots(Nfigs, Nfigs, figsize=(8*Nfigs, 8*Nfigs))
    new_bins = [np.arange(4,10,.2), np.arange(0,2.5,.2), np.arange(0,.5,.02),
                np.arange(0,20,.25), np.arange(-1,4,.2)]
    new_bins_mid = [(x[1:]+x[:-1])/2. for x in new_bins]

    #for i, (ybins, ybins_mid, ylabel, ylim) in enumerate(zip(all_bins, all_bins_mid, prop_labels, prop_xlims)):
    for i, (ybins, ybins_mid, ylabel, ylim) in enumerate(zip(new_bins, new_bins_mid, prop_labels, prop_xlims)):
        # First do i == j
        ax = axes[i,i]
        h,x,p = ax.hist(df[cols_to_plot[i]], bins=ybins)
        ax.set_xlabel(ylabel)
        ax.set_xlim(ylim)

        #for j, (xbins, xbins_mid, xlabel, xlim) in enumerate(zip(all_bins, all_bins_mid, prop_labels, prop_xlims)):
        for j, (xbins, xbins_mid, xlabel, xlim) in enumerate(zip(new_bins, new_bins_mid, prop_labels, prop_xlims)):
            if j==i: continue
            #if j>i: ax.set_frame_on(False); continue
            ax = axes[i,j]
            H, xe, ye = np.histogram2d(df[cols_to_plot[j]], df[cols_to_plot[i]],
                                       bins=[xbins,ybins])
            H = np.log10(H)
            im = ax.imshow(H.T, origin="lower", vmin=0, vmax=3.5, cmap=global_cmap,
                           extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]], aspect="auto")
                           #extent=[xlim[0],xlim[1],ylim[0],ylim[1]], aspect="auto")

            ax.set_xlim(xlim)
            ax.set_xlabel(xlabel)
            ax.set_ylim(ylim)
            ax.set_ylabel(ylabel)

    fig.subplots_adjust(right=.85)
    cax = fig.add_axes([.87,.15,.03,.7])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("log N")

    fig.savefig("haloprops_2d_{}_z{}.pdf".format(which_to_plot,zin))
    return fig

def plot_2d_psurv(zin, which_to_plot):
    hpaths = dm.get_hpaths(field=False, lx=14)
    all_dfs = []
    for hpath in hpaths:
        hid = haloutils.get_parent_hid(hpath)
        data = pd.DataFrame(np.load("UFDSEARCH_Z0/{}_z{}haloprops.npy".format(
                    haloutils.hidstr(hid), zin)))
        data["hid"] = hid
        all_dfs.append(data)
    df = pd.concat(all_dfs, ignore_index=True)
    df["logmvir"] = np.log10(df["mvir"]/h0)
    tdf = df[df[which_to_plot]]
    #if which_to_plot != "all":
    #    df = df[df[which_to_plot]]
    
    cols_to_plot = ["logmvir","T/|U|","spin","conc","logD"]

    Nfigs = len(all_bins)
    fig, axes = plt.subplots(Nfigs, Nfigs, figsize=((8*Nfigs)*1.1, 8*Nfigs))
    fig.subplots_adjust(right=.9)
    for i, (ybins, ybins_mid, ylabel, ylim) in enumerate(zip(all_bins, all_bins_mid, prop_labels, prop_xlims)):
        # First do i == j
        ax = axes[i,i]
        h,x,p = ax.hist(df[cols_to_plot[i]], bins=ybins)
        ax.set_xlabel(ylabel)
        ax.set_xlim(ylim)

        for j, (xbins, xbins_mid, xlabel, xlim) in enumerate(zip(all_bins, all_bins_mid, prop_labels, prop_xlims)):
            if j==i: continue
            #if j>i: ax.set_frame_on(False); continue
            ax = axes[i,j]
            Hall, xe, ye = np.histogram2d(df[cols_to_plot[j]], df[cols_to_plot[i]],
                                          bins=[xbins,ybins])
            Hall = Hall.astype(float)
            
            H, xe, ye = np.histogram2d(tdf[cols_to_plot[j]], tdf[cols_to_plot[i]],
                                       bins=[xbins,ybins])
            

            im = ax.imshow(np.log10((H/Hall).T), origin="lower", vmin=-2.5, vmax=0, cmap=global_cmap,
                           extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]], aspect="auto")
                           #extent=[xlim[0],xlim[1],ylim[0],ylim[1]], aspect="auto")

            ax.set_xlim(xlim)
            ax.set_xlabel(xlabel)
            ax.set_ylim(ylim)
            ax.set_ylabel(ylabel)

    fig.subplots_adjust(right=.85)
    cax = fig.add_axes([.87,.15,.03,.7])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("log P({})".format(which_to_plot))

    fig.savefig("halosurv_2d_{}_z{}.pdf".format(which_to_plot,zin))
    return fig

if __name__=="__main__":
    #for zin in [4,6,8,10,12]:
    for zin in [4]:
        #for which_to_plot in ["all","surv","maxm","h14m","h14r","h14i"]:
        for which_to_plot in ["h14m"]:
            start = time.time()
            fig = plot_2d_hist(zin, which_to_plot)
            plt.close(fig)
            if which_to_plot != "all":
                fig = plot_2d_psurv(zin, which_to_plot)
                plt.close(fig)

            print "{} {} {:.1f}".format(zin, which_to_plot, time.time()-start)

def plot_quick():
    zin = 8
    hpaths = dm.get_hpaths(field=False, lx=14)
    hpath = hpaths[0]
    
    out = load_one_halo_data(zin, hpath)
    zrobjs = out[0]
    
    df = pd.DataFrame(zrobjs)
    df['logM'] = out[2]
    df['conc'] = df['rvir']/df['rs']
    df['logD'] = out[5]
    # Manual clipping to increase contrast
    df = df[df['logD'] > 0]
    df = df[df['conc'] <= 20]

    g = sns.PairGrid(df, vars=["logM","T/|U|","conc","spin","logD"],
                     despine=False,size=8)
    g = g.map_diag(plt.hist)
    #g = g.map_offdiag(plt.scatter, marker='.', color='k')
    #g = g.map_lower(sns.kdeplot, cmap="Blues_d")
    g = g.map_offdiag(sns.kdeplot, cmap="Blues_d")
    g.savefig("jointvars_all_kde_z{}.png".format(zin))
    
    plug = AlexExtantDataPlugin()
    df2 = plug.read(hpath)
    zx_ = "z{}_".format(zin)
    df2["logM"] = np.log10(df2[zx_+"mvir"]/h0)
    df2 = df2[df2["logM"] >= 6.5]
    df2["T/|U|"] = df2[zx_+"T/|U|"]
    df2["conc"] = df2[zx_+"rvir"]/df2[zx_+"rs"]
    df2["spin"] = df2[zx_+"spin"]
    print "logM",np.min(df2["logM"]), np.max(df2["logM"])
    print "T/|U|",np.min(df2["T/|U|"]), np.max(df2["T/|U|"])
    print "conc",np.min(df2["conc"]), np.max(df2["conc"])
    print "spin",np.min(df2["spin"]), np.max(df2["spin"])
    
    plug2 = MassAccrPlugin()
    mb = plug2.read(hpath)
    ii = (mb['snap'] == np.unique(df2[zx_+"snap"])[-1]) # avoid -1
    assert np.sum(ii) == 1, "ERROR_A: {} has {} (want {} have {})".format(
        hid,np.sum(ii),np.unique(df2[zx_+"snap"]), mb['snap'])
    hostpos = np.array(mb[ii][['x','y','z']]).view(np.float).reshape(-1,3)
    halopos = df2[[zx_+pos for pos in ['posX','posY','posZ']]].as_matrix()
    df2["logD"] = np.log10(np.sqrt(np.sum((halopos - hostpos)**2,axis=1))) + 3
    print "logD",np.min(df2["logD"]), np.max(df2["logD"])
    
    # Manual cuts
    df2 = df2[df2["T/|U|"] < 5]
    df2 = df2[df2["conc"] < 20]
    df2 = df2[df2["logD"] > 0]
    df2 = df2[df2["spin"] < .3]

    g = sns.PairGrid(df2, vars=["logM","T/|U|","conc","spin","logD"],
                     despine=False,size=8)
    g = g.map_diag(plt.hist)
    #g = g.map_offdiag(plt.scatter, marker='.', color='k')
    g = g.map_offdiag(sns.kdeplot, cmap="Blues_d")
    g.savefig("jointvars_surv_kde_z{}.png".format(zin))
    
    
