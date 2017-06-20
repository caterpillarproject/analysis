import numpy as np
import time,sys
import pandas as pd
import cPickle as pickle
import matplotlib.pyplot as plt

import sys
sys.path.append("./greg_dwarfs")
import DwarfMethods as dm
import abundance_matching

import haloutils
import caterpillaranalysis; reload(caterpillaranalysis)

from select_z8_objects import zin_to_zr_snapr
from classify_z8_objects import calculate_total_list, medianscatter, logMbins, logMbinsmid, \
    allufdtypes, ufdtypescolors, ufdlinestyles

all_zin = [8, 10, 12]
all_colors = [(0.0, 0.4470588235294118, 0.6980392156862745),
              (0.0, 0.6196078431372549, 0.45098039215686275),
              (0.8352941176470589, 0.3686274509803922, 0.0),
              (0.8, 0.4745098039215686, 0.6549019607843137),
              (0.9411764705882353, 0.8941176470588236, 0.25882352941176473),
              (0.33725490196078434, 0.7058823529411765, 0.9137254901960784)]
all_linestyles = ['-','--',':']

h0 = .6711

logMpeakbins = np.arange(5,12.5,.2)
logMpeakbinsmid = (logMpeakbins[1:]+logMpeakbins[:-1])/2.

def histogram_objects(zin):
    am = abundance_matching.GarrisonKimmel()
    logLMmin, logLMmax = np.log10(am.stellar_to_halo_mass([1000., 2.e5]))
    plug_mtkey2mpeak = caterpillaranalysis.MaxMassPlugin()

    hpaths = dm.get_hpaths(field=False, lx=14)

    ## Load Fig 1 UFD progenitor histograms
    with open("UFDSEARCH_Z0/summary_data_z{}.pkl".format(zin),"r") as fp:
        all_hids, all_hists_m, all_hists_v, all_hists_maxm, all_num_missing = \
            pickle.load(fp)
    ### Load histograms of objects at
    #with open("UFDSEARCH_Z0/histograms_of_z{}.pkl".format(zin),"r") as fp:
    #    alldata = pickle.load(fp)
    #all_total = calculate_total_list(alldata)
    #all_hs = []
    #all_fracs = []
    #for data,total in zip(alldata, all_total):
    #    h = data[0]
    #    frac = h.astype(float)/total
    #    all_hs.append(h)
    #    all_fracs.append(frac)
    
    all_weighted_hists = []
    all_z0logM_hists = []
    for hpath,all_h_m in zip(hpaths, all_hists_m):
        print hpath
        hid = haloutils.hidstr(haloutils.get_parent_hid(hpath))
        zoomid = haloutils.load_zoomid(hpath)
        
        ## load p(UFD | M_at_zr)
        h_all, h_surv, h_maxm, h_h14m, h_h14r, h_h14i = all_h_m
        all_p_ufd = [h/h_surv.astype(float) for h in [h_maxm, h_h14m, h_h14r, h_h14i]]
        for p_ufd in all_p_ufd: # fix cases where h_surv=0
            p_ufd[np.isnan(p_ufd)] = 0.0

        ## all objects at zr
        zrobjs = np.load("UFDSEARCH_Z0/{}_z{}halos.npy".format(haloutils.hidstr(hid),zin))
        logMzr = np.log10(zrobjs['mvir']/h0)
        ## compute p(UFD) for each object
        indices = np.digitize(logMzr, bins=logMbins)
        over_bin = indices>=len(logMbins)-1
        if np.any(over_bin):
            print "{} over bins (max {:.1f})".format(np.sum(over_bin), np.max(logMzr))
            indices[over_bin] = len(logMbins)-2
        allweights = [p_ufd[indices] for p_ufd in all_p_ufd]

        ## map between z=0 rsid (mtkey) and Mpeak
        ## Greg doesn't have everything so I recomputed it to include small halos
        mtkey2mpeak = plug_mtkey2mpeak.read(hpath)
        logmpeak_of_objs = np.log10(mtkey2mpeak.ix[zrobjs['mtkey']]) #h0 included already
        #print hid, np.sum(pd.isnull(mtkey2mpeak.ix[zrobjs['mtkey']]))
        h_z0_logMpeak, _ = np.histogram(np.log10(mtkey2mpeak), bins = logMpeakbins)
        
        labelfmt = "{} ({}/{})"
        all_Mpeak_of_ufd_hists = []
        fig,ax = plt.subplots(figsize=(8,8))
        fig2,ax2 = plt.subplots(figsize=(8,8))
        for ufdtype, weights, color, ls in zip(allufdtypes, allweights, ufdtypescolors, ufdlinestyles):
            h, _ = np.histogram(logmpeak_of_objs, bins=logMpeakbins, weights=weights)
            ii = np.logical_and(logMpeakbinsmid >=logLMmin, logMpeakbinsmid <= logLMmax)
            Nufd = int(round(np.sum(h[ii]),0))
            Nall = int(round(np.sum(h),0))
            ax.plot(logMpeakbinsmid, h, label=labelfmt.format(ufdtype, Nufd, Nall), drawstyle='steps-mid', color=color, ls=ls)
            all_Mpeak_of_ufd_hists.append(h)

            ax2.plot(logMpeakbinsmid, h/h_z0_logMpeak.astype(float), label=ufdtype, drawstyle='steps-mid', color=color, ls=ls)
        ax.set_xlim(5,12.5)
        ax.set_xlabel(r"Mpeak of objects at z=0")
        ax.set_ylabel(r"Number of UFDs in objects of that Mpeak")
        ax.set_yscale('log')
        ax.set_ylim(.9,1e4)
        ax.add_patch(plt.Rectangle((logLMmin,.9),(logLMmax-logLMmin),(300-.9), color='k', alpha=.2))
        ax.legend(loc='upper left',fontsize=16)
        fig.savefig("UFDSEARCH_Z0/Mpeak_histograms/{}_z{}_Mpeakhists.pdf".format(hid,zin))
        
        ax2.set_xlim(5,12.5)
        ax2.set_xlabel(r"Mpeak of objects at z=0")
        ax2.set_ylabel(r"Average number of UFDs in objects of that Mpeak")
        ax2.set_yscale('log')
        ax2.set_ylim(.9,1e4)
        ax2.add_patch(plt.Rectangle((logLMmin,.9),(logLMmax-logLMmin),(300-.9), color='k', alpha=.2))
        ax2.legend(loc='upper left',fontsize=16)
        fig2.savefig("UFDSEARCH_Z0/Mpeak_histograms/{}_z{}_Mpeakavg.pdf".format(hid,zin))
        
        #with open("UFDSEARCH_Z0/Mpeak_histograms/{}_z{}.pkl".format(hid,zin),"w") as fp:
        #    pickle.dump([logMpeakbins, logMpeakbinsmid, all_Mpeak_of_ufd_hists, h_z0_logMpeak], fp)

        all_weighted_hists.append(all_Mpeak_of_ufd_hists)
        all_z0logM_hists.append(h_z0_logMpeak)
        plt.close(fig)
    with open("UFDSEARCH_Z0/Mpeak_histograms/all_z{}.pkl".format(zin),"w") as fp:
        pickle.dump([logMpeakbins, logMpeakbinsmid, all_weighted_hists, all_z0logM_hists], fp)

def plot_combined(zin):
    fig1, ax1 = plt.subplots(figsize=(8,8))
    fig2, ax2 = plt.subplots(figsize=(8,8))
    with open("UFDSEARCH_Z0/Mpeak_histograms/all_z{}.pkl".format(zin),"r") as fp:
        logMpeakbins, logMpeakbinsmid, all_weighted_hists, all_z0logM_hists = pickle.load(fp)
    hpaths = dm.get_hpaths(field=False, lx=14)
    
    all_counts_of_all_types = [[],[],[],[]] # 4x32xM
    all_avgs_of_all_types = [[],[],[],[]]
    for hpath, all_Mpeak_of_ufd_hists, h_z0_logMpeak in zip(hpaths, all_weighted_hists, all_z0logM_hists):
        for ufdtype, all_counts, all_avgs, h_count in zip(allufdtypes, all_counts_of_all_types,
                                                          all_avgs_of_all_types, all_Mpeak_of_ufd_hists):
            all_counts.append(h_count)
            avg = h_count/h_z0_logMpeak.astype(float)
            avg[np.isnan(avg)] = 0.0
            all_avgs.append(avg)
    
    for ufdtype, all_counts, all_avgs, color, ls in zip(allufdtypes, all_counts_of_all_types,
                                                        all_avgs_of_all_types, ufdtypescolors, ufdlinestyles):
        toplot_counts = medianscatter(np.array(all_counts), axis=0)
        toplot_averages = medianscatter(np.array(all_avgs), axis=0)

        ax1.plot(logMpeakbinsmid, toplot_counts[3], color=color, ls=ls, label=ufdtype)
        ax1.fill_between(logMpeakbinsmid, toplot_counts[1], toplot_counts[5], color=color, facecolor=color, alpha=.2)
        ax2.plot(logMpeakbinsmid, toplot_averages[3], color=color, ls=ls, label=ufdtype)
        ax2.fill_between(logMpeakbinsmid, toplot_averages[1], toplot_averages[5], color=color, facecolor=color, alpha=.2)

    ax1.set_title(r"$z_r = {}$".format(zin))
    ax1.set_xlabel(r'$\log M_{peak}/M_\odot$')
    ax1.set_ylabel('Number')
    ax1.set_xlim(5, 12.5)
    ax1.set_yscale('log')
    ax1.set_ylim(.9,1e4)
    ax1.legend(loc='upper left')
    fig1.savefig("UFDSEARCH_Z0/Mpeak_histograms/all_z{}_Mpeak_counts.pdf".format(zin))

    ax2.set_title(r"$z_r = {}$".format(zin))
    ax2.set_xlabel(r'$\log M_{peak}/M_\odot$')
    ax2.set_ylabel('Average Number per halo')
    ax2.set_xlim(5, 12.5)
    ax2.set_yscale('log')
    ax2.set_ylim(.9,1e4)
    ax2.legend(loc='upper left')
    fig2.savefig("UFDSEARCH_Z0/Mpeak_histograms/all_z{}_Mpeak_avgs.pdf".format(zin))
    return fig1, fig2

if __name__=="__main__":
#    histogram_objects(8)
#    histogram_objects(10)
#    histogram_objects(12)

    plot_combined(8)
    plot_combined(10)
    plot_combined(12)
