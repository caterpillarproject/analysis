import numpy as np
import time,sys
import pandas as pd
import cPickle as pickle
import matplotlib.pyplot as plt

from select_z8_objects import zin_to_zr_snapr
from classify_z8_objects import calculate_total_list, medianscatter, logMbins, logMbinsmid

## Global logMbins
logMbins = np.arange(4,10,.1)
logMbinsmid = (logMbins[1:]+logMbins[:-1])/2.
allufdtypes = ["maxm","h14m","h14r","h14i"]
ufdtypescolors = ['b','g','m','c']
ufdlinestyles = ['-','--',':','-.']

all_zin = [8, 10, 12]
all_colors = [(0.0, 0.4470588235294118, 0.6980392156862745),
              (0.0, 0.6196078431372549, 0.45098039215686275),
              (0.8352941176470589, 0.3686274509803922, 0.0),
              (0.8, 0.4745098039215686, 0.6549019607843137),
              (0.9411764705882353, 0.8941176470588236, 0.25882352941176473),
              (0.33725490196078434, 0.7058823529411765, 0.9137254901960784)]
all_linestyles = ['-','--',':']

import seaborn as sns
sns.set(context='poster',style='ticks',font='serif',palette='colorblind', font_scale=1.5)
sns.set_style({"xtick.direction":"in","ytick.direction":"in"})

if __name__=="__main__":
    all_alldata = []
    for zin in all_zin:
        with open("UFDSEARCH_Z0/histograms_of_z{}.pkl".format(zin),"r") as fp:
            alldata = pickle.load(fp)
        all_alldata.append(alldata)
    
    plotlabelfmt = r"$z_r={}$"

    fig, ax = plt.subplots(figsize=(8,8))
    for zin, alldata, color, ls in zip(all_zin, all_alldata, all_colors, all_linestyles):
        all_total = calculate_total_list(alldata)
        all_hs = []
        all_fracs = []
        for data,total in zip(alldata, all_total):
            h = data[0]
            frac = h.astype(float)/total
            all_hs.append(h)
            all_fracs.append(frac)
        
        toplot = medianscatter(all_fracs, axis=0)
        ax.plot(logMbinsmid, toplot[3], color=color, ls=ls, label=plotlabelfmt.format(zin))
        ax.fill_between(logMbinsmid, toplot[1], toplot[5], color=color, facecolor=color, alpha=.3)
    ax.legend(loc='upper left')
    ax.set_xlabel(r'$\log M/M_\odot$')
    ax.set_ylabel('fraction merged into host')
    ax.set_xlim(6.5,9.0)
    fig.savefig("frac_of_objs_in_host_allz_maxm.pdf", bbox_inches='tight')
    fig.savefig("frac_of_objs_in_host_allz_maxm.png", bbox_inches='tight')

    fig, ax = plt.subplots(figsize=(8,8))
    for zin, alldata, color, ls in zip(all_zin, all_alldata, all_colors, all_linestyles):
        all_total = calculate_total_list(alldata)
        all_hs = []
        all_fracs = []
        for data,total in zip(alldata, all_total):
            h = data[1][0]
            frac = h.astype(float)/total
            all_hs.append(h)
            all_fracs.append(frac)
        
        toplot = medianscatter(all_fracs, axis=0)
        ax.plot(logMbinsmid, toplot[3], color=color, ls=ls, label=plotlabelfmt.format(zin))
        ax.fill_between(logMbinsmid, toplot[1], toplot[5], color=color, facecolor=color, alpha=.3)
    ax.legend(loc='upper left')
    ax.set_xlabel(r'$\log M/M_\odot$')
    ax.set_ylabel('fraction survived as z=0 UFD (maxm)')
    ax.set_xlim(6.5,9.0)
    fig.savefig("frac_of_objs_in_ufds_allz_maxm.pdf", bbox_inches='tight')
    fig.savefig("frac_of_objs_in_ufds_allz_maxm.png", bbox_inches='tight')

    fig, ax = plt.subplots(figsize=(8,8))
    for zin, alldata, color, ls in zip(all_zin, all_alldata, all_colors, all_linestyles):
        all_total = calculate_total_list(alldata)
        all_hs = []
        all_fracs = []
        for data,total in zip(alldata, all_total):
            h = data[2][0]
            frac = h.astype(float)/total
            all_hs.append(h)
            all_fracs.append(frac)
        
        toplot = medianscatter(all_fracs, axis=0)
        ax.plot(logMbinsmid, toplot[3], color=color, ls=ls, label=plotlabelfmt.format(zin))
        ax.fill_between(logMbinsmid, toplot[1], toplot[5], color=color, facecolor=color, alpha=.3)
    ax.legend(loc='upper left')
    ax.set_xlabel(r'$\log M/M_\odot$')
    ax.set_ylabel('fraction merged into z=0 subhalo (maxm)')
    ax.set_xlim(6.5,9.0)
    fig.savefig("frac_of_objs_in_subs_allz_maxm.pdf", bbox_inches='tight')
    fig.savefig("frac_of_objs_in_subs_allz_maxm.png", bbox_inches='tight')

    #plt.show()
