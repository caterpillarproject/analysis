import numpy as np
import haloutils
import time,sys
import pandas as pd
import matplotlib.pyplot as plt
import cPickle as pickle

from select_z8_objects import zin_to_zr_snapr
from classify_z8_objects import logMbins, logMbinsmid, allufdtypes, \
    ufdtypescolors, ufdlinestyles, medianscatter
    
import seaborn as sns
sns.set(context='poster',style='ticks',font='serif',palette='colorblind', font_scale=1.5)
sns.set_style({"xtick.direction":"in","ytick.direction":"in"})

def plot_total(zin):
    with open("UFDSEARCH_Z0/summary_data_z{}.pkl".format(zin),"r") as fp:
        all_hids, all_hists_m, all_hists_v, all_hists_maxm, all_num_missing = pickle.load(fp)
    
    fig, ax = plt.subplots(figsize=(8,8))
    Marr = np.array(all_hists_m).astype(float)
    logMbinsmid = (logMbins[1:]+logMbins[:-1])/2.
    colors = ufdtypescolors
    labels = allufdtypes
    for j in range(4):
        med_max = np.percentile(Marr[:,j+2,:]/Marr[:,1,:], 100., axis=0)
        med_p2s = np.percentile(Marr[:,j+2,:]/Marr[:,1,:], 50+95/2., axis=0)
        med_p1s = np.percentile(Marr[:,j+2,:]/Marr[:,1,:], 50+68/2., axis=0)
        median  = np.percentile(Marr[:,j+2,:]/Marr[:,1,:], 50, axis=0)
        med_m1s = np.percentile(Marr[:,j+2,:]/Marr[:,1,:], 50-68/2., axis=0)
        med_m2s = np.percentile(Marr[:,j+2,:]/Marr[:,1,:], 50-95/2., axis=0)
        med_min = np.percentile(Marr[:,j+2,:]/Marr[:,1,:], 0., axis=0)
        ax.plot(logMbinsmid, median, color=colors[j], drawstyle='steps-mid', label=labels[j],
                ls=ufdlinestyles[j])
        #ax.fill_between(logMbinsmid, med_m1s, med_p1s, color=colors[j], facecolor=colors[j], alpha=.5)
        ax.fill_between(logMbinsmid, med_m2s, med_p2s, color=colors[j], facecolor=colors[j], alpha=.3, step='mid')
        #ax.fill_between(logMbinsmid, med_min, med_max, color=colors[j], facecolor=colors[j], alpha=.05)
        #ax.plot(logMbinsmid, med_min, ':', color=colors[j], lw=.5, alpha=1)
        #ax.plot(logMbinsmid, med_max, ':', color=colors[j], lw=.5, alpha=1)

    ax.set_xlim(5.5,8.5)
    ax.set_ylim(0,1)
    ax.legend(loc='upper left', fontsize=16)
    #ax.set_title('$z_r$={}'.format(zin))
    ax.set_xlabel(r'$\log M/M_\odot$ at $z_r$={}'.format(zin))
    ax.set_ylabel('fraction')
    #fig.savefig("UFDSEARCH_Z0/total_z{}.png".format(zin), bbox_inches='tight')
    fig.savefig("UFDSEARCH_Z0/total_z{}.pdf".format(zin), bbox_inches='tight')
        

def plot_hists_one(zin, halo_ix=0):
    fig,ax = plt.subplots(figsize=(8,8))
    
    with open("UFDSEARCH_Z0/summary_data_z{}.pkl".format(zin),"r") as fp:
        all_hids, all_hists_m, all_hists_v, all_hists_maxm, all_num_missing = pickle.load(fp)
    
    h_all, h_surv, h_maxm, h_h14m, h_h14r, h_h14i = all_hists_m[halo_ix]
    ax.plot(logMbinsmid, h_all,  label='all',  drawstyle='steps-mid', color='k', ls=':')
    ax.plot(logMbinsmid, h_surv, label='surv', drawstyle='steps-mid', color='k')
    h_to_plot = [h_maxm, h_h14m, h_h14r, h_h14i]
    for j in range(4):
        ax.plot(logMbinsmid, h_to_plot[j], label=allufdtypes[j], drawstyle='steps-mid',
                ls=ufdlinestyles[j], color=ufdtypescolors[j])

    ax.set_xlabel(r'$\log M/M_\odot$ at $z_r$={}'.format(zin))
    #ax.set_xlabel("M200c at z={}".format(zin))
    ax.set_yscale('log')
    ax.set_ylabel('number')
    #ax.text(.95,.93,haloutils.hidstr(hid), transform=ax.transAxes, ha='right')
    #ax.text(.95,.87,"{:.2e}".format(float(haloutils.load_haloprops(hpath)[0])), transform=ax.transAxes, ha='right')
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.add_patch(plt.Rectangle((xlim[0],ylim[0]),6-xlim[0],ylim[1]-ylim[0], alpha=.2, color='k'))
    fig.savefig("UFDSEARCH_Z0/ufd_histogram_{:02}.pdf".format(halo_ix))

if __name__=="__main__":
    plot_total(8)
    plot_total(10)
    plot_total(12)
    plot_hists_one(8)

