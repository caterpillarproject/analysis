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

import seaborn as sns
sns.set(context='poster',style='ticks',font='serif',palette='muted')
sns.set_style({"xtick.direction":"in","ytick.direction":"in"})

from select_z8_objects import zin_to_zr_snapr
from classify_z8_objects import medianscatter, load_one_halo_data #, histogram_one_halo
from classify_z8_objects import logMbins, logMbinsmid, allufdtypes, ufdtypescolors, ufdlinestyles
from classify_z8_objects import logVmaxbins, logVmaxbinsmid, concbins, concbinsmid
from classify_z8_objects import logDbins, logDbinsmid, spinbins, spinbinsmid
from classify_z8_objects import all_bins, all_bins_mid, prop_labels, prop_xlims
from trace_z0_ufds_to_zr import AlexExtantDataPlugin

#logmass_bigbins = np.arange(6.5,10.1,0.5)
logmass_bigbins = np.arange(6.5,9.1,0.5)
logmass_bigbinsmid = (logmass_bigbins[1:]+logmass_bigbins[:-1])/2.

h0=.6711

"""
def histogram_halos(zin):
    hpaths = dm.get_hpaths(field=False, lx=14)

    Nhalos = len(hpaths)
    Nbigbins = len(logmass_bigbinsmid)
    Nprops = 4

    allpdfs = [np.zeros((Nhalos, Nbigbins, len(all_bins_mid[j+1]))) for j in range(Nprops)]
    allufdpdfs = [[np.zeros((Nhalos, Nbigbins, len(all_bins_mid[j+1]))) \
                       for k in range(len(allufdtypes))] for j in range(Nprops)]
    for i, hpath in enumerate(hpaths):
        start = time.time()
        # Load data
        hid = haloutils.get_parent_hid(hpath)
        zrobjs, rscatobjs, zrlogmass, zrlogvmax, zrspin, zrlogD, conc = load_one_halo_data(zin, hpath)
        zrTU = zrobjs["T/|U|"]
        with open("UFDSEARCH_Z0/{}_ufdids.pkl".format(haloutils.hidstr(hid)),"r") as fp:
            allufdids = pickle.load(fp)
        
        # Make histograms
        props_to_hist = [zrTU, zrspin, zrlogD, conc]
        assert len(props_to_hist) == Nprops
        bigbin_indices = np.digitize(zrlogmass, bins=logmass_bigbins)
        for bigbin_ix in range(Nbigbins):
            bigbin_ii = bigbin_indices == (bigbin_ix + 1)
            for j in range(Nprops):
                prop_ix = j+1 # skip mvir
                bins = all_bins[prop_ix]
                binsmid = all_bins_mid[prop_ix]
    
                ii = bigbin_ii
                x = props_to_hist[j][ii]
                x = x[np.isfinite(x)]
                h, _ = np.histogram(x, bins=bins)
                y = h/float(len(x))
                dx = bins[1]-bins[0]
                y /= dx
                allpdfs[j][i,bigbin_ix,:] = y
                
                for k, (ufdids, ufdtype) in enumerate(zip(allufdids, allufdtypes)):
                    ufdids = np.array(ufdids).astype(int)
                    ii_ufd = np.array(map(lambda x: x in ufdids, zrobjs['mtkey']))
                    
                    ii = bigbin_ii & ii_ufd
                    x = props_to_hist[j][ii]
                    x = x[np.isfinite(x)]
                    h, _ = np.histogram(x, bins=bins)
                    y = h/float(len(x))
                    y /= dx
    
                    allufdpdfs[j][k][i,bigbin_ix,:] = y
        print "{} took {:.1f}".format(hid, time.time()-start)
    
    with open("UFDSEARCH_Z0/hist_by_mass_bin_z{}.pkl".format(zin),"w") as fp:
        pickle.dump([allpdfs, allufdpdfs], fp)

def plot_massbin_histograms(zin):
    with open("UFDSEARCH_Z0/hist_by_mass_bin_z{}.pkl".format(zin),"r") as fp:
        allpdfs, allufdpdfs = pickle.load(fp)

    hpaths = dm.get_hpaths(field=False, lx=14)
    Nhalos = len(hpaths)
    Nbigbins = len(logmass_bigbinsmid)
    Nprops = 4

    ## Make plot
    fig, axes = plt.subplots(Nbigbins, Nprops, figsize=(8*Nprops, 8*Nbigbins))

    for bigbin_ix in range(Nbigbins):
        for j in range(Nprops):
            prop_ix = j+1
            x = all_bins_mid[prop_ix]
            
            ax = axes[bigbin_ix, j]
            pdfs = allpdfs[j][:,bigbin_ix,:]
            
            y1,y2,y3 = medianscatter(pdfs, percentiles=[50-95/2.,50.,50+95/2.], axis=0)
            ax.plot(x, y2, 'k', drawstyle='steps-mid', label="All halos")
            ax.fill_between(x,y1,y3,color='k',facecolor='k',alpha=.3,step='mid')

            for k, (ufdtype,color,ls) in enumerate(zip(allufdtypes,ufdtypescolors,ufdlinestyles)):
                pdfs = allufdpdfs[j][k][:,bigbin_ix,:]
                y1,y2,y3 = medianscatter(pdfs, percentiles=[50-95/2.,50.,50+95/2.], axis=0)
                ax.plot(x, y2, color=color, ls=ls, drawstyle='steps-mid', label=ufdtypescolors)
                ax.fill_between(x, y1, y3, color=color, facecolor=color, alpha=.3, step='mid')

            ax.set_title("logM={:.1f}-{:.1f}".format(logmass_bigbins[bigbin_ix], logmass_bigbins[bigbin_ix+1]))
            ax.set_xlabel(prop_labels[prop_ix])
            ax.set_xlim(prop_xlims[prop_ix])
            ax.set_yscale('log')
        
    fig.savefig("haloprops_by_massbin_z{}.pdf".format(zin))
    plt.close(fig)
"""

def histogram_extant_over_all_halos_by_massbin(zin):
    hpaths = dm.get_hpaths(field=False, lx=14)

    Nhalos = len(hpaths)
    Nbigbins = len(logmass_bigbinsmid)
    Nprops = 4
    plug = AlexExtantDataPlugin()
    mbplug = MassAccrPlugin()
    zx_ = "z{}_".format(zin)

    allpdfs = [np.zeros((Nhalos, Nbigbins, len(all_bins_mid[j+1]))) for j in range(Nprops)]
    allufdpdfs = [[np.zeros((Nhalos, Nbigbins, len(all_bins_mid[j+1]))) \
                       for k in range(len(allufdtypes))] for j in range(Nprops)]
    allufdratios = [[np.zeros((Nhalos, Nbigbins, len(all_bins_mid[j+1]))) \
                         for k in range(len(allufdtypes))] for j in range(Nprops)]
    for i, hpath in enumerate(hpaths):
        start = time.time()
        # Load data
        hid = haloutils.get_parent_hid(hpath)
        zrobjs, rscatobjs, zrlogmass, zrlogvmax, zrspin, zrlogD, conc = load_one_halo_data(zin, hpath)
        zrTU = zrobjs["T/|U|"]
        with open("UFDSEARCH_Z0/{}_ufdids.pkl".format(haloutils.hidstr(hid)),"r") as fp:
            allufdids = pickle.load(fp)

        extant = plug.read(hpath)
        extant[zx_+"logmass"] = np.log10(extant[zx_+"mvir"]/h0)
        # Only use objects with masses at the relevant redshift
        extant = extant[pd.notnull(extant[zx_+"logmass"])]
        #extant[zx_+"logvmax"] = np.log10(extant[zx_+"vmax"])
        extant[zx_+"conc"] = extant[zx_+"rvir"]/extant[zx_+"rs"]
        
        mb = mbplug.read(hpath)
        z_r, snap_r = zin_to_zr_snapr(zin, verbose=False)
        iizin = mb['snap'] == snap_r
        assert np.sum(iizin) == 1, "{} has {}".format(hid,np.sum(iizin))
        halopos = extant[[zx_+"pos"+pos for pos in ['X','Y','Z']]].as_matrix()
        hostpos = np.array(mb[iizin][['x','y','z']]).view(np.float).reshape(-1,3)
        extant[zx_+"dist"] = np.sqrt(np.sum((halopos-hostpos)**2, axis=1)) * 1000.
        extant[zx_+"logD"] = np.log10(extant[zx_+"dist"])
        
        # Make histograms
        ### ERROR OOOOOPS: SWAPPED ORDER OF LOGD AND CONC HERE COMPARED TO LABELS
        props_to_hist = [zrTU, zrspin, zrlogD, conc]
        extant_props_to_hist = [extant[zx_+propname] for propname in ["T/|U|","spin","logD","conc"]]
        assert len(props_to_hist) == Nprops
        bigbin_indices = np.digitize(zrlogmass, bins=logmass_bigbins)
        extant_bigbin_indices = np.digitize(extant[zx_+"logmass"], bins=logmass_bigbins)
                                
        for bigbin_ix in range(Nbigbins):
            bigbin_ii = bigbin_indices == (bigbin_ix + 1)
            extant_bigbin_ii = extant_bigbin_indices == (bigbin_ix + 1)
            for j in range(Nprops):
                prop_ix = j+1 # skip mvir
                bins = all_bins[prop_ix]
                binsmid = all_bins_mid[prop_ix]
    
                ii = bigbin_ii
                x = props_to_hist[j][ii]
                x = x[np.isfinite(x)]
                h1, _ = np.histogram(x, bins=bins)
                y = h1/float(len(x))
                dx = bins[1]-bins[0]
                y /= dx
                allpdfs[j][i,bigbin_ix,:] = y
                
                for k, (ufdids, ufdtype) in enumerate(zip(allufdids, allufdtypes)):
                    ufdids = np.array(ufdids).astype(int)
                    # Filter objects in this mass bin
                    x = extant_props_to_hist[j][extant_bigbin_ii]
                    # Only grab extant UFD progenitors
                    x = x.ix[ufdids].as_matrix()
                    x = x[np.isfinite(x)]
                    h2, _ = np.histogram(x, bins=bins)
                    y = h2/float(len(x))
                    y /= dx
    
                    allufdpdfs[j][k][i,bigbin_ix,:] = y
                    allufdratios[j][k][i,bigbin_ix,:] = h2/h1.astype(float)
        print "{} took {:.1f}".format(hid, time.time()-start)
    
    with open("UFDSEARCH_Z0/hist_extant_by_mass_bin_z{}.pkl".format(zin),"w") as fp:
        pickle.dump([allpdfs, allufdpdfs, allufdratios], fp)

def plot_extant_massbin_histograms(zin, plot_type="count", log=False):
    """ note count is actually a PDF """
    assert plot_type in ["count","ratio"], plot_type
    
    with open("UFDSEARCH_Z0/hist_extant_by_mass_bin_z{}.pkl".format(zin),"r") as fp:
        allpdfs, allufdpdfs, allufdratios = pickle.load(fp)

    hpaths = dm.get_hpaths(field=False, lx=14)
    Nhalos = len(hpaths)
    Nbigbins = len(logmass_bigbinsmid)
    Nprops = 4

    ## Make plot
    fig, axes = plt.subplots(Nbigbins, Nprops, figsize=(8*Nprops, 8*Nbigbins))

    ### ERROR OOOOOPS: SWAPPED ORDER OF LOGD AND CONC HERE COMPARED TO ACTUAL DATA IN HISTOGRAMS
    for bigbin_ix in range(Nbigbins):
        for j in range(Nprops):
            prop_ix = j+1
            x = all_bins_mid[prop_ix]
            
            ax = axes[bigbin_ix, j]
            
            if plot_type == "count":
                pdfs = allpdfs[j][:,bigbin_ix,:]
                y1,y2,y3 = medianscatter(pdfs, percentiles=[50-95/2.,50.,50+95/2.], axis=0)
                ax.plot(x, y2, 'k', drawstyle='steps-mid', label="All halos")
                ax.fill_between(x,y1,y3,color='k',facecolor='k',alpha=.3,step='mid')

            for k, (ufdtype,color,ls) in enumerate(zip(allufdtypes,ufdtypescolors,ufdlinestyles)):
                if plot_type == "count":
                    pdfs = allufdpdfs[j][k][:,bigbin_ix,:]
                elif plot_type == "ratio":
                    pdfs = allufdratios[j][k][:,bigbin_ix,:]
                y1,y2,y3 = medianscatter(pdfs, percentiles=[50-95/2.,50.,50+95/2.], axis=0)
                ax.plot(x, y2, color=color, ls=ls, drawstyle='steps-mid', label=ufdtypescolors)
                ax.fill_between(x, y1, y3, color=color, facecolor=color, alpha=.3, step='mid')

            ax.set_title("logM={:.1f}-{:.1f}".format(logmass_bigbins[bigbin_ix], logmass_bigbins[bigbin_ix+1]))
            ax.set_xlabel(prop_labels[prop_ix])
            ax.set_xlim(prop_xlims[prop_ix])
            if plot_type == "count" and log:
                ax.set_yscale('log')
            if plot_type == "ratio":
                ax.set_yscale('log')
                ax.set_ylim(.001,1)
        
    if plot_type == "count":
        if log:
            fig.savefig("extant_haloprops_by_massbin_log_z{}.pdf".format(zin))
        else:
            fig.savefig("extant_haloprops_by_massbin_z{}.pdf".format(zin))
    elif plot_type == "ratio":
        fig.savefig("extant_halopropsratio_by_massbin_z{}.pdf".format(zin))
    return fig

if __name__=="__main__":
    for zin in [4,6,8,10,12]:
        #histogram_extant_over_all_halos_by_massbin(zin)
        fig = plot_extant_massbin_histograms(zin, plot_type="count", log=False)
        plt.close(fig)
        fig = plot_extant_massbin_histograms(zin, plot_type="count", log=True)
        plt.close(fig)
        #plot_extant_massbin_histograms(zin, plot_type="ratio")
