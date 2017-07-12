import numpy as np
import haloutils
import time,sys
import pandas as pd
import matplotlib.pyplot as plt
import cPickle as pickle

import MTanalysis3 as mta
import MTaddition as mtadd
from FindMiniHalos import mcrit
from caterpillaranalysis import PluginBase, MassAccrPlugin

sys.path.append("./greg_dwarfs")
import DwarfMethods as dm
import abundance_matching

from select_z8_objects import zin_to_zr_snapr, all_zin
from classify_z8_objects import medianscatter #, load_one_halo_data, histogram_one_halo
from classify_z8_objects import logMbins, logMbinsmid, allufdtypes, ufdtypescolors, ufdlinestyles
from classify_z8_objects import logVmaxbins, logVmaxbinsmid, concbins, concbinsmid
from classify_z8_objects import logDbins, logDbinsmid, spinbins, spinbinsmid, TUbins, TUbinsmid
from classify_z8_objects import all_bins, all_bins_mid, prop_labels, prop_xlims

class AlexExtantDataPlugin(PluginBase):
    """
    This essentially recalculates some of Greg's catalogs but customized for me.
    It is basically a subset of the main branch 
    Properties I store:
    - mvir, rvir, rs, vrms, vmax, x, y, z, spin, T/|U|, snap, origid, phantom
    For each halo in the MTC (at z=0), I calculate:
    - properties at zin = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14] from the merger tree
    - properties at maxmass
    - properties at maxvmax
    
    - properties at infall [NOT DONE, takes a long time]
    """
    def __init__(self):
        super(AlexExtantDataPlugin, self).__init__()
        self.filename='AlexExtantData.npy'
        self.autofigname = 'NAN'

        self.property_colnames = ['mvir','rvir','rs','vrms','vmax','posX','posY','posZ','spin','T/|U|','snap','origid','phantom']

        self.zin_prefixes = ['z'+str(zin) for zin in all_zin]
        self.prefixes = self.zin_prefixes+["maxmass","maxvmax"]

    def _analyze(self, hpath):
        if not haloutils.check_mergertree_exists(hpath,autoconvert=True):
            raise IOError("No Merger Tree")
    
        hid = haloutils.get_parent_hid(hpath)
        
        property_colnames = self.property_colnames
        float_colnames = ['mvir','rvir','rs','vrms','vmax','posX','posY','posZ','spin','T/|U|']
        int_colnames = ['snap','origid','phantom']
    
        zin_prefixes = self.zin_prefixes
        prefixes = self.prefixes
        
        #am = abundance_matching.GarrisonKimmel()
        #LMmin, LMmax = am.stellar_to_halo_mass([1000., 2.e5])
        
        start = time.time()
        mtc = haloutils.load_zoom_mtc(hpath, indexbyrsid=True)
        Ntrees = len(mtc.Trees)
        print "Time to load zoom MTC with {} trees: {:.2f}".format(Ntrees, time.time()-start)
        property_dtypes = [mtc.fmttype[prop_name] for prop_name in property_colnames]
        
        
        all_colnames = ["mtkey"]
        all_dtypes = [np.dtype(np.int)]
        for prefix in prefixes:
            for prop_name,prop_dtype in zip(property_colnames,property_dtypes):
                all_colnames.append(prefix+'_'+prop_name)
                all_dtypes.append(prop_dtype)
        Ncolumns = len(all_colnames)
    
        dtypes = np.dtype([(col,dtype) for col,dtype in zip(all_colnames,all_dtypes)])
        output = np.zeros(Ntrees, dtype=dtypes)
    
        def fill_data(i,prefix,mb,snap):
            """
            row i of output has columns according to prefix
            Fill in from mb occording to snap
            """
            ii = np.where(mb['snap'] == snap)[0]
            colnames = [prefix+"_"+prop_name for prop_name in property_colnames]
            if len(ii) == 0:
                # Fill empty
                for col in colnames:
                    if col.split("_")[1] in float_colnames:
                        output[col][i] = np.nan
                    elif col.split("_")[1] in int_colnames:
                        output[col][i] = -1
                    else:
                        raise RuntimeError
            else:
                assert len(ii)==1, ii
                mbrow = mb[ii][0]
                for col,propcol in zip(colnames,property_colnames):
                    output[col][i] = mbrow[propcol]
            pass
    
        start = time.time()
        for i,(mtkey, mt) in enumerate(mtc.Trees.iteritems()):
            output[i][0] = mtkey
            if len(mt.data) < 1e5:
                mb = mt.getMainBranch()
            else:
                mmp_map = mt.get_mmp_map()
                mb = mt.getMainBranch(mmp_map=mmp_map)
            
            # zin
            for zin,zin_prefix in zip(all_zin,zin_prefixes):
                zr, snap_zin = zin_to_zr_snapr(zin, verbose=False)
                fill_data(i,zin_prefix,mb,snap_zin)
            
            # maxmass snap
            ix = np.argmax(mb['mvir'])
            snap_maxmass = mb[ix]['snap']
            fill_data(i,"maxmass",mb,snap_maxmass)
            
            # maxvmax snap
            ix = np.argmax(mb['vmax'])
            snap_maxvmax = mb[ix]['snap']
            fill_data(i,"maxvmax",mb,snap_maxvmax)
            if i % 1000 == 0:
                print "  {:5} {:.1f}".format(i,time.time()-start)
                start = time.time()
        np.save(self.get_outfname(hpath), output)
    def _read(self, hpath):
        output = np.load(self.get_filename(hpath))
        index = output["mtkey"]
        return pd.DataFrame(output,index=index)

def run_plugin():
    hpaths = dm.get_hpaths(field=False, lx=14)
    try:
        whichone = int(sys.argv[1])
    except: # First argument doesn't work
        pass
    else:
        assert whichone in [0,1,2,3], whichone
        hpaths = hpaths[(whichone*8):((whichone+1)*8)]

    #plug = AlexExtantDataPlugin()
    plug = MassAccrPlugin()
    for hpath in hpaths:
        start = time.time()
        plug.analyze(hpath, recalc=True)
        #data = plug.read(hpath)
        print hpath, time.time()-start
        sys.stdout.flush()

def histogram_one_halo_extant(zin, hpath, full_output=False):
    h0 = .6711
    prefix = "z{}_".format(zin)

    ## UFD IDs
    hid = haloutils.get_parent_hid(hpath)
    with open("UFDSEARCH_Z0/{}_ufdids.pkl".format(haloutils.hidstr(hid)),"r") as fp:
        allufdids = pickle.load(fp)
    
    ## Extant data at z=1-14 and maxmass/maxvmax
    plug = AlexExtantDataPlugin()
    extantdf = plug.read(hpath)
    
    ## IMPORTANT: Here I make a cut to only consider halos with M > 10**6.5 to match my other catalog
    ## This cuts out a good number of halos
    extantdf = extantdf[(extantdf[prefix+"mvir"]/h0) > 10**6.5]

    allprophistssurv = []
    allprophistsufds = []
    #props_to_histogram = ["mvir","vmax","spin","conc","dist"]
    #log_or_not = [True, True, False, False, True]
    for i in range(len(all_bins_mid)):
        bins = all_bins[i]
        binsmid = all_bins_mid[i]
        
        # Ugly copy/paste but whatever
        if i==0:
            xx = np.log10(extantdf[prefix+"mvir"]/h0)
        elif i==1:
            #xx = np.log10(extantdf[prefix+"vmax"])
            xx = extantdf[prefix+"T/|U|"]
        elif i==2:
            xx = extantdf[prefix+"spin"]
        elif i==3:
            xx = extantdf[prefix+"rvir"]/extantdf[prefix+"rs"]
        elif i==4:
            plug2 = MassAccrPlugin()
            mb = plug2.read(hpath)
            ii = (mb['snap'] == np.unique(extantdf[prefix+"snap"])[-1]) # avoid -1
            assert np.sum(ii) == 1, "ERROR_A: {} has {} (want {} have {})".format(
                hid,np.sum(ii),np.unique(extantdf[prefix+"snap"]), mb['snap'])
            hostpos = np.array(mb[ii][['x','y','z']]).view(np.float).reshape(-1,3)
            halopos = extantdf[[prefix+pos for pos in ['posX','posY','posZ']]].as_matrix()
            xx = np.log10(np.sqrt(np.sum((halopos - hostpos)**2,axis=1))) + 3
        xx = np.array(xx)
        xx = xx[np.isfinite(xx)]
        h, _ = np.histogram(xx, bins=bins)
        allprophistssurv.append(h)

        thishistsufds = []
        for j in range(len(allufdtypes)):
            ufdids = allufdids[j]
            tdf = extantdf.ix[ufdids]
            # Ugly copy/paste but whatever
            if i==0:
                xx = np.log10(tdf[prefix+"mvir"]/h0)
            elif i==1:
                #xx = np.log10(tdf[prefix+"vmax"])
                xx = tdf[prefix+"T/|U|"]
            elif i==2:
                xx = tdf[prefix+"spin"]
            elif i==3:
                xx = tdf[prefix+"rvir"]/tdf[prefix+"rs"]
            elif i==4:
                plug2 = MassAccrPlugin()
                mb = plug2.read(hpath)
                ii = mb['snap'] == np.unique(extantdf[prefix+"snap"])[-1] # avoid -1
                assert np.sum(ii) == 1, "ERROR_B: {} has {} (want {} have {})".format(
                    hid,np.sum(ii),np.unique(extantdf[prefix+"snap"]), mb['snap'])
                hostpos = np.array(mb[ii][['x','y','z']]).view(np.float).reshape(-1,3)
                halopos = tdf[[prefix+pos for pos in ['posX','posY','posZ']]].as_matrix()
                xx = np.log10(np.sqrt(np.sum((halopos - hostpos)**2,axis=1))) + 3
            xx = np.array(xx)
            xx = xx[np.isfinite(xx)]
            h, _ = np.histogram(xx, bins=bins)
            thishistsufds.append(h)
        thishistsufds = thishistsufds
        allprophistsufds.append(thishistsufds)
    return allprophistssurv, allprophistsufds
def histogram_all_halos_extant(zin):
    num_hists = len(all_bins_mid)
    all_data = [[] for i in range(num_hists)]
    all_data2= [[] for i in range(num_hists)]

    hpaths = dm.get_hpaths(field=False, lx=14)
    for hpath in hpaths:
        start = time.time()
        try:
            h_out,hufd_out = histogram_one_halo_extant(zin, hpath)
            for x1,x2,h1,h2s in zip(all_data, all_data2, h_out, hufd_out):
                x1.append(h1)  # 32 x Nbin
                x2.append(h2s) # 32 x 4 x Nbin
        except Exception as e:
            print e
        print "{}: {:.2f}".format(haloutils.hidstr(haloutils.get_parent_hid(hpath)), time.time()-start)
            
    to_save = []
    for m,h1,h2s in zip(all_bins_mid,all_data,all_data2):
        to_save.append([m,np.array(h1),np.array(h2s)])
    with open("UFDSEARCH_Z0/haloprops_extant_at_z{}.pkl".format(zin),"w") as fp:
        pickle.dump(to_save, fp)
def plot_histograms_extant(zin, plot_type="count"):
    assert plot_type in ["count","ratio"], plot_type

    with open("UFDSEARCH_Z0/haloprops_extant_at_z{}.pkl".format(zin),"r") as fp:
        all_data = pickle.load(fp)
    fig, axes = plt.subplots(3,2,figsize=(8*2,8*3))
    labels = prop_labels
    xlims = prop_xlims
    for data,ax,label,xlim in zip(all_data,axes.flat[0:5],labels,xlims):
        x, h1, h2s = data
        y1,y2,y3 = medianscatter(h1, percentiles=[50-95/2.,50.,50+95/2.], axis=0)
        if plot_type == "count":
            ax.plot(x, y2, 'k', drawstyle='steps-mid', label="All surviving halos")
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
    ## Make legend in 6th axis
    ax = axes[0,0]
    handles,labels = ax.get_legend_handles_labels()
    lax = axes[2,1]
    lax.legend(handles,labels,loc='center')

    if plot_type == "count":
        fig.savefig("extant_haloprops_at_z{}.pdf".format(zin))
    elif plot_type == "ratio":
        fig.savefig("extant_halopropsratio_at_z{}.pdf".format(zin))
    plt.close(fig)


if __name__=="__main__":
    hpath = haloutils.get_hpath_lx(1631506, 14)
    plug = AlexExtantDataPlugin()
    data = plug.read(hpath)
    data = data[pd.notnull(data["z8_T/|U|"])]
    data = data[data["z8_mvir"] > 1e7]
    cols = ["z{}_T/|U|".format(zin) for zin in all_zin]
    
    fig, ax = plt.subplots()
    for i,row in data.iterrows():
        ax.plot(all_zin, row[cols], lw=.5)
    plt.savefig("test.pdf")

def tmp():
    #run_plugin()

    for zin in all_zin:
        print zin
        histogram_all_halos_extant(zin)
        plot_histograms_extant(zin, plot_type="ratio")
        plot_histograms_extant(zin, plot_type="count")
