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
    - properties at halfmaxmass (~formation time)
    - properties at infall (to MW host)
    - properties at firstinfall (to any host)
    """
    def __init__(self):
        super(AlexExtantDataPlugin, self).__init__()
        self.filename='AlexExtantData.npy'
        self.autofigname = 'NAN'

        self.float_colnames = ['scale','mvir','rvir','rs','vrms','vmax','posX','posY','posZ','spin','T/|U|']
        self.int_colnames = ['snap','origid','phantom']
        self.property_colnames = self.float_colnames + self.int_colnames

        self.zin_prefixes = ['z'+str(zin) for zin in all_zin]
        self.prefixes = self.zin_prefixes+["maxmass","maxvmax","halfmaxmass","infall","firstinfall"]

    def _analyze(self, hpath):
        if not haloutils.check_mergertree_exists(hpath,autoconvert=True):
            raise IOError("No Merger Tree")
    
        hid = haloutils.get_parent_hid(hpath)
        
        property_colnames = self.property_colnames
        float_colnames = self.float_colnames
        int_colnames = self.int_colnames
    
        zin_prefixes = self.zin_prefixes
        prefixes = self.prefixes
        
        #am = abundance_matching.GarrisonKimmel()
        #LMmin, LMmax = am.stellar_to_halo_mass([1000., 2.e5])
        
        start = time.time()
        mtc = haloutils.load_zoom_mtc(hpath, indexbyrsid=True)
        Ntrees = len(mtc.Trees)
        print "Time to load zoom MTC with {} trees: {:.2f}".format(Ntrees, time.time()-start)
        property_dtypes = [mtc.fmttype[prop_name] for prop_name in property_colnames]
        
        # Grab the host main branch
        zoomid = haloutils.load_zoomid(hpath)
        hostmt = mtc[zoomid]
        mmp_map = hostmt.get_mmp_map()
        hostmb = hostmt.getMainBranch(mmp_map=mmp_map)
        
        # Setup output data
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
            colnames = [prefix+"_"+prop_name for prop_name in property_colnames]
            def fill_empty():
                for col in colnames:
                    if col.split("_")[1] in float_colnames:
                        output[col][i] = np.nan
                    elif col.split("_")[1] in int_colnames:
                        output[col][i] = -1
                    else:
                        raise RuntimeError
            
            if snap is None:
                fill_empty()
                return
            ii = np.where(mb['snap'] == snap)[0]
            if len(ii) == 0:
                fill_empty()
            else:
                assert len(ii)==1, ii
                mbrow = mb[ii][0]
                for col,propcol in zip(colnames,property_colnames):
                    output[col][i] = mbrow[propcol]
            return
    
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
            
            # halfmaxmass = zhalf of Mpeak
            snap_halfmaxmass = get_Mpeak_zhalf_snap(mb)
            fill_data(i,"halfmaxmass",mb,snap_halfmaxmass)
            
            # infall (to main MW)
            snap_infall = get_infall_to_host_snap(mb, hostmb)
            fill_data(i,"infall",mb,snap_infall)
            
            # firstinfall (to main MW)
            snap_firstinfall = get_infall_to_any_snap(mb)
            fill_data(i,"firstinfall",mb,snap_firstinfall)
            
            if i % 1000 == 0:
                print "  {:5} {:.1f}".format(i,time.time()-start)
                start = time.time()
        np.save(self.get_outfname(hpath), output)
    def _read(self, hpath):
        output = np.load(self.get_filename(hpath))
        index = output["mtkey"]
        return pd.DataFrame(output,index=index)

def get_Mpeak_zhalf_snap(mb):
    """ zhalf = first time you go above half of Mpeak """
    assert mb[0]['snap'] > mb[-1]['snap'], (mb[0]['snap'], mb[-1]['snap'])
    ix = np.argmax(mb["mvir"])
    mhalf = mb[ix]["mvir"]/2.
    # largest index (earlist time) where you go above half of mpeak
    ixhalf = np.max(np.where(mb["mvir"] > mhalf)[0])
    snap = mb[ixhalf]['snap']
    return snap

def get_infall_to_host_snap(submb, hostmb):
    """
    Snap before first time it is subhalo in main MW host
    Based on getInfall from MTanalysis3
    """
    Nsub = len(submb)
    Nhost = len(hostmb)
    if Nhost < Nsub:
        still_sub = np.where(hostmb["id"] == submb["upid"][0:Nhost])[0]
    else:
        still_sub = np.where(hostmb["id"][0:Nsub] == submb["upid"])[0]
    # If never actually a sub
    if len(still_sub)==0:
        return None
    # If began life as a subhalo
    if still_sub[-1] == Nsub-1:
        return None
    # Otherwise grab the snap right before it falls into the host
    ix = still_sub[-1] + 1
    # Greg did a check for phantoms, but I will ignore that
    return submb[ix]['snap']
    
def get_infall_to_any_snap(mb):
    """
    Snap before first time it is subhalo of ANY halo
    """
    still_sub = np.where(mb["upid"] != -1)[0]
    if len(still_sub) == 0:
        return None
    ix = np.max(still_sub)
    if ix == len(mb)-1:
        return None
    return mb[ix+1]['snap']

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
    hpaths = dm.get_hpaths(field=False, lx=14)
    alldata = []
    plug = AlexExtantDataPlugin()
    start = time.time()
    for hpath in hpaths:
        #hpath = haloutils.get_hpath_lx(1631506, 14)
        data = plug.read(hpath)
        alldata.append(data)
    am = abundance_matching.GarrisonKimmel()
    LMmin1, LMmax1 = am.stellar_to_halo_mass([1000., 2.e5])
    LMmin2, LMmax2 = am.stellar_to_halo_mass([2.e5, 5e7])
    ii1 = np.logical_and(data["maxmass_mvir"] >= LMmin1, data["maxmass_mvir"] <= LMmax1)
    ii2 = np.logical_and(data["maxmass_mvir"] >= LMmin2, data["maxmass_mvir"] <= LMmax2)
    zinfall = 1./data["infall_scale"] - 1
    finite_infall = np.isfinite(zinfall)
    zfirstinfall = 1./data["firstinfall_scale"] - 1
    finite_firstinfall = np.isfinite(zfirstinfall)
    fig,axes = plt.subplots(1,2)
    bins = np.arange(0,15,.5)
    ax = axes[0]
    ax.hist(zinfall[ii1 & finite_infall], label='UFDs', bins=bins, normed=True, histtype='step', color='k')
    ax.hist(zinfall[ii2 & finite_infall], label='UFD < gal < SMC', bins=bins, normed=True, histtype='step', color='r')
    ax.set_xlabel('z_infall MW')
    print float(np.sum(zinfall < 1))/np.sum(finite_infall)
    print float(np.sum((zinfall < 1) & ii1))/np.sum(finite_infall & ii1)
    print float(np.sum((zinfall < 1) & ii2))/np.sum(finite_infall & ii2)
    ax.legend(loc='upper right',fontsize=12)
    ax = axes[1]
    ax.hist(zfirstinfall[ii1 & finite_firstinfall], label='UFDs', bins=bins, normed=True, histtype='step', color='k')
    ax.hist(zfirstinfall[ii2 & finite_firstinfall], label='UFD < gal < SMC', bins=bins,normed=True, histtype='step', color='r')
    ax.set_xlabel('z_infall any')
    print float(np.sum(zfirstinfall < 1))/np.sum(finite_firstinfall)
    print float(np.sum((zfirstinfall < 1) & ii1))/np.sum(finite_firstinfall & ii1)
    print float(np.sum((zfirstinfall < 1) & ii2))/np.sum(finite_firstinfall & ii2)
    fig.savefig("test.png", bbox_inches='tight')

def tmp():
    data = pd.concat(alldata, ignore_index=True)
    print time.time()-start
    #data = data[pd.notnull(data["z8_T/|U|"])]
    data = data[data["z8_mvir"] > 1e7]
    #cols = ["z{}_T/|U|".format(zin) for zin in all_zin]
    #cols = ["z{}_rs".format(zin) for zin in all_zin]
    #cols = ["z{}_rvir".format(zin) for zin in all_zin]
    cols1 = ["z{}_rvir".format(zin) for zin in all_zin]
    cols2 = ["z{}_rs".format(zin) for zin in all_zin]
    
    vals = np.zeros((len(data.index),len(all_zin)))
    for i,(index,row) in enumerate(data.iterrows()):
        #vals[i,:] = np.array(row[cols])
        vals[i,:] = np.array(row[cols1])/np.array(row[cols2])
    meds = np.nanmedian(vals,axis=0)
    p0,p1,p2,p3,p4 = np.nanpercentile(vals,[2.5,16,50,84,97.5],axis=0)

    fig, ax = plt.subplots()
    #for i,row in data.iterrows():
    #    ax.plot(all_zin, row[cols], lw=.5)
    #    #ax.plot(all_zin, np.array(row[cols1])/np.array(row[cols2]), lw=.5)
    ax.plot(all_zin, p2, lw=2, color='k')
    ax.plot(all_zin, p3, 'k:', lw=2)
    ax.plot(all_zin, p1, 'k:', lw=2)
    ax.plot(all_zin, p4, 'k--', lw=2)
    ax.plot(all_zin, p0, 'k--', lw=2)
    ax.set_xlabel('z')
    
    ax.set_ylim(0,30)
    ax.set_ylabel('conc')
    plt.savefig("extant_zevol_conc.pdf")
    
    #ax.set_ylabel('rs')
    #ax.set_ylim(0,10)
    #plt.savefig("extant_zevol_rs.pdf")

    #ax.set_ylabel('rvir')
    #ax.set_ylim(0,30)
    #plt.savefig("extant_zevol_rvir.pdf")

def tmp():
    #run_plugin()

    for zin in all_zin:
        print zin
        histogram_all_halos_extant(zin)
        plot_histograms_extant(zin, plot_type="ratio")
        plot_histograms_extant(zin, plot_type="count")
