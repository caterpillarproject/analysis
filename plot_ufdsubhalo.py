from plot_ufd_surv_figs import *
from classify_z8_objects import load_one_halo_data

def plot_sub_vs_not(zin):
    hids = load_all_hids()
    hpaths = dm.get_hpaths(field=False, lx=14)
    #allzrobjs = []
    #allzrlogD = []
    #for hpath,hid in zip(hpaths,hids):
    #    out = load_one_halo_data(zin, hpath)
    #    zrobjs = out[0]
    #    zrlogD = out[5]
    #    #zrobjs = np.load("UFDSEARCH_Z0/{}_z{}halos.npy".format(haloutils.hidstr(hid),zin))
    #    allzrobjs.append(zrobjs)
    #    allzrlogD.append(zrlogD)
    #allzrobjs = pd.DataFrame(np.concatenate(allzrobjs))
    #allzrlogD = np.concatenate(allzrlogD)
    #allzrobjs["logD"] = allzrlogD
    #allzrobjs["conc"] = allzrobjs["rvir"]/allzrobjs["rs"]
    #allzrobjs["logmvir"] = np.log10(allzrobjs["mvir"]/h0)
    allzrobjs = load_all_data(zin_to_use=[zin], use_phantoms=global_use_phantoms)
    
    ii_subs = allzrobjs["pid"] != -1
    fig, axes = plt.subplots(5,1, figsize=(8,8*5))
    for ax, bins, bins_mid, xlabel, xlim, propcol, in \
            zip(axes.flatten(), all_bins, all_bins_mid, 
                prop_labels, prop_xlims, prop_cols):
        dx = bins[1]-bins[0]

        for ii,label in zip([np.ones(len(allzrobjs),dtype=bool), ii_subs, ~ii_subs],["All", "Subs","Not Subs"]):
            vals = allzrobjs[ii][propcol]
            vals2 = vals[np.isfinite(vals)]
            h, _ = np.histogram(vals2, bins=bins)
            y = h.astype(float)/dx/len(vals)
            ax.plot(bins_mid, y, label=label)
        ax.set_xlabel(xlabel)
        ax.set_xlim(xlim)
    axes[0].legend()
    fig.savefig("Sub_vs_not_z{}.pdf".format(zin))

if __name__=="__main__":
    for zin in [4,6,8,10,12]:
        plot_sub_vs_not(zin)
