import numpy as np
import pylab as plt

from minihalo_finder import MinihaloFinderPlugin, LWMinihaloFinderPlugin
import haloutils

def count_minihalos(hpath):
    rscat = haloutils.load_rscat(hpath,haloutils.get_lastsnap(hpath))
    zoomid= haloutils.load_zoomid(hpath)
    host = rscat.ix[zoomid]
    hpos = np.array(host[['posX','posY','posZ']])
    subs = rscat.get_all_subhalos_within_halo(zoomid)
    apos = rscat[['posX','posY','posZ']]
    subids = set(subs['id'].astype(int))

    #plug = MinihaloFinderPlugin()
    plug = LWMinihaloFinderPlugin()
    MHs = plug.read(hpath,autocalc=False)
    if MHs == None:
        return None
    N_MH = len(MHs)
    norm_Mhost = rscat.ix[zoomid]['mgrav']/rscat.h0/1e12
    norm_N = N_MH/norm_Mhost

    base_rsid = MHs['base_rsid']
    def where_is_halo(rsid):
        if rsid==zoomid: return -1.0
        if rsid in subids: return -2.0
        try:
            pos = np.array(apos.ix[rsid])
            return np.sqrt(np.sum((pos-hpos)**2))
        except KeyError:
            return -3.0
    halo_location = np.array(map(where_is_halo,np.array(base_rsid)))
    return N_MH,norm_N,halo_location,MHs
    
if __name__=="__main__":
    all_data = {}
    for cid in range(1,25):
        hid = haloutils.cid2hid[cid]
        shid= haloutils.hidstr(hid)
        hpath = haloutils.get_hpath_lx(hid,14)
        output = count_minihalos(hpath)
        if output==None:
            print shid,"not complete"
            continue
        all_data[hid] = output
        N_MH,norm_N,halo_location,MHs = output
        n_host = np.sum(halo_location==-1.0)
        n_dwarf= np.sum(halo_location==-2.0)
        n_other= np.sum(halo_location > 0)
        n_bad  = np.sum(halo_location==-3.0)
        print "{:8} N={:5}, N/M={:5.0f}, {:5}/{:5}/{:5} ({:5} bad)".format(shid,N_MH,norm_N,n_host,n_dwarf,n_other,n_bad)
