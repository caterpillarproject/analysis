import numpy as np
import haloutils
import FindMiniHalos
import time,sys
from numpy.lib.recfunctions import append_fields
from caterpillaranalysis import PluginBase

class MinihaloFinderPlugin(PluginBase):
    def __init__(self,verbose=False,Tvir=2000):
        super(MinihaloFinderPlugin,self).__init__()
        self.filename = 'minihalo_array.npy'
        if Tvir != 2000:
            self.filename = 'minihalo_array_{}.npy'.format(int(Tvir))
        self.verbose = verbose
        self.Tvir = Tvir

    def is_above_threshold(self,scale,mass):
        z = 1./scale-1.0
        mcrit = FindMiniHalos.mcrit(self.Tvir,z)
        return mass > mcrit
    def is_above_threshold_and_not_phantom(self,scale,mass,phantom):
        threshold = self.is_above_threshold(scale,mass)
        not_phantom = phantom == 0
        return threshold * not_phantom
    
    def search_one_tree(self,mt,h0=0.6711):
        start = time.time()
        # This step is really slow but not bad enough to warrant fixing yet
        mt.data = append_fields(mt.data,['base_rsid','row'],data=[np.zeros(len(mt.data),dtype=int)+mt.rockstar_id,np.arange(len(mt.data),dtype=int)])
        #if self.verbose: print "Time to append: {:.1f}".format(time.time()-start)
    
        threshold = self.is_above_threshold_and_not_phantom(mt['scale'],mt['mvir']/h0,mt['phantom'])
        minihalo_candidates = mt.data[threshold]
        is_minihalo = np.zeros(len(minihalo_candidates),dtype=bool)
        for i,mh in enumerate(minihalo_candidates):
            dfid_base = mh['dfid']; dfid_last = mh['lastprog_dfid']
            mask = np.logical_and(mt.data['dfid']>=dfid_base, mt.data['dfid']<=dfid_last)
            if np.sum(threshold[mask]) == 1:
                is_minihalo[i] = True
        minihalos = minihalo_candidates[is_minihalo]
        if self.verbose: print "Time to find {} minihalos in tree {}: {:.1f}".format(len(minihalos),mt.rockstar_id,time.time()-start)
        return minihalos

    def _analyze(self,hpath):
        mtc = haloutils.load_zoom_mtc(hpath,indexbyrsid=True)
        all_minihalos = []
        start = time.time()
        for base_rsid,mt in mtc.Trees.iteritems():
            all_minihalos.append(self.search_one_tree(mt))
        print "Total time: {:.1f}".format(time.time()-start)
        all_minihalos = np.concatenate(all_minihalos)
        #assert not all_minihalos.masked,"If masked there is probably some data corruption!"
        all_minihalos = np.array(all_minihalos)
        np.save(self.get_outfname(hpath),all_minihalos)
    
    def _read(self,hpath):
        return np.load(self.get_outfname(hpath))

if __name__=="__main__":
    assert len(sys.argv)==3
    hid = int(sys.argv[1])
    lx = int(sys.argv[2])
    assert lx==14
    plug = MinihaloFinderPlugin(verbose=True)
    hpath = haloutils.get_hpath_lx(hid,lx)
    MHs = plug.read(hpath,recalc=True)
    
