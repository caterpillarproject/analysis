import numpy as np
import time,sys
import pandas as pd
import cPickle as pickle
import matplotlib.pyplot as plt

import haloutils
from caterpillaranalysis import MassAccrPlugin, PluginBase

sys.path.append("./greg_dwarfs")
import DwarfMethods as dm

from select_z8_objects import zin_to_zr_snapr
from classify_z8_objects import load_one_halo_data
from classify_z8_objects import allufdtypes
from trace_z0_ufds_to_zr import AlexExtantDataPlugin
    
if __name__=="__main__":
    zin_list = [4,6,8,10,12]
    hpaths = dm.get_hpaths(field=False, lx=14)
    extantplug = AlexExtantDataPlugin()
    mbplug = MassAccrPlugin()
    
    #head_columns = ["hid","mtkey","zin"]
    #head_columns = ["mtkey"]
    #prop_columns = ["mvir","vmax","conc","T/|U|","spin","logD"]
     #surv_columns = ["surv","maxm","h14m","h14r","h14i"]

    for hpath in hpaths:
        start = time.time()
        hid = haloutils.get_parent_hid(hpath)
        print "Starting {}".format(hid)
        
        mb = mbplug.read(hpath)
        
        # Load UFD and extant/surv IDs
        with open("UFDSEARCH_Z0/{}_ufdids.pkl".format(haloutils.hidstr(hid)),"r") as fp:
            allufdids = pickle.load(fp)
        extant = extantplug.read(hpath)
        
        for zin in zin_list:
            start2 = time.time()
            zrobjs = pd.DataFrame(np.load("UFDSEARCH_Z0/{}_z{}halos.npy".format(haloutils.hidstr(hid),zin)))
            zx_ = "z{}_".format(zin)
    
            # Store some extra fields
            # logD
            z,snap = zin_to_zr_snapr(zin)
            ii = (mb['snap'] == snap)
            assert np.sum(ii) == 1, "{} has {}".format(hid,np.sum(ii))
            halopos = np.array(zrobjs[['posX','posY','posZ']])#.view(np.float).reshape(-1,3)
            hostpos = np.array(mb[ii][['x','y','z']]).view(np.float).reshape(-1,3)
            zrobjs["logD"] = np.log10(np.sqrt(np.sum((halopos - hostpos)**2,axis=1))) + 3
            # conc
            zrobjs["conc"] = zrobjs["rvir"]/zrobjs["rs"]
    
            ## Tag objects that are main branch progenitors of extant/surv UFDs
            survids_at_zin = set(np.array(extant[zx_+"origid"]))
            zrobjs["surv"] = map(lambda x: x in survids_at_zin, zrobjs["origid"])
            #allufdids_at_zin = []
            for ufdtype, ufdids in zip(allufdtypes, allufdids):
                ufdids_at_zin = set(np.array(extant.ix[ufdids][zx_+"origid"]))
                #allufdids_at_zin.append(ufdids_at_zin)
                zrobjs[ufdtype] = map(lambda x: x in ufdids_at_zin, zrobjs["origid"])
            
            output = zrobjs[["mtkey",
                             "mvir", "vmax", "vrms", "conc", "T/|U|", "spin", "logD",
                             "surv", "maxm", "h14m", "h14r", "h14i",
                             "pid", "upid",
                             "origid","phantom"]]
    
            ## Save output as npy structured array
            np.save("UFDSEARCH_Z0/{}_z{}haloprops.npy".format(haloutils.hidstr(hid),zin),
                    output.to_records(index=False))
            print "    z={} {:.1f}".format(zin, time.time()-start2)
        print "Took {:.1f}".format(time.time()-start)
