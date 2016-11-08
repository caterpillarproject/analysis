import numpy as np
import haloutils
import time,sys
from caterpillaranalysis import PluginBase

import readhalos.RSDataReader as RDR
import mergertrees.MTCatalogue as MTC
import readsnapshots.readsnapHDF5_greg as rsg


if __name__=="__main__":
    hpath = haloutils.get_hpath_lx(1387186, 13)
    z_r = 8
    snap_r = haloutils.get_snap_z(hpath, [z_r])
    
    
