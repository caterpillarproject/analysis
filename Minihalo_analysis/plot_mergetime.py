import haloutils as htils
import sys,os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

hpath = htils.get_hpath_lx(1387186,14)
minihalos = np.load(hpath+"/analysis/minihalo_array.npy")
form_snap = np.load(hpath+'/analysis/form_snap.npy')
merge_snap = np.load(hpath+'/analysis/merge_snap.npy')
#merge_snap=np.array(merge_snap)
#form_snap = np.array(form_snap)

mask = merge_snap != 320
diff = merge_snap[mask]-form_snap[mask]
plt.hist(diff,bins=np.arange(310))
plt.xlabel('Num Snaps')
plt.ylabel('Num halos')
plt.savefig('merge_time')
