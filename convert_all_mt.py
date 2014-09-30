import numpy as np
from haloutils import get_foldername,find_halo_paths
import mergertrees.MTCatalogue as MTC
import os

def convert_mt_in_dir(outpath,version=4):
    if os.path.exists(outpath+'/halos/trees/tree_0_0_0.dat'):
        if not os.path.exists(outpath+'/halos/trees/tree.bin'):
            print '---converting: '+get_foldername(outpath)
            MTC.convertmt(outpath+'/halos/trees',version=version)
        else:
            print '---already converted '+get_foldername(outpath)
    else:
        print '---MISSING! '+get_foldername(outpath)
    
def convert_all_mt():
    halopathlist = find_halo_paths(require_rockstar=True)
    for outpath in halopathlist:
        convert_mt_in_dir(outpath)

if __name__=="__main__":
    #convert_all_mt()
    path = "/bigbang/data/AnnaGroup/caterpillar/halos/H1327707/H1327707_BB_Z127_P7_LN7_LX14_O4_NV4/halos/trees"
    MTC.convertmt(path,version=4,filenamein=path+'/tree_1_1_1.dat')
