import numpy as np
from findhalos.haloutils import get_numsnaps,get_foldername,find_halo_paths
import mergertrees.MTCatalogue as MTC
import os

def convert_mt_in_dir(outpath,version=2):
    numsnaps = get_numsnaps(outpath); snapstr = str(numsnaps-1)
    if os.path.exists(outpath+'/halos/halos_'+snapstr+'/halos_'+snapstr+'.0.bin'):
        if os.path.exists(outpath+'/halos/trees/tree_0_0_0.dat'):
            if not os.path.exists(outpath+'/halos/trees/tree.bin'):
                print '---converting: '+get_foldername(outpath)
                MTC.convertmt(outpath+'/halos/trees',version=version)
            else:
                print '---already converted '+get_foldername(outpath)
        else:
            print '---MISSING! '+get_foldername(outpath)
    
def convert_all_mt():
    halopathlist = find_halo_paths()
    for outpath in halopathlist:
        convert_mt_in_dir(outpath)

if __name__=="__main__":
    convert_all_mt()
