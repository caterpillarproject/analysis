import numpy as np
import readhalos.RSDataReader as RDR
import asciitable
import glob
from optparse import OptionParser

import haloutils


if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("-w","--write-summary",dest="write_summary",
                      action="store_true",default=False,
                      help="flag to write a summary file in /bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index.txt")
    (options, args) = parser.parse_args()

    levellist = [11,12,13,14]
    nrvirlist = [4]

    hlist = haloutils.find_halo_paths(require_rockstar=True,
                            nrvirlist=nrvirlist,levellist=levellist)
    print "Number of halos with rockstar:",len(hlist)

    hindex = []
    
    for hpath in hlist:
        parenthid = haloutils.get_parent_hid(hpath)
        ictype,lx,nv = haloutils.get_zoom_params(hpath)
        lastsnap = haloutils.get_numsnaps(hpath) - 1
        hcat = RDR.RSDataReader(hpath+'/halos',lastsnap,version=6)
        zoomid = hcat.data.index[hcat['npart'].argmax()]
        zoommass = hcat.ix[zoomid]['mvir']
        print "Halo %7i Snap %3i has ID %6i: hM = %3.2e" % (parenthid, lastsnap, zoomid, zoommass)
        hindex.append([parenthid,lx,nv,zoomid])
    
    if options.write_summary:
        asciitable.write(hindex,"/bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index.txt",names=['parentid','LX','NV','zoomid'])
