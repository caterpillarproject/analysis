import numpy as np
import MTaddition_field as mtaf
import haloutils
import sys
import DwarfMethods as dm
lx=14
#hpaths = haloutils.get_all_halo_paths_lx(lx)
#hpaths = [haloutils.catnum_hpath(36,lx),haloutils.catnum_hpath(37,lx), haloutils.catnum_hpath(40,lx),haloutils.catnum_hpath(53,lx)]

hpaths = dm.get_hpaths(field=True)  


for hpath in hpaths[17:]:
    hid = haloutils.get_parent_hid(hpath)
    print hid, 'hid that is running'
    sys.stdout.flush()
    
    field = mtaf.ExtantDataReionizationField()
    try:
        data = field.analyze(hpath,recalc=True)
        print 'done with extant first pass', hid, 'Cat-', haloutils.hid_catnum(hid)
        sys.stdout.flush()

    except:
        print 'halo did not run', hid
        print ' '
        sys.stdout.flush()
    
