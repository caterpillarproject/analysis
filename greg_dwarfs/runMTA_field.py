import numpy as np
import MTanalysis_field as mtaf
import haloutils
import sys
import DwarfMethods as dm
lx=14
#hpaths = haloutils.get_all_halo_paths_lx(lx)
#hpaths = [haloutils.catnum_hpath(36,lx),haloutils.catnum_hpath(37,lx), haloutils.catnum_hpath(40,lx),haloutils.catnum_hpath(53,lx)]

hpaths = dm.get_hpaths(field=False)[31:]


for hpath in hpaths:
    hid = haloutils.get_parent_hid(hpath)
    print hid, 'hid that is running'
    sys.stdout.flush()
    
    #obj = mtaf.HostHaloM350()
    #data = obj.analyze(hpath,recalc=True)
    field = mtaf.FieldHaloSubstructureFirstPass()
    data = field.analyze(hpath,recalc=True)
    field_all = mtaf.AllExtantFieldData()
    data_all = field_all.analyze(hpath,recalc=True)
    print 'done with', hid

    #except:
     #   print 'halo did not run'
      #  sys.stdout.flush()
    
