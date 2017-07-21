import numpy as np
import MTanalysis_field as mtaf
FieldData = mtaf.FieldHaloSubstructureFirstPass()
import DwarfMethods as dm
hpaths = dm.get_hpaths(True)

n_field =0

nsub_field=0
nsub_host=0

for hpath in hpaths:
    dataE = FieldData.read(hpath)
    n_field += len(np.unique(dataE['field_rsid']))

    dataE_host = dm.get_extant_data(hpath, False)
    dataE_field = dm.get_extant_data(hpath, True)
    nsub_field += len(dataE_field)
    nsub_host += len(dataE_host)

    
print n_field
print nsub_field
print nsub_host
