import DwarfMethods as dm
import numpy as np
import haloutils as htils


hpaths = dm.get_hpaths(field=False,lx=14)
hpath = hpaths[0]
data = dm.get_extant_data(hpath)

data.dtypes

snap_z0 = htils.get_numsnaps(hpath)-1
cat=htils.load_rscat(hpath,snap_z0,rmaxcut=False)
hostID = htils.load_zoomid(hpath)
subhalos = cat.get_subhalos_within_halo(hostID)


common = np.in1d(np.array(data['rsid']), np.array(subhalos['id']))
print len(data)
print np.sum(common)



# I need the z=0 pid to tell if it is a subhalo or not
# in extant data, I should add upid.
# do I then need to run the 2nd pass? AllExtantData?
# Can I take the AllExtantData and append it to the first pass?
# Do I need to re-run the MTadd?
# For the field halos, do they need 


data[common]   # these are all the halos


import DwarfMethods as dm
hpath = dm.get_hpaths(field=False,lx=14)[0]
%run MTanalysis3.py
AE = AllExtantData()
alldata = AE.read(hpath)
ed = extract_dataE()
data2 = ed.read(hpath)
