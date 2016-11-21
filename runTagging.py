import haloutils
from caterpillaranalysis import *
from caterpillarplot import *
from MTanalysis import *


TagExtant = TagExtantPlugin()
TagDestroyed = TagDestroyedPlugin()
TagMass = TagMass()

lx = 11
haloidlist = get_haloidlist(1)
for hid in haloidlist[1:]:
    hpath = haloutils.get_hpath_lx(hid,lx)
    TagExtant.analyze(hpath)
    TagDestroyed.analyze(hpath)
    TagDestroyed.combinefiles(hpath)
    TagMass.analyze(hpath)

"""
lx = 12
haloidlist = get_haloidlist(1)
for hid in haloidlist:
    hpath = haloutils.get_hpath_lx(hid,lx)
    TagExtant.analyze(hpath)
    TagDestroyed.analyze(hpath)
    TagDestroyed.combinefiles(hpath)
    TagMass.analyze(hpath)
"""
