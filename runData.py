import haloutils
from caterpillaranalysis import *
from caterpillarplot import *
from MTanalysis2 import *


def runCode(lx,start,end):
    Extant = ExtantDataFirstPass()
    Destroyed = DestroyedDataFirstPass()
    AE = AllExtantData()
    AD = AllDestroyedData()
    halo_paths = haloutils.find_halo_paths(levellist=[lx],require_mergertree=True,require_subfind=False,verbose=False) # length is 18 April 7th
    for hpath in halo_paths[start:end]:
        print hpath, 'STARTING'
        Extant.analyze(hpath,recalc=True)
        AE.analyze(hpath,recalc=True)
        Destroyed.analyze(hpath,recalc=True)
        AD.analyze(hpath,recalc=True)
        #Destroyed.combinefiles(hpath) # now part of analyze
        print hpath, 'FINISHED'

if __name__ == "__main__":
    import sys
    lx = int(sys.argv[1])
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    runCode(lx,start,end)
    

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

