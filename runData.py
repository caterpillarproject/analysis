import haloutils
from caterpillaranalysis import *
from caterpillarplot import *
from FindMiniHalos import *



"""
lx = 14
halo_paths = haloutils.get_paper_paths_lx(lx)
Destroyed = DestroyedDataFirstPass()
hpath = halo_paths[23]
#Destroyed.analyze(hpath,recalc=True)
AD = AllDestroyedData()
TM = TagMass()
AD.analyze(hpath,recalc=True)
TM.analyze(hpath,recalc=True)
dataD = Destroyed.read(hpath)
hid = haloutils.hpath_name(hpath)
"""

#hpath = '/bigbang/data/AnnaGroup/caterpillar/halos/H796175/H796175_EB_Z127_P7_LN7_LX13_O4_NV4'

#hpath = '/bigbang/data/AnnaGroup/caterpillar/halos/H1599988/H1599988_EX_Z127_P7_LN7_LX14_O4_NV4'


#Extant = ExtantDataFirstPass()
#Destroyed = DestroyedDataFirstPass()
#AE = AllExtantData()
#AD = AllDestroyedData()
#TM = TagMass()

#Extant.analyze(hpath,recalc=True)
#Destroyed.analyze(hpath,recalc=True)
#AE.analyze(hpath,recalc=True)
#AD.analyze(hpath,recalc=True)
#TM.analyze(hpath,recalc=True)

def runSingle(hid,lx):
    hpath = haloutils.get_hpath_lx(hid,lx)
    print hpath, "Starting"
    Extant = ExtantDataFirstPass()
    Destroyed = DestroyedDataFirstPass()
    AE = AllExtantData()
    AD = AllDestroyedData()
    TM = TagMass()
    #Extant.analyze(hpath,recalc=True)
    AE.analyze(hpath,recalc=True)
    Destroyed.analyze(hpath,recalc=True)
    AD.analyze(hpath,recalc=True)
    TM.analyze(hpath,recalc=True)
    print hpath, 'FINISHED'


#runSingle("H388476",14)
#runSingle("H388476",13)

def runCode(lx,start,end):
    Extant = ExtantDataFirstPass()
    Destroyed = DestroyedDataFirstPass()
    AE = AllExtantData()
    AD = AllDestroyedData()
    TM = TagMass()

    halo_paths = haloutils.get_paper_paths_lx(lx)
    for hpath in halo_paths[start:end]:
        if not hpath:
            print hpath, 'NOT FOUND'
            continue
        hid = haloutils.hpath_name(hpath)
        print hid, 'STARTING'
        #Extant.analyze(hpath,recalc=True)
        #AE.analyze(hpath,recalc=True)
        Destroyed.analyze(hpath,recalc=True)
        AD.analyze(hpath,recalc=True)
        TM.analyze(hpath,recalc=True)
        print hpath, 'FINISHED'

#    halo_paths = haloutils.find_halo_paths(levellist=[lx],require_mergertree=True,require_subfind=False,verbose=False) # length is 18 April 7th

if __name__ == "__main__":
    import sys
    lx = int(sys.argv[1])
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    runCode(lx,start,end)
 

# H1422331 starting

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

