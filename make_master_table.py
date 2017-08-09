import numpy as np
import time,sys
import pandas as pd
import cPickle as pickle
import matplotlib.pyplot as plt
from scipy import optimize

import haloutils
from caterpillaranalysis import MassAccrPlugin, PluginBase

sys.path.append("./greg_dwarfs")
import DwarfMethods as dm

from select_z8_objects import zin_to_zr_snapr
from classify_z8_objects import load_one_halo_data
from classify_z8_objects import allufdtypes
from trace_z0_ufds_to_zr import AlexExtantDataPlugin
    
def _solve_for_vmaxconc(val):
    """
    Eqn 9 of Dooley et al. 2014 (from Springel et al. 2008)
    c^3/(ln(1+c)-c/(1+c)) = 0.21639 * (Vmax/(H0 Rvmax))**2
    
    The RHS here is val so we numerically solve the equation.
    """
    def myfn(c):
        return c**3/(np.log(1.+c) - c/(1+c)) - val
    try:
        x0, r = optimize.brentq(myfn, 1.0, 100.0, full_output=True)
        if r.converged: return x0
        return np.nan
    except:
        return np.nan

solve_for_vmaxconc = np.vectorize(_solve_for_vmaxconc)

if __name__=="__main__":
    zin_list = [4,6,8,10,12]
    hpaths = dm.get_hpaths(field=False, lx=14)
    extantplug = AlexExtantDataPlugin()
    mbplug = MassAccrPlugin()
    
    #head_columns = ["hid","mtkey","zin"]
    #head_columns = ["mtkey"]
    #prop_columns = ["mvir","vmax","conc","T/|U|","spin","logD"]
    #surv_columns = ["surv","maxm","h14m","h14r","h14i"]

    for hpath in hpaths:
        start = time.time()
        hid = haloutils.get_parent_hid(hpath)
        print "Starting {}".format(hid)
        
        mb = mbplug.read(hpath)
        
        # Load UFD and extant/surv IDs
        with open("UFDSEARCH_Z0/{}_ufdids.pkl".format(haloutils.hidstr(hid)),"r") as fp:
            allufdids = pickle.load(fp)
        extant = extantplug.read(hpath)
        
        for zin in zin_list:
            z,snap = zin_to_zr_snapr(zin)
            start2 = time.time()
            #zrobjs = pd.DataFrame(np.load("UFDSEARCH_Z0/{}_z{}halos.npy".format(haloutils.hidstr(hid),zin)))
            output = load_one_halo_data(zin, hpath, use_vmaxconc=True)
            zrobjs = pd.DataFrame(output[0])
            vmaxconc = output[-1]
            nbad = np.sum(np.isnan(vmaxconc))
            nall = len(vmaxconc)
            print "   z={} vmaxconc missing {}/{} objects ({:.3f})".format(zin, nbad, nall, float(nbad)/nall)
            start3 = time.time()
            actual_vmaxconc = solve_for_vmaxconc(vmaxconc)
            print "   Took {:.1f} to solve for vmaxconc".format(time.time()-start3)
            zrobjs["vmaxconc"] = actual_vmaxconc

            zx_ = "z{}_".format(zin)
    
            # Store some extra fields
            # logD
            ii = (mb['snap'] == snap)
            assert np.sum(ii) == 1, "{} has {}".format(hid,np.sum(ii))
            halopos = np.array(zrobjs[['posX','posY','posZ']])#.view(np.float).reshape(-1,3)
            hostpos = np.array(mb[ii][['x','y','z']]).view(np.float).reshape(-1,3)
            zrobjs["logD"] = np.log10(np.sqrt(np.sum((halopos - hostpos)**2,axis=1))) + 3
            # conc
            zrobjs["conc"] = zrobjs["rvir"]/zrobjs["rs"]
    
            ## Tag objects that are main branch progenitors of extant/surv UFDs
            survids_at_zin = set(np.array(extant[zx_+"origid"]))
            zrobjs["surv"] = map(lambda x: x in survids_at_zin, zrobjs["origid"])
            #allufdids_at_zin = []
            for ufdtype, ufdids in zip(allufdtypes, allufdids):
                ufdids_at_zin = set(np.array(extant.ix[ufdids][zx_+"origid"]))
                #allufdids_at_zin.append(ufdids_at_zin)
                zrobjs[ufdtype] = map(lambda x: x in ufdids_at_zin, zrobjs["origid"])
            
            ## Add some quantities of the extant halos
            mtkeys = zrobjs["mtkey"].as_matrix()
            zrobjs["extant_maxmass"] = np.array(extant.ix[mtkeys]["maxmass_mvir"])
            zrobjs["halfmaxmass_scale"] = np.array(extant.ix[mtkeys]["halfmaxmass_scale"])
            zrobjs["infall_scale"] = np.array(extant.ix[mtkeys]["infall_scale"])
            zrobjs["infall_mvir"] = np.array(extant.ix[mtkeys]["infall_mvir"])
            zrobjs["firstinfall_scale"] = np.array(extant.ix[mtkeys]["firstinfall_scale"])
            zrobjs["firstinfall_mvir"] = np.array(extant.ix[mtkeys]["firstinfall_mvir"])
            #"halfmaxmass_scale"
            #"infall_scale"
            #"infall_mvir"
            #"firstinfall_scale"
            #"firstinfall_mvir"

            #u'scale', u'id', u'desc_id', u'num_prog', u'pid', u'upid', u'phantom',
            #u'sam_mvir', u'mvir', u'rvir', u'rs', u'vrms', u'mmp',
            #u'scale_of_last_MM', u'vmax', u'posX', u'posY', u'posZ', u'pecVX',
            #u'pecVY', u'pecVZ', u'Jx', u'Jy', u'Jz', u'spin', u'bfid', u'dfid',
            #u'origid', u'lastprog_dfid', u'm200c_all', u'm200b', u'xoff', u'voff',
            #u'spin_bullock', u'b_to_a(500c)', u'c_to_a(500c)', u'A[x](500c)',
            #u'A[y](500c)', u'A[z](500c)', u'T/|U|', u'snap', u'mtkey'
            output = zrobjs[["mtkey",
                             "mvir", "vmax", "vrms", "conc", "T/|U|", "spin", "logD",
                             "rvir","rs","vmaxconc","spin_bullock",
                             "scale","scale_of_last_MM",
                             "posX","posY","posZ","pecVX","pecVY","pecVZ","Jx","Jy","Jz",
                             "extant_maxmass","halfmaxmass_scale",
                             "infall_scale","infall_mvir",
                             "firstinfall_scale","firstinfall_mvir",
                             "surv", "maxm", "h14m", "h14r", "h14i",
                             "pid", "upid","mmp",
                             "origid","phantom"]]
    
            ## Save output as npy structured array
            np.save("UFDSEARCH_Z0/{}_z{}haloprops.npy".format(haloutils.hidstr(hid),zin),
                    output.to_records(index=False))
            print "    z={} {:.1f}".format(zin, time.time()-start2)
        print "Took {:.1f}".format(time.time()-start)
