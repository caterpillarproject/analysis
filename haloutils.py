import os,sys,platform
from collections import defaultdict
# Allow plot creation on antares
if 'compute-0-' in platform.node():
    import matplotlib
    matplotlib.use('Agg')

import numpy as np
import asciitable
#import pickle
import cPickle as pickle
import pandas as pd
import warnings
import glob
from multiprocessing import Pool
import time,subprocess,itertools
#import seaborn as sns

import readsnapshots.readsnapHDF5_greg as rsg
import readhalos.RSDataReader as RDR
import readhalos.readsubf as RSF
import mergertrees.MTCatalogue as MTC
import brendanlib.conversions as bconversions
import pylab as plt

def determinebasepath(node):
    if node == "csr-dyn-150.mit.edu":
        basepath = '/Users/griffen/Desktop/'
    elif node == "Brendans-MacBook-Pro.local":
        basepath = '/Users/griffen/Desktop/'
    elif node == "spacebase":
        basepath = '/bigbang/data/AnnaGroup/'
    elif node == "bigbang.mit.edu":
        basepath = '/bigbang/data/AnnaGroup/'
    elif node == "blender.mit.edu":
        basepath = '/bigbang/data/AnnaGroup/'
    elif node == "antares":
        basepath = '/bigbang/data/AnnaGroup/'
    elif 'compute-0-' in node:
        basepath = '/bigbang/data/AnnaGroup/'
    elif 'Griffs-Macbook.local' in node:
        basepath = '-'
    else:
        raise ValueError(node+" is not a valid node")
        
    return basepath

cdict = {'red'  :  ((0., 0., 0.),     (0.3,0,0),     (0.6, 0.8, 0.8), (1., 1., 1.)),
'green':  ((0., 0., 0.),     (0.3,0.3,0.3), (0.6, 0.4, 0.4), (1., 1.0, 1.0)),
'blue' :  ((0., 0.05, 0.05), (0.3,0.5,0.5), (0.6, 0.6, 0.6), (1.0, 1.0, 1.0))}

cmap = plt.matplotlib.colors.LinearSegmentedColormap('dmdens_cmap', cdict, 1024)
plt.cm.register_cmap(name='caterpillar', cmap=cmap)


global_basepath = os.path.normpath(determinebasepath(platform.node()))
global_halobase = global_basepath+'/caterpillar/halos'
global_prntbase = global_basepath+'/caterpillar/parent/gL100X10'


cid2hid = {1: 1631506,2:  264569,3: 1725139,
           4:  447649,5:    5320,6:  581141,
           7:   94687,8: 1130025,9: 1387186,
           10: 581180,11:1725372,12:1354437,
           13:1725272,14:1195448,15:1292085,
           16: 796175,17: 388476,18:1079897,
           19:  94638,20:  95289,21:1232164,
           22:1422331,23: 196589,24:1268839,
           25:1599988,26:1195075,27:1631582,
           28: 831159,29:1422429,30:1848643,
           31:  65777,32:1269057,33:1232423,
           34:  94629,35:1232313,36: 196078,
           37:1599902,38: 196078,39:1232127,
           40: 795802,41:1161558,42: 616647,
           43: 795912,44:1080116,45:1354579,
           46:  94105,47:  94562,48: 581380,
           49: 650016,50: 795405,51: 795648,
           52: 831279,53:1104787,54:1105005,
           55:1161589,56:1194734,57:1194785,
           58:1354730,59:1422249,60:1476003,
           61:1507285,62:1697899,63:1725470,
           64:1764462,65:1940047}
#cid2hid = {1:1631506,
#2:264569,
#3:1725139,
#4:447649,
#5:5320,
#6:581141,
#7:94687,
#8:1130025,
#9:1387186,
#10:581180,
#11:1725372,
#12:1354437,
#13:1725272,
#14:1195448,
#25:1599988,
#16:796175,
#17:388476,
#18:1079897,
#19:94638,
#20:95289,
#21:1232164,
#22:1422331,
#23:196589,
#24:1268839,
#15:1292085,
#26:1195075,
#27:1631582,
#28:831159,
#29:1422429,
#30:1848643,
#31:65777,
#32:1269057,
#33:1232423,
#34:94629,
#35:1232313,
#36:196078,
#37:1599902,
#38:1161351,
#39:1232127,
#40:795802,
#41:1161558,
#42:616647,
#43:795912,
#44:1080116,
#45:1354579,
#46:94105,
#47:94562,
#48:581380,
#49:650016,
#50:795405,
#51:795648,
#52:831279,
#53:1104787,
#54:1105005,
#55:1161589,
#56:1194734,
#57:1194785,
#58:1354730,
#59:1422249,
#60:1476003,
#61:1507285,
#62:1697899,
#63:1725470,
#64:1764462,
#65:1940047}

hid2name = {}
hid2catnum = {}
for k,v in cid2hid.items():
    hid2name[v] = 'Cat-'+str(k)
    hid2catnum[v] = k

#for k, v in cid2hid.iteritems():
#    hid2hpath[v].append(k)

# functions for switching between hid,hpath,catnum
def catnum_hpath(catnum,lx):
    hid = cid2hid[catnum]
    hpath = hid_hpath_lx(hid,lx)
    return hpath
def hid_name(hid):
    return hid2name[hidint(hid)]
def hpath_name(hpath):
    return hid_name(get_parent_hid(hpath))
def hid_catnum(hid):
    return hid2catnum[hidint(hid)]
def hpath_catnum(hpath):
    return hid_catnum(get_parent_hid(hpath))
def hid_hpath_lx(hid,do_lx):
    lxpaths = get_lxlist(hid,gethpaths=True)
    for hpath in lxpaths:
        if 'LX'+str(do_lx) in hpath: return hpath
    return None
def get_paper_paths_lx(do_lx):
    return [hid_hpath_lx(cid2hid[i+1],do_lx) for i in np.arange(24)] #[hid_hpath_lx(hid,do_lx) for hid in hid2name.keys()]
def get_all_halo_paths_lx(lx):
    hpaths = []
    for i in np.arange(60):
        try:
            path = hid_hpath_lx(cid2hid[i+1],lx)
            if path != None: hpaths.append(path)
        except KeyError:
            pass
    return hpaths
def get_paper_paths():
    return [global_halobase+"/H"+str(hid) for hid in hid2name.keys()]
def get_good_paper_paths():
    hpaths = [global_halobase+"/H"+str(hid) for hid in hid2name.keys()]
    hpaths.remove(global_halobase+"/H95289")
    return hpaths

def hidint(hid):
    """ converts halo ID to int """
    if type(hid)==int or type(hid)==np.int64: return hid
    if type(hid)==str:
        if hid[0]=='H': return int(hid[1:])
        return int(hid)
    raise ValueError("hid must be int or str, is "+str(type(hid)))
def hidstr(hid):
    """ converts halo ID to str Hxxxxxx """
    if type(hid)==int or type(hid)==np.int64: return 'H'+str(hid)
    if type(hid)==str:
        if hid[0]=='H': return hid
        return 'H'+hid
    raise ValueError("hid must be int or str, is "+str(type(hid)))

def get_parent_zoom_index(filename=global_halobase+"/parent_zoom_index.txt"):
    return asciitable.read(filename, Reader=asciitable.FixedWidth)
def get_numsnaps(outpath):
    """
    Uses hpath/ExpansionList to get the number of snaps in this halo
    """
    if os.path.exists(outpath+'/ExpansionList'):
        return sum(1 for line in open(outpath+'/ExpansionList'))
    else:
        warnings.warn(outpath+"/ExpansionList not found, using default (256)")
        return 256
def get_lastsnap(outpath):
    return get_numsnaps(outpath)-1
def get_foldername(outpath):
    return os.path.basename(os.path.normpath(outpath))
def get_parent_hid(outpath):
    hidstr = get_foldername(outpath).split('_')[0]
    return int(hidstr[1:])
def get_contamtype(outpath):
    return get_zoom_params(outpath)[0]
def get_zoom_params(outpath):
    """ return ictype, LX, NV """
    split = get_foldername(outpath).split('_')
    return split[1],int(split[5][2:]),int(split[7][2:])
def get_outpath(haloid,ictype,lx,nv,contamtype=None,halobase=global_halobase,check=True):
    haloid = hidstr(haloid); ictype = ictype.upper()
    outpath = halobase+'/'+haloid+'/'+haloid+'_'+ictype+'_'+'Z127_P7_LN7_LX'+str(lx)+'_O4_NV'+str(nv)
    if contamtype != None:
        outpath += '_'+str(contamtype)
    if check and not os.path.exists(outpath):
        raise IOError("Invalid hpath")
    return outpath
def get_hpath(haloid,ictype,lx,nv,contamtype=None,halobase=global_halobase,check=True):
    return get_outpath(haloid,ictype,lx,nv,contamtype=contamtype,halobase=global_halobase,check=check)


def get_hpath_lx(hid,do_lx):
    lxpaths = get_lxlist(hid,gethpaths=True)
    for hpath in lxpaths:
        if 'LX'+str(do_lx) in hpath: return hpath
    return None


def get_scale_snap(hpath,snaps):
    snaps = np.ravel(snaps)
    if type(snaps[0]) == np.float64 or type(snaps[0])==np.float:
        badii = np.isnan(snaps)
    else:
        badii = snaps < 0
    goodii = ~badii
    assert np.all(snaps[goodii].astype(int)==snaps[goodii]),'Snaps must be integers'
    snaps = snaps.astype(int)
    snaps[badii] = -1

    numsnaps = get_numsnaps(hpath)
    assert np.all(snaps[goodii] < numsnaps) and np.all(snaps[goodii] >= 0), "Snaps must be between 0 and {0}".format(numsnaps-1)
    with open(hpath+'/ExpansionList','r') as f:
        lines = f.readlines()
    #assert(len(lines))==numsnaps
    def _get_scale_snap(snap):
        if snap == -1: return np.nan
        return float(lines[snap].split()[0])
    return np.array(map(_get_scale_snap,snaps))
def get_z_snap(hpath,snap):
    scale = get_scale_snap(hpath,snap)
    return (1./scale) - 1.0
def get_t_snap(hpath,snap,OmegaM=.3175,OmegaL=.6825,h=.6711):
    scale = get_scale_snap(hpath,snap)
    return bconversions.GetTime(scale,OmegaM=OmegaM,OmegaL=OmegaL,h=h)
def get_snap_z(hpath, zs):
    scales = get_scale_list(hpath)
    snaps = [np.argmin(np.abs(scales-1./(1+z))) for z in np.ravel(zs)]
    return np.array(snaps)
def get_scale_list(hpath):
    with open(hpath+'/ExpansionList','r') as f:
        lines = f.readlines()
    scales = np.array([float(line.split()[0]) for line in lines])
    return scales

def get_available_hpaths(hid,contam=False,
                         checkgadget=True,
                         onlychecklastsnap=True,
                         checkallblocks=False,
                         hdf5=True,verbose=False,
                         basepath=global_halobase,
                         hires=False):
    hidpath = basepath+'/'+hidstr(hid)
    if contam: hidpath += '/contamination_suite'
    if not os.path.exists(hidpath):
        raise IOError("Invalid hid: "+hidpath)
    hpathlist = []
    for foldername in os.listdir(hidpath):
        if foldername[0] != 'H': continue
        hpath = hidpath+'/'+foldername
        if not os.path.isdir(hpath): continue
        try:
            if checkgadget and not gadget_finished(hpath,onlychecklastsnap=onlychecklastsnap,
                                                   checkallblocks=checkallblocks,hdf5=hdf5,
                                                   verbose=verbose,hires=hires):
                continue
        except IOError:
            continue
        hpathlist.append(hpath)
    return hpathlist
def get_lxlist(hid,gethpaths=False):
    outlist = []
    availablehpaths = get_available_hpaths(hid)
    for lx in [14,13,12,11]:
        for hpath in availablehpaths:
            if 'LX'+str(lx) in hpath: 
                if gethpaths:
                    outlist.append(hpath)
                else: 
                    outlist.append(lx)
    assert len(outlist)==len(np.unique(outlist))
    return outlist

def check_last_subfind_exists(outpath):
    numsnaps = get_numsnaps(outpath)
    lastsnap = numsnaps - 1; snapstr = str(lastsnap).zfill(3)
    group_tab = os.path.exists(outpath+'/outputs/groups_'+snapstr+'/group_tab_'+snapstr+'.0')
    subhalo_tab = os.path.exists(outpath+'/outputs/groups_'+snapstr+'/subhalo_tab_'+snapstr+'.0')
    return group_tab and subhalo_tab

def check_rockstar_exists(outpath,snap,boundbin=True,fullbin=False,particles=False):
    snapstr = str(snap)
    if fullbin:
        halo_exists = os.path.exists(outpath+'/halos/halos_'+snapstr+'/halos_'+snapstr+'.0.fullbin')
    elif boundbin:
        halo_exists = os.path.exists(outpath+'/halos_bound/halos_'+snapstr+'/halos_'+snapstr+'.0.boundbin')
    else:
        halo_exists = os.path.exists(outpath+'/halos/halos_'+snapstr+'/halos_'+snapstr+'.0.bin')
    if not particles:
        return halo_exists
    part_exists = os.path.exists(outpath+'/halos/halos_'+snapstr+'/halos_'+snapstr+'.0.particles')
    return halo_exists and part_exists

def check_last_rockstar_exists(outpath,boundbin=True,fullbin=False,particles=False):
    numsnaps = get_numsnaps(outpath)
    lastsnap = numsnaps - 1; snapstr = str(lastsnap)
    return check_rockstar_exists(outpath,lastsnap,boundbin=boundbin,fullbin=fullbin,particles=particles)

def check_mergertree_exists(outpath,autoconvert=False,boundbin=True,treedir='trees'):
    if boundbin: halodir = 'halos_bound'
    else: halodir = 'halos'
    ascii_exists = os.path.exists(outpath+'/'+halodir+'/'+treedir+'/tree_0_0_0.dat')
    binary_exists = os.path.exists(outpath+'/'+halodir+'/'+treedir+'/tree.bin')
    if autoconvert and ascii_exists and not binary_exists:
        print "---check_mergertree_exists: Automatically converting ascii to binary"
        MTC.convertmt(outpath+'/'+halodir+'/'+treedir,version=4)
        binary_exists = os.path.exists(outpath+'/'+halodir+'/'+treedir+'/tree.bin')
    return ascii_exists and binary_exists

def check_is_sorted(outpath,snap=0,hdf5=True):
    #TODO: option to check all snaps
    snap = str(snap).zfill(3)
    filename = outpath+'/outputs/snapdir_'+snap+'/snap_'+snap+'.0'
    if hdf5: filename += '.hdf5'
    h = rsg.snapshot_header(filename)
    try:
        if h.sorted=='yes': return True
    except:
        return False

def gadget_finished(outpath,
                    onlychecklastsnap=False,
                    checkallblocks=False,
                    hdf5=True,verbose=False,
                    hires=False):
    numsnaps = get_numsnaps(outpath)
    gadgetpath = outpath+'/outputs'
    if hires: 
        numsnaps = 320
        gadgetpath = outpath+'/outputs_hires'
    if (not os.path.exists(gadgetpath)):
        if verbose: print "  Gadget folder not present in "+get_foldername(outpath)
        return False
    if onlychecklastsnap: #only check last snap
        snapstr = str(numsnaps-1).zfill(3)
        snappath = gadgetpath+"/snapdir_"+snapstr+"/snap_"+snapstr+".0"
        if hdf5: snappath += ".hdf5"
        if (not os.path.exists(snappath)):
            if verbose: print "  Snap "+snapstr+" not in "+get_foldername(outpath)
            return False
        else:
            return True
    for snap in xrange(numsnaps): # check that all snaps are there
        snapstr = str(snap).zfill(3)
        snappath = gadgetpath+"/snapdir_"+snapstr+"/snap_"+snapstr+".0"
        if hdf5: snappath += ".hdf5"
        if (not os.path.exists(snappath)):
            if verbose: print "  Snap "+snapstr+" not in "+get_foldername(outpath)
            return False
        if checkallblocks:
            for snapfile in glob.glob(gadgetpath+"/snapdir_"+snapstr+'/*'):
                if (os.path.getsize(snapfile) <= 0):
                    if verbose: print snapfile,"has no data (skipping)"
                    return False
    return True

def restrict_halopaths(halopathlist,
                       require_rockstar=False,
                       require_subfind=False,
                       require_sorted=False,
                       require_mergertree=False,
                       autoconvert_mergertree=False,
                       use_fullbin_rockstar=False):
    if require_rockstar:
        newhalopathlist = []
        for outpath in halopathlist:
            if use_fullbin_rockstar:
                if check_last_rockstar_exists(outpath,fullbin=True,boundbin=False):
                    newhalopathlist.append(outpath) 
            else:
                if check_last_rockstar_exists(outpath,boundbin=True):
                    newhalopathlist.append(outpath) 
        halopathlist = newhalopathlist
    if require_subfind:
        newhalopathlist = []
        for outpath in halopathlist:
            if check_last_subfind_exists(outpath):
                newhalopathlist.append(outpath) 
        halopathlist = newhalopathlist
    if require_sorted:
        newhalopathlist = []
        for outpath in halopathlist:
            if check_is_sorted(outpath,snap=get_numsnaps(outpath)-1):
                newhalopathlist.append(outpath) 
        halopathlist = newhalopathlist
    if require_mergertree:
        newhalopathlist = []
        for outpath in halopathlist:
            if check_mergertree_exists(outpath,autoconvert=autoconvert_mergertree):
                newhalopathlist.append(outpath) 
        halopathlist = newhalopathlist
    return halopathlist

def find_halo_paths(basepath=global_halobase,
                    nrvirlist=[3,4,5,6],levellist=[11,12,13,14],
                    ictypelist=["BA","BB","BC","BD","EA","EB","EC","EX","CA","CB","CC"],
                    contamsuite=False,
                    require_rockstar=False,require_subfind=False,
                    require_mergertree=False,autoconvert_mergertree=False,
                    require_sorted=False,
                    checkallblocks=False,
                    onlychecklastsnap=False,verbose=False,hdf5=True,
                    use_fullbin_rockstar=False,hires=False):
    """ Returns a list of paths to halos that have gadget completed/rsynced
        with the specified nrvirlist/levellist/ictype """
    if verbose:
        print "basepath:",basepath
        print "nrvirlist:",nrvirlist
        print "levellist:",levellist
        print "ictypelist:",ictypelist

    halopathlist = []
    haloidlist = []
    for filename in os.listdir(basepath):
        if filename[0] == "H":
            haloidlist.append(filename)
    for haloid in haloidlist:
        try:
            hpathlist = get_available_hpaths(haloid,contam=contamsuite,checkgadget=False,
                                             basepath=basepath)
        except IOError as e:
            print "ERROR: skipping",haloid
            continue
        for hpath in hpathlist:
            ictype,levelmax,nrvir = get_zoom_params(hpath)
            if (int(levelmax) in levellist and int(nrvir) in nrvirlist and ictype in ictypelist):
                try:
                    if gadget_finished(hpath,
                                       onlychecklastsnap=onlychecklastsnap,
                                       checkallblocks=checkallblocks,
                                       hdf5=hdf5,verbose=verbose,
                                       hires=hires):
                        halopathlist.append(hpath)
                except IOError as e:
                    print "ERROR: skipping",hpath

    halopathlist = restrict_halopaths(halopathlist,
                                      require_rockstar=require_rockstar,
                                      require_subfind=require_subfind,
                                      require_sorted=require_sorted,
                                      require_mergertree=require_mergertree,
                                      autoconvert_mergertree=autoconvert_mergertree,
                                      use_fullbin_rockstar=use_fullbin_rockstar)
    return halopathlist

def _load_index_row(hpath,filename=global_halobase+"/parent_zoom_index.txt"):
    haloid = get_parent_hid(hpath)
    ictype,lx,nv = get_zoom_params(hpath)
    htable = get_parent_zoom_index()
    haloid = hidint(haloid); lx = int(lx); nv = int(nv)

    idmask = htable['parentid']==haloid
    icmask = htable['ictype']==ictype.upper()
    lxmask = htable['LX']==lx
    nvmask = htable['NV']==nv
    maskall = idmask & icmask & lxmask & nvmask
    if np.sum(maskall) == 0:
        raise ValueError("no such halo in index for %s" % (hpath))
    if np.sum(maskall) > 1:
        print "FATAL ERROR: duplicate row in index"
        exit()
    row = htable[maskall]
    if row['badflag']+row['badsubf'] > 0:
        if (lx != 14) or (lx==14 and row['badflag']>0):
            print "WARNING: potentially bad halo match for H%i %s LX%i NV%i" % (haloid,ictype,lx,nv)
    return row
def load_zoomid(hpath,filename=global_halobase+"/parent_zoom_index.txt",snap=None):
    """
    @param hpath: halo path to load zoom id
    @param snap: default None (automatically picks snap with get_numsnaps)
    @return: rockstar id of host halo associated with hpath and snap

    IMPORTANT: Uses MassAccrPlugin to get main branch MT for snap != last snap.
    """
    if snap==None:
        snap = get_numsnaps(hpath)-1
    if snap==(get_numsnaps(hpath)-1):
        try:
            row = _load_index_row(hpath,filename=filename)
        except ValueError:
            if check_last_rockstar_exists(hpath):
                print "WARNING: halo is not in index, using halo with most particles (npart)"
                rscat = load_rscat(hpath,get_numsnaps(hpath),rmaxcut=False)
                # For pandas < 0.13.0, np.argmax returns the array index rather than the pandas index
                pdversion = tuple([int(x) for x in pd.version.version.split('.')])
                badversion = (0,13,0)
                bestid = np.argmax(rscat['npart'])
                if pdversion < badversion:
                    bestid = rscat.data.index[bestid]
                return bestid
            else:
                raise ValueError("No rockstar catalogue for {0}!".format(get_foldername(hpath)))
        return row['zoomid'][0]
    else:
        from caterpillaranalysis import MassAccrPlugin
        plug = MassAccrPlugin()
        tab = plug.read(hpath)
        snaplist = tab['snap']
        ii = snaplist == snap
        if np.sum(ii) != 1: 
            raise ValueError("{0} snap {1} does not have a valid main branch rsid ({2} indices match snap)".format(get_foldername(hpath),snap,np.sum(ii)))
        return tab['origid'][ii][0]
    
def load_haloprops(hpath,filename=global_halobase+"/parent_zoom_index.txt"):
    row = _load_index_row(hpath,filename=filename)
    mvir = float(row['mgrav']) # physical Msun
    rvir = float(row['rvir'])  # physical kpc
    vvir = np.sqrt(4.34e-6 * mvir/rvir) # physical km/s
    return mvir,rvir,vvir

def load_pcatz0(old=False):
    if old:
        return RDR.RSDataReader(global_basepath+"/caterpillar/parent/RockstarData",63,version=2)
    else:
        return RDR.RSDataReader(global_prntbase+"/rockstar",127,version=6)

def load_scat(hpath):
    snap = get_lastsnap(hpath)
    if "LX14" in hpath:
        try:
            return RSF.subfind_catalog(hpath+'/outputs',snap,double=True)
        except ValueError:
            return RSF.subfind_catalog(hpath+'/outputs',snap)
    else: 
        return RSF.subfind_catalog(hpath+'/outputs',snap)

def load_rscat(hpath,snap,verbose=True,halodir='halos_bound',unboundfrac=None,minboundpart=None,version=None,rmaxcut=True,AllParticles=False):
    if version != None:
        rscat = RDR.RSDataReader(hpath+'/'+halodir,snap,version=version,digits=1,unboundfrac=unboundfrac,minboundpart=minboundpart,AllParticles=AllParticles)
    else:
        try:
            rscat = RDR.RSDataReader(hpath+'/'+halodir,snap,version=10,digits=1,unboundfrac=unboundfrac,minboundpart=minboundpart,AllParticles=AllParticles)
        except IOError as e: #try to identify a unique valid rockstar version
            print e
            versionlist = [2,3,4,5,6,7,8,9]
            testlist = []
            for version in versionlist:
                try:
                    rscat = RDR.RSDataReader(hpath+'/'+halodir,snap,version=version,digits=1,unboundfrac=unboundfrac,minboundpart=minboundpart,AllParticles=AllParticles)
                    testlist.append(True)
                except KeyError:
                    testlist.append(False)
                except IOError:
                    testlist.append(False)
            if sum(testlist) != 1:
                raise RuntimeError("Can't determine what version to use {0}".format(get_foldername(hpath)))
            else:
                version = np.array(versionlist)[np.array(testlist)][0]
                if verbose:
                    print "Using version "+str(version)+" for "+get_foldername(hpath)
                    rscat = RDR.RSDataReader(hpath+'/'+halodir,snap,version=version,digits=1,unboundfrac=unboundfrac,minboundpart=minboundpart,AllParticles=AllParticles)

    if rmaxcut:
        zoomid = load_zoomid(hpath,snap=snap)
        hpos = np.array(rscat.ix[zoomid][['posX','posY','posZ']])
        spos = rscat[['posX','posY','posZ']]
        dr = np.sqrt(np.sum((spos-hpos)**2,1))
        rscat['dr'] = dr*1000.
        badii = rscat['dr']<rscat['rvmax']
        badii[zoomid] = False
        goodii = ~badii
        rscat.badhalos = rscat.data.ix[badii]
        rscat.numbad = len(rscat.badhalos)
        rscat.data = rscat.data.ix[goodii]
        rscat.ix = rscat.data.ix
        rscat.index = rscat.data.index
        rscat.num_halos = len(rscat.data)
    return rscat

def load_rsboundindex(hpath,snap):
    return RDR.load_rsboundindex(hpath,snap)
    
def load_mtc_general(hpath,version,verbose=True,halodir='halos_bound',treedir='trees',**kwargs):
    return MTC.MTCatalogue(hpath+'/'+halodir+'/'+treedir,version=version,**kwargs)
def load_mtc(hpath,verbose=True,halodir='halos_bound',treedir='trees',**kwargs):
    return MTC.MTCatalogue(hpath+'/'+halodir+'/'+treedir,version=4,**kwargs)
def load_zoom_mtc(hpath,verbose=True,halodir='halos_bound',treedir='trees',**kwargs):
    return MTC.MTCatalogue(hpath+'/'+halodir+'/'+treedir,version=4,haloids=[load_zoomid(hpath)],**kwargs)
def make_mtindex_key(snap,origid):
    return str(snap).zfill(3)+","+str(origid)
def make_mt_snapid_to_baseidrow(hpath,recalc=False,halodir='halos_bound',treedir='trees'):
    def make_key(arg):
        return make_mtindex_key(arg[0],arg[1])
    def make_val(arg):
        return arg
    #def make_mtindex_val(base_rsid,row):
    #    return (base_rsid,row)
    indexpath = hpath+'/'+halodir+'/'+treedir+'/snap_id_to_baseid_row.p'
    if os.path.exists(indexpath) and (not recalc): return

    start = time.time()
    mtc = load_mtc(hpath,indexbyrsid=True)
    allkeys = []
    allvals = []
    print "Load Time: {0:.1f} sec".format(time.time()-start)

    start = time.time()
    for base_rsid,mt in mtc.Trees.iteritems():
        snaps = mt['snap']
        origids = mt['origid']
        keys = map(make_key,zip(snaps,origids))
        vals = map(make_val,zip(itertools.repeat(base_rsid,len(keys)),itertools.count()))
        
        allkeys += keys
        allvals += vals
    
    index = dict(zip(allkeys,allvals))
    with open(indexpath,'w') as f:
        pickle.dump(index,f)
    print "Convert Time: {0:.1f} sec".format(time.time()-start)
    subprocess.call(['chmod g+rwx '+indexpath],shell=True)
    subprocess.call(['chgrp annaproj '+indexpath],shell=True)
def load_mt_snapid_to_baseidrow(hpath,halodir='halos_bound',treedir='trees'):
    indexpath = hpath+'/'+halodir+'/'+treedir+'/snap_id_to_baseid_row.p'
    assert os.path.exists(indexpath)
    with open(indexpath,'r') as f:
        index = pickle.load(f)
    return index

def load_pmtc(hpath=global_prntbase,verbose=True,halodir='rockstar',treedir='trees',**kwargs):
    return MTC.MTCatalogue(hpath+'/'+halodir+'/'+treedir,version=3,**kwargs)

def load_partblock(hpath,snap,block,parttype=-1,ids=-1,hdf5=True):
    #assert check_is_sorted(hpath,snap=snap,hdf5=hdf5),"snap is sorted"
    snapstr = str(snap).zfill(3)
    snappath = hpath+'/outputs/snapdir_'+snapstr+'/snap_'+snapstr
#    if "14" in hpath:
#        return rsg.read_block(snappath,block,parttype=parttype,ids=ids,doubleprec=True)
#    else:
    return rsg.read_block(snappath,block,parttype=parttype,ids=ids)

def load_soft(hpath):
    """ plummer equivalent grav. softening = h/2.8 """
    try:
        fname = hpath+'/param.txt-usedvalues'
        if not os.path.exists(fname): raise IOError("Could not find file "+fname)
        forceres=-1
        f = open(fname,'r')
        for line in f:
            s = line.split()
            if s[0]=="SofteningHaloMaxPhys" or s[0]=="SofteningMaxPhysType1":
                forceres = float(s[1])
                break
        f.close()
        if forceres==-1: raise IOError("Could not find force resolution")
    except IOError as e:
        print "WARNING:",e
        ictype,lx,nv = get_zoom_params(hpath)
        forceres = 100./2.^lx/80.
    return forceres

def load_aqcat(whichAq,snap):
    assert whichAq in ['A','B','C','D','E','F']
    if snap > 127: 
        raise ValueError("Aquarius is snaps 0-127")
    rspath = global_basepath+'/aquarius/Aq-'+whichAq+'/2/halos'
    return RDR.RSDataReader(rspath,snap,version=7)

def get_quant_zoom(halo_path,quant):
    htable = get_parent_zoom_index()
    halo_split = halo_path.split("_")
    haloid = int(halo_path.split("/H")[-1].split("_")[0].strip("H"))
    geom,lx,nrvir = get_zoom_params(halo_path.split("/")[-1])
    mask = (haloid == htable['parentid']) & \
           (geom == htable['ictype']) & \
           (int(lx) == htable['LX']) & \
           (int(nrvir) == htable['NV']) 

    if len(htable[mask][quant])>1:
        return htable[mask][quant]
    else:
        return float(htable[mask][quant])

def get_main_branch(hpath):
    return pickle.load( open( hpath+"/analysis/main_branch.p", "rb" ) )

def get_halo_header(hpath,snap=None):
    if snap==None: snap = get_lastsnap(hpath)
    return rsg.snapshot_header(hpath+"/outputs/snapdir_"+str(snap)+"/snap_"+str(snap))

def get_colors(ncolors=12):
    colors = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
              (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
              (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
              (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
              (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),
              (32, 159, 117), (43, 166, 241),  (6, 115,   1), (203,  56,  62), (255,255,255)] 
    assert len(colors)==len(hid2name)
    
    for i in range(ncolors):  
        r, g, b = colors[i]  
        colors[i] = (r / 255., g / 255., b / 255.)  

    return colors

def get_halo_colors(cmap='jet',nhalos=24):
    _numfirstbatch = 24
    #colors = sns.color_palette("husl",_numfirstbatch)
    #print cm
    cm = plt.cm.get_cmap(cmap)
    ncolors = 24
    colors=[]
    for i in range(ncolors):
        colors.append(cm(1.*i/float(ncolors)))

    index = {}
#    for i in range(nhalos):
    for i,haloi in cid2hid.iteritems():
#        index[hid2name.keys()[i]] = colors[i]
	index[haloi] = colors[i-1]

    return index
    

def get_colors_for_halos(nhalos=len(hid2name)):
    colors = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
              (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
              (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
              (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
              (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),
              (32, 159, 117), (43, 166, 241),  (6, 115,   1), (203,  56,  62)]
#, (255,255,255)] 
    assert len(colors)==len(hid2name)
    
    for i in range(nhalos):  
        r, g, b = colors[i]  
        colors[i] = (r / 255., g / 255., b / 255.)  

    index = {}
    for i in range(nhalos):
        index[hid2name.keys()[i]] = colors[i]

    return index

def tabulate(tabfn,lx=14,hids=None,exclude_hids=None,savefile=None,numprocs=1,usecid=False):
    """
    @param tabfn: a function whose only argument is hpath. If successful, returns data,names,formats; data is a tuple of the values to go into the array, names is a list/tuple of the column names, formats is the data types (used for np.dtype). All three variables should be of the same length. If unsuccessful, tabfn should return None (e.g., in cases of no rockstar data), then tabulate() will make the DataFrame row be marked with missing data. TODO can I make this happen with exceptions in a nice way? Raising exceptions is better than returning None.
    @param lx: which LX to tabulate (default 14)
    @param hids: list of hids to tabulate (default, everything in cid2hid)
    @param exclude_hids: list of hids to exclude
    @param savefile: name of file to save df as a csv to
    @param numprocs: if larger than 1, uses multiprocessing.Pool.map() to tabulate
    @param usecid: DataFrame index is by Cat-ID rather than hid
    @return tab: pandas DataFrame, indexed by hid (or Cat-ID if usecid)
    """
    if hids==None: hids = cid2hid.values()
    if exclude_hids != None:
        for ex_hid in exclude_hids:
            ex_hid = hidint(ex_hid)
            if ex_hid in hids: 
                hids.remove(ex_hid)
            else:
                print "WARNING: H{0} not in hids, not removing"

    if numprocs==1:
        datalist = map(tabfn,[hid_hpath_lx(hid,lx) for hid in hids])
    else:
        pool = Pool(numprocs)
        datalist = pool.map(tabfn,[hid_hpath_lx(hid,lx) for hid in hids])
        pool.close()

    for item in datalist:
        if item != None:
            data,names,formats = item
            first_dtype = np.dtype({'names':names,'formats':formats})
            break
    rows = []
    invalid = [False for item in datalist]
    dtype = first_dtype
    for i,item in enumerate(datalist):
        if item == None:
            data = tuple([0 for x in dtype.names])
            invalid[i] = True
        else:
            assert len(item)==3,"tabfn must return data,colnames,types"
            data,names,formats = item
            assert len(data)==len(names) and len(data)==len(formats),"data,names,formats must be same length"
            data = tuple(data)
            dtype = np.dtype({'names':names,'formats':formats})
            assert dtype==first_dtype,"Not all dtypes the same"
        rows.append(data)

    dataarr = np.array(rows,dtype=dtype)

    if usecid:
        myindex = [hid2catnum[hid] for hid in hids]
    else:
        myindex = hids

    df = pd.DataFrame(dataarr,index=myindex)
    for i,hid in enumerate(df.index):
        if invalid[i]: df.iloc[i] = None

    if savefile != None:
        df.to_csv(path_or_buf=savefile)

    return df
