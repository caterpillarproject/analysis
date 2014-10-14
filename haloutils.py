import numpy as np
import os,sys,platform
import subprocess
import asciitable

import readsnapshots.readsnapHDF5_greg as rsg
import readhalos.RSDataReader as RDR
import readhalos.readsubf as RSF
import mergertrees.MTCatalogue as MTC
from brendanlib.grifflib import determinebasepath

global_basepath = determinebasepath(platform.node())
global_halobase = global_basepath+'/caterpillar/halos'
global_prntbase = global_basepath+'/caterpillar/parent/gL100X10'

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
    return sum(1 for line in open(outpath+'/ExpansionList'))
def get_foldername(outpath):
    return os.path.basename(os.path.normpath(outpath))
def get_parent_hid(outpath):
    hidstr = get_foldername(outpath).split('_')[0]
    return int(hidstr[1:])
def get_contamtype(outpath):
    contamtype = get_foldername(outpath).split('_')[-1]
    if 'NV' in contamtype: raise ValueError("outpath has no contamtype")
    return contamtype
def get_zoom_params(outpath):
    """ return ictype, LX, NV """
    split = get_foldername(outpath).split('_')
    return split[1],int(split[5][2:]),int(split[7][2:])
def get_outpath(haloid,ictype,lx,nv,contamtype=None,halobase=global_halobase):
    haloid = hidstr(haloid); ictype = ictype.upper()
    if contamtype==None:
        return halobase+'/'+haloid+'/'+haloid+'_'+ictype+'_'+'Z127_P7_LN7_LX'+str(lx)+'_O4_NV'+str(nv)
    return halobase+'/'+haloid+'/'+haloid+'_'+ictype+'_'+'Z127_P7_LN7_LX'+str(lx)+'_O4_NV'+str(nv)+'_'+str(contamtype)
def get_hpath(haloid,ictype,lx,nv,contamtype=None,halobase=global_halobase):
    return get_outpath(haloid,ictype,lx,nv,contamtype=contamtype,halobase=global_halobase)

def check_last_subfind_exists(outpath):
    numsnaps = get_numsnaps(outpath)
    lastsnap = numsnaps - 1; snapstr = str(lastsnap).zfill(3)
    group_tab = os.path.exists(outpath+'/outputs/groups_'+snapstr+'/group_tab_'+snapstr+'.0')
    subhalo_tab = os.path.exists(outpath+'/outputs/groups_'+snapstr+'/subhalo_tab_'+snapstr+'.0')
    return group_tab and subhalo_tab

def check_last_rockstar_exists(outpath,fullbin=True,particles=False):
    numsnaps = get_numsnaps(outpath)
    lastsnap = numsnaps - 1; snapstr = str(lastsnap)
    if fullbin:
        halo_exists = os.path.exists(outpath+'/halos/halos_'+snapstr+'/halos_'+snapstr+'.0.fullbin')
    else:
        halo_exists = os.path.exists(outpath+'/halos/halos_'+snapstr+'/halos_'+snapstr+'.0.bin')
    if not particles:
        return halo_exists
    part_exists = os.path.exists(outpath+'/halos/halos_'+snapstr+'/halos_'+snapstr+'.0.particles')
    return halo_exists and part_exists

def check_mergertree_exists(outpath,autoconvert=False):
    ascii_exists = os.path.exists(outpath+'/trees/tree_0_0_0.dat')
    binary_exists = os.path.exists(outpath+'/trees/tree.bin')
    if not binary_exists and autoconvert:
        print "---check_mergertree_exists: Automatically converting ascii to binary"
        MTC.convertmt(outpath+'/trees',version=4)
        binary_exists = os.path.exists(outpath+'/trees/tree.bin')
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

def find_halo_paths(basepath=global_halobase,
                    nrvirlist=[3,4,5,6],levellist=[11,12,13,14],ictype="BB",
                    contamsuite=False,
                    require_rockstar=False,require_subfind=False,
                    require_mergertree=False,autoconvert_mergertree=False,
                    require_sorted=False,
                    checkallblocks=False,
                    onlychecklastsnap=False,verbose=False,hdf5=True):
    """ Returns a list of paths to halos that have gadget completed/rsynced
        with the specified nrvirlist/levellist/ictype """
    if verbose:
        print "basepath:",basepath
        print "nrvirlist:",nrvirlist
        print "levellist:",levellist
        print "ictype:",ictype
    def gadget_finished(outpath):
        numsnaps = get_numsnaps(outpath)
        gadgetpath = outpath+'/outputs'
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

    halopathlist = []
    haloidlist = []
    for filename in os.listdir(basepath):
        if filename[0] == "H":
            haloidlist.append(filename)
    for haloid in haloidlist:
        if contamsuite:
            subdirnames = basepath + "/" + haloid + '/contamination_suite'
        else:
            subdirnames = basepath + "/" + haloid
        halosubdirlist = []
        try:
            for filename in os.listdir(subdirnames):
                if filename=='contamination_suite': continue
                halosubdirlist.append(filename)
                thisictype,levelmax,nrvir = get_zoom_params(filename)
                haloid = hidstr(get_parent_hid(filename))
                if (int(levelmax) in levellist and int(nrvir) in nrvirlist and ictype==thisictype):
                    if contamsuite:
                        outpath = basepath+"/"+haloid+"/contamination_suite/"+filename
                    else:
                        outpath = basepath+"/"+haloid+"/"+filename
                    try:
                        if gadget_finished(outpath): halopathlist.append(outpath)
                    except IOError as e:
                        print "ERROR: skipping",outpath
                        continue
        except IOError as e:
            print "IOError: skipping",subdirnames
            continue
        except OSError as e:
            print "OSError: skipping",subdirnames
            continue

    if require_rockstar:
        newhalopathlist = []
        for outpath in halopathlist:
            if check_last_rockstar_exists(outpath):
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
            if check_is_sorted(outpath,snap=255):
                newhalopathlist.append(outpath) 
        halopathlist = newhalopathlist
    if require_mergertree:
        newhalopathlist = []
        for outpath in halopathlist:
            if check_mergertree_exists(outpath,autoconvert=autoconvert_mergertree):
                newhalopathlist.append(outpath) 
        halopathlist = newhalopathlist
    return halopathlist

def load_zoomid(hpath,filename=global_halobase+"/parent_zoom_index.txt"):
    haloid = get_parent_hid(hpath)
    ictype,lx,nv = get_zoom_params(hpath)
    htable = get_parent_zoom_index()
    haloid = hidint(haloid); lx = int(lx); nv = int(nv)
    #if lx==14 and haloid==1327707: return 188661
    if lx==14 and haloid==706754: return 269650
    if lx==14 and haloid==649524: return 481480
    if lx==14 and haloid==1725139: return 135325

    idmask = htable['parentid']==haloid
    icmask = htable['ictype']==ictype.upper()
    lxmask = htable['LX']==lx
    nvmask = htable['NV']==nv
    maskall = idmask & icmask & lxmask & nvmask
    if np.sum(maskall) == 0:
        raise ValueError("no such halo in index")
    if np.sum(maskall) > 1:
        print "FATAL ERROR: duplicate row in index"
        exit()
    row = htable[maskall]
    if row['badflag']+row['badsubf'] > 0:
        print "WARNING: potentially bad halo match for H%i %s LX%i NV%i" % (haloid,ictype,lx,nv)
    return row['zoomid'][0]

def load_pcatz0(old=False):
    if old:
        return RDR.RSDataReader(global_basepath+"/caterpillar/parent/RockstarData",63,version=2)
    else:
        return RDR.RSDataReader(global_prntbase+"/rockstar",127,version=6)

def load_scat(hpath):
    return RSF.subfind_catalog(hpath+'/outputs',255)

def load_rscat(hpath,snap,verbose=True):
    try:
        rcat = RDR.RSDataReader(hpath+'/halos',snap,version=7,digits=1)
    except IOError:
        versionlist = [2,3,4,5,6,7]
        testlist = []
        for version in versionlist:
            try:
                rcat = RDR.RSDataReader(hpath+'/halos',snap,version=version)
                testlist.append(True)
            except KeyError:
                testlist.append(False)
        if sum(testlist) != 1:
            raise RuntimeError("Can't determine what version to use")
        else:
            version = np.array(versionlist)[np.array(testlist)][0]
            if verbose:
                print "Using version "+str(version)+" for "+get_foldername(hpath)
            rcat = RDR.RSDataReader(hpath+'/halos',snap,version=version)
    return rcat

def load_mtc(hpath,verbose=True,**kwargs):
    return MTC.MTCatalogue(hpath+'/halos/trees',version=4,**kwargs)

def load_partblock(hpath,snap,block,parttype=-1,ids=-1,hdf5=True):
    #assert check_is_sorted(hpath,snap=snap,hdf5=hdf5),"snap is sorted"
    snapstr = str(snap).zfill(3)
    snappath = hpath+'/outputs/snapdir_'+snapstr+'/snap_'+snapstr
    return rsg.read_block(snappath,block,parttype=parttype,ids=ids)

def load_soft(hpath):
    try:
        fname = hpath+'/param.txt-usedvalues'
        if not os.path.exists(fname): raise IOError("Could not find file "+fname)
        forceres=-1
        f = open(fname,'r')
        for line in f:
            s = line.split()
            if s[0]=="SofteningHaloMaxPhys":
                forceres = float(s[1])
                break
        f.close()
        if forceres==-1: raise IOError("Could not find force resolution")
    except IOError as e:
        print "WARNING:",e
        ictype,lx,nv = haloutils.get_zoom_params(hpath)
        forceres = 100./2.^lx/40.
    return forceres

def load_aqcat(whichAq,snap):
    assert whichAq in ['A','B','C','D','E','F']
    if snap > 127: 
        raise ValueError("Aquarius is snaps 0-127")
    rspath = global_basepath+'/aquarius/Aq-'+whichAq+'/2/halos'
    return RDR.RSDataReader(rspath,snap,version=7)
