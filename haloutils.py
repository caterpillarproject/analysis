import os,sys,platform
if 'compute-0-' in platform.node():
    import matplotlib
    matplotlib.use('Agg')

import numpy as np
import asciitable
import pickle
import pandas as pd

import readsnapshots.readsnapHDF5_greg as rsg
import readhalos.RSDataReader as RDR
import readhalos.readsubf as RSF
import mergertrees.MTCatalogue as MTC
import glob



def determinebasepath(node):
    if node == "csr-dyn-150.mit.edu":
        basepath = '/Users/griffen/Desktop/'
    elif node == "Brendans-MacBook-Pro.local":
        basepath = '/Users/griffen/Desktop/'
    elif node == "spacebase":
        basepath = '/bigbang/data/AnnaGroup/'
    elif node == "bigbang.mit.edu":
        basepath = '/bigbang/data/AnnaGroup/'
    elif node == "antares":
        basepath = '/bigbang/data/AnnaGroup/'
    elif 'compute-0-' in node:
        basepath = '/bigbang/data/AnnaGroup/'
    else:
        raise ValueError(node+" is not a valid node")
        
    return basepath

global_basepath = os.path.normpath(determinebasepath(platform.node()))
global_halobase = global_basepath+'/caterpillar/halos'
global_prntbase = global_basepath+'/caterpillar/parent/gL100X10'


cid2hid = {1:1631506,
           2:264569,
           3:1725139,
           4:447649,
           5:5320,
           6:581141,
           7:94687,
           8:1130025,
           9:1387186,
           10:581180,
           11:1725372,
           12:1354437,
           13:1725272,
           14:1195448,
           15:1599988,
           16:796175,
           17:388476,
           18:1079897,
           19:94638,
           20:95289,
           21:1232164,
           22:1422331,
           23:196589,
           24:1268839}

hid2name = {}
hid2catnum = {}
for k,v in cid2hid.items():
    hid2name[v] = 'Cat-'+str(k)
    hid2catnum[v] = k

def hid_name(hid):
    return hid2name[hidint(hid)]
def hpath_name(hpath):
    return hid_name(get_parent_hid(hpath))
def hid_catnum(hid):
    return hid2catnum[hidint(hid)]
def hpath_catnum(hpath):
    return hid_catnum(get_parent_hid(hpath))

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
    if os.path.exists(outpath+'/ExpansionList'):
        return sum(1 for line in open(outpath+'/ExpansionList'))
    else:
        print "WARNING: "+outpath+"/ExpansionList not found, using default (256)"
        return 256
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

def get_paper_paths_lx(do_lx):
    return [get_hpath_lx(hid,do_lx) for hid in hid2name.keys()]

def get_paper_paths():
    return [global_halobase+"/H"+str(hid) for hid in hid2name.keys()]


def get_available_hpaths(hid,contam=False,
                         checkgadget=True,
                         onlychecklastsnap=True,
                         checkallblocks=False,
                         hdf5=True,verbose=False,
                         basepath=global_halobase):
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
                                                   checkallblocks=checkallblocks,hdf5=hdf5,verbose=verbose):
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

def check_mergertree_exists(outpath,autoconvert=False,boundbin=True):
    if boundbin: halodir = 'halos_bound'
    else: halodir = 'halos'
    ascii_exists = os.path.exists(outpath+'/'+halodir+'/trees/tree_0_0_0.dat')
    binary_exists = os.path.exists(outpath+'/'+halodir+'/trees/tree.bin')
    if autoconvert and ascii_exists and not binary_exists:
        print "---check_mergertree_exists: Automatically converting ascii to binary"
        MTC.convertmt(outpath+'/'+halodir+'/trees',version=4)
        binary_exists = os.path.exists(outpath+'/'+halodir+'/trees/tree.bin')
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
                    hdf5=True,verbose=False):
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

def find_halo_paths(basepath=global_halobase,
                    nrvirlist=[3,4,5,6],levellist=[11,12,13,14],
                    ictypelist=["BA","BB","BC","BD","EA","EB","EC","EX","CA","CB","CC"],
                    contamsuite=False,
                    require_rockstar=False,require_subfind=False,
                    require_mergertree=False,autoconvert_mergertree=False,
                    require_sorted=False,
                    checkallblocks=False,
                    onlychecklastsnap=False,verbose=False,hdf5=True,
                    use_fullbin_rockstar=False):
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
        hpathlist = get_available_hpaths(haloid,contam=contamsuite,checkgadget=False,
                                         basepath=basepath)
        for hpath in hpathlist:
            ictype,levelmax,nrvir = get_zoom_params(hpath)
            if (int(levelmax) in levellist and int(nrvir) in nrvirlist and ictype in ictypelist):
                try:
                    if gadget_finished(hpath,
                                       onlychecklastsnap=onlychecklastsnap,
                                       checkallblocks=checkallblocks,
                                       hdf5=hdf5,verbose=verbose):
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
def load_zoomid(hpath,filename=global_halobase+"/parent_zoom_index.txt"):
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
    
def load_haloprops(hpath,filename=global_halobase+"/parent_zoom_index.txt"):
    row = _load_index_row(hpath,filename=filename)
    mvir = float(row['mvir']) # physical Msun
    rvir = float(row['rvir']) # physical kpc
    vvir = np.sqrt(4.34e-6 * mvir/rvir) # physical km/s
    return mvir,rvir,vvir

def load_pcatz0(old=False):
    if old:
        return RDR.RSDataReader(global_basepath+"/caterpillar/parent/RockstarData",63,version=2)
    else:
        return RDR.RSDataReader(global_prntbase+"/rockstar",127,version=6)

def load_scat(hpath):
    if "LX14" in hpath:
        return RSF.subfind_catalog(hpath+'/outputs',255,double=True)
    else: 
        return RSF.subfind_catalog(hpath+'/outputs',255)

def load_rscat(hpath,snap,verbose=True,halodir='halos_bound',unboundfrac=None,minboundpart=None,version=None,rmaxcut=True):
    if version != None:
        rscat = RDR.RSDataReader(hpath+'/'+halodir,snap,version=version,digits=1,unboundfrac=unboundfrac,minboundpart=minboundpart)
    else:
        try:
            rscat = RDR.RSDataReader(hpath+'/'+halodir,snap,version=10,digits=1,unboundfrac=unboundfrac,minboundpart=minboundpart)
        except IOError as e: #try to identify a unique valid rockstar version
            print e
            versionlist = [2,3,4,5,6,7,8,9]
            testlist = []
            for version in versionlist:
                try:
                    rscat = RDR.RSDataReader(hpath+'/'+halodir,snap,version=version,digits=1,unboundfrac=unboundfrac,minboundpart=minboundpart)
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
                    rscat = RDR.RSDataReader(hpath+'/'+halodir,snap,version=version,digits=1,unboundfrac=unboundfrac,minboundpart=minboundpart)

    if rmaxcut:
        zoomid = load_zoomid(hpath)
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

def load_mtc(hpath,verbose=True,halodir='halos_bound',treedir='trees',**kwargs):
    return MTC.MTCatalogue(hpath+'/'+halodir+'/'+treedir,version=4,**kwargs)
def load_zoom_mtc(hpath,verbose=True,halodir='halos_bound',treedir='trees',**kwargs):
    return MTC.MTCatalogue(hpath+'/'+halodir+'/'+treedir,version=4,haloids=[load_zoomid(hpath)],**kwargs)
    
def load_pmtc(hpath,verbose=True,halodir='rockstar',treedir='trees',**kwargs):
    return MTC.MTCatalogue(hpath+'/'+halodir+'/'+treedir,version=4,**kwargs)

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
            if s[0]=="SofteningHaloMaxPhys":
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

def get_halo_header(hpath,snap=255):
    return rsg.snapshot_header(hpath+"/outputs/snapdir_"+str(snap)+"/snap_"+str(snap))

def get_colors(ncolors=12):
    colors = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
              (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
              (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
              (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
              (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),
              (32, 159, 117), (43, 166, 241),  (6, 115,   1), (203,  56,  62)] 
    assert len(colors)==len(hid2name)
    
    for i in range(ncolors):  
        r, g, b = colors[i]  
        colors[i] = (r / 255., g / 255., b / 255.)  

    return colors

def get_colors_for_halos(nhalos=len(hid2name)):
    colors = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
              (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
              (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
              (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
              (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),
              (32, 159, 117), (43, 166, 241),  (6, 115,   1), (203,  56,  62)] 
    assert len(colors)==len(hid2name)
    
    for i in range(nhalos):  
        r, g, b = colors[i]  
        colors[i] = (r / 255., g / 255., b / 255.)  

    index = {}
    for i in range(nhalos):
        index[hid2name.keys()[i]] = colors[i]

    return index

