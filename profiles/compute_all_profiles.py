import numpy as np
import profilefit
from densityprofile import densityprofile
from findhalos.haloutils import get_numsnaps,get_foldername,find_halo_paths
from findhalos.haloutils import check_last_rockstar_exists

import readhalos.RSDataReader as RDR
import readhalos.readsubf as readsubf
import readsnapshots.readids as readids
import readsnapshots.readsnap as rs

from sets import Set
from operator import itemgetter

def get_best_halo_id(outpath,cat):
    fileparts = get_foldername(outpath).split('_')
    hid = int(fileparts[0][1:])
    lx = int(fileparts[5][2:4])
    nv = int(fileparts[7][2])
    
    if hid==241932:
        if lx == 11 and nv == 3:
            return 5563
        if lx == 12 and nv == 3:
            return 36854
    if hid==121869:
        if lx == 11 and nv == 3:
            return 1221
        if lx == 12 and nv == 3:
            return 12232
    if hid==268422:
        if lx == 11 and nv == 3:
            return 1911
        if lx == 12 and nv == 3:
            return 24939
    if hid==21047:
        if lx == 11 and nv == 3:
            return 2338
        if lx == 12 and nv == 3:
            return 17844
    #default heuristic:
    print "  get_best_halo_id: guessing based on cat['npart']"
    return int(cat[cat['npart']==np.max(cat['npart'])]['id'])

def auto_rarr(rvir,dr=1):
    #return np.arange(0,(1.5*rvir)+dr,dr)/1000. #kpc to Mpc
    return np.logspace(-5,0,100)

def compute_all_profiles(subfind=False,rockstar_get_all=False,rarr=-1,ictype="BB",levellist=[11,12,13,14],nrvirlist=[3]):
    if rarr == -1:
        autorflag = True
    else:
        autorflag = False
    halopathlist = find_halo_paths(ictype=ictype,levellist=levellist,nrvirlist=nrvirlist)
    rspathlist = []
    for outpath in halopathlist:
        if check_last_rockstar_exists(outpath):
            hname = get_foldername(outpath)[0:7]
            if hname == 'H268422' or hname == 'H241932' or \
                    hname[:-1] == 'H21047' or hname == 'H121869':
                rspathlist.append(outpath)
    for outpath in rspathlist:
        print "----------------------------------"
        print "computing profile for "+get_foldername(outpath)
        lastsnap = get_numsnaps(outpath)-1
        snapstr = str(lastsnap).zfill(3)
        snapfile = outpath+'/outputs/snapdir_'+snapstr+'/snap_'+snapstr
        header = rs.snapshot_header(snapfile+'.0')
        snapIDs = rs.read_block(snapfile,"ID  ",parttype=1)
        argsorted = np.argsort(snapIDs)
        snapPOS = rs.read_block(snapfile,"POS ",parttype=1,doubleprec=False)
        if subfind:
            cat = readsubf.subfind_catalog(outpath+'/outputs',lastsnap)
            bestgroup = min(enumerate(cat.group_contamination_count[0:10]),key=itemgetter(1))[0]
            substart = cat.group_firstsub[bestgroup]
            subnum = cat.group_nsubs[bestgroup]
            
            subii = np.arange(substart,substart+subnum)
            submass = cat.sub_mass[subii]*10**10/header.hubble
            subpos = cat.sub_pos[subii]
            subvmax = cat.sub_vmax[subii]
            subnpart= cat.sub_len[subii]
            sub_id_list = np.where((submass <= 1e13) & (submass >= 5e11))[0]
            if len(sub_id_list) == 1:
                subid = sub_id_list[0]+substart
            else:
                print "More than one candidate in relevant mass range, skipping..."
                continue
            haloparts = readids.subid_file(outpath+'/outputs',lastsnap,cat.group_offset[bestgroup],cat.sub_len[cat.group_firstsub[bestgroup]])
            haloparts = haloparts.SubIDs
            halopos = cat.sub_pos[subid,:]
            halorvir = cat.group_r_crit200[bestgroup]*1000.
            print "bestgroup %i subid %i mvir %3.2e" % (bestgroup,subid,cat.sub_mass[subid]*10**10/header.hubble)
        else:
            cat = RDR.RSDataReader(outpath+'/halos',lastsnap,AllParticles=True)
            haloid = get_best_halo_id(outpath,cat)
            print "rsid %i mvir %3.2e" % (haloid,cat.ix[haloid]['mvir']/cat.h0)
            if rockstar_get_all:
                haloparts = cat.get_all_particles_from_halo(haloid)
            else:
                haloparts = cat.get_particles_from_halo(haloid)
            halopos = np.array(cat.ix[haloid][['posX','posY','posZ']])
            halorvir = float(cat.ix[haloid]['rvir']) #h^-1 kpc; rvir is 18pi^2
        if autorflag:
            rarr = auto_rarr(halorvir)
        try:
            rhoarr,p03rmin = densityprofile(rarr,snapPOS,argsorted,header,haloparts,halopos,power03=True)
        except IndexError as e:
            print "  ",e
            print "  ---contamination in "+get_foldername(outpath)
            print "  removing those particles to compute density profile..."
            nall = len(haloparts)
            nbad,goodparts = remove_contaminated_particles(snapfile,haloparts, snapIDs,count_all=True)
            print "  "+str(nbad)+" contamination particles out of "+str(nall)
            partmasstab = rs.snapshot_header(snapfile+'.0').massarr
            print "  "+str((np.sum(nbad[2:] * partmasstab[2:]))/(np.sum(nbad * partmasstab)))
            rhoarr,p03rmin = densityprofile(rarr,snapPOS,argsorted,header,goodparts,halopos,power03=True)
        if subfind:
            f = open(outpath+'/subf-halo-profile.dat','w')
        elif rockstar_get_all:
            f = open(outpath+'/rs-halo-profile-allpart.dat','w')
        else:
            f = open(outpath+'/rs-halo-profile.dat','w')
        f.write(str(p03rmin)+" "+str(halorvir)+"\n")
        for r,rho in zip(rarr[2:],rhoarr[1:]):
            f.write(str(r)+" "+str(rho)+"\n")
        f.close()

def remove_contaminated_particles(snapfile,haloparts, snapIDs,count_all = False):
    hpartset = Set(haloparts); snapset = Set(snapIDs)
    goodparts = list(hpartset & snapset)
    if count_all:
        badparts = hpartset - snapset
        npart = [0,len(goodparts),0,0,0,0]
        for parttype in [2,3,4,5]:
            contamIDs = Set(rs.read_block(snapfile,"ID  ",parttype=parttype))
            theseparts = badparts & contamIDs
            npart[parttype] = len(theseparts)
            badparts = badparts - theseparts
        if len(badparts) != 0:
            npart[0] = len(badparts)
            print "    remove_contaminated_particles: ERROR not all IDs are part 1-5"
        return npart,np.array(goodparts)
    else:
        nbad = len(haloparts)-len(goodparts)
        return nbad,np.array(goodparts)

if __name__=="__main__":
    levellist = [11,12]; nrvirlist = [3]
    #compute_all_profiles(levellist=levellist,nrvirlist=nrvirlist)
    compute_all_profiles(levellist=levellist,nrvirlist=nrvirlist,subfind=True)
    #compute_all_profiles(levellist=levellist,nrvirlist=nrvirlist,rockstar_get_all=True)
