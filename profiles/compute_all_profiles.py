import numpy as np
import profilefit
from densityprofile import densityprofile
from findhalos.haloutils import get_numsnaps,get_foldername,find_halo_paths
from findhalos.haloutils import check_last_rockstar_exists

import readhalos.RSDataReader as RDR
import readsnapshots.readsnap as rs

from sets import Set
['H21047', 'BB', 'Z127', 'P7', 'LN7', 'LX12', 'O4', 'NV3']

def get_best_halo_id(outpath,cat):
    fileparts = get_foldername(outpath).split('_')
    hid = int(fileparts[0][1:])
    lx = int(fileparts[5][2:4])
    nv = int(fileparts[7][2])

    if hid==241932:
        if lx == 11:
            return 5563
        if lx == 12:
            return 36854
    if hid==121869:
        if lx == 11:
            return 1221
        if lx == 12:
            return 12232
    return int(cat[cat['npart']==np.max(cat['npart'])]['id'])

def auto_rarr(rvir,dr=1):
    return np.arange(0,rvir+dr,dr)/1000. #kpc to Mpc

def compute_all_profiles(rarr=-1,ictype="BB",levellist=[11,12,13,14],nrvirlist=[3]):
    if rarr == -1:
        autorflag = True
    else:
        autorflag = False
    halopathlist = find_halo_paths(ictype=ictype,levellist=levellist,nrvirlist=nrvirlist)
    rspathlist = []
    for outpath in halopathlist:
        if check_last_rockstar_exists(outpath):
            #if get_foldername(outpath)[0:7] != 'H242183':# and \
                    #get_foldername(outpath)[0:7] != 'H140666': #this guy is messed up
            hname = get_foldername(outpath)[0:7]
            if hname == 'H268422' or hname == 'H241932' or \
                    hname[:-1] == 'H21047' or hname == 'H121869':
                rspathlist.append(outpath)
    for outpath in rspathlist:
        print "computing profile for "+get_foldername(outpath)
        lastsnap = get_numsnaps(outpath)-1
        snapstr = str(lastsnap).zfill(3)
        snapfile = outpath+'/outputs/snapdir_'+snapstr+'/snap_'+snapstr
        header = rs.snapshot_header(snapfile+'.0')
        snapIDs = rs.read_block(snapfile,"ID  ",parttype=1)
        argsorted = np.argsort(snapIDs)
        snapPOS = rs.read_block(snapfile,"POS ",parttype=1,doubleprec=False)
        cat = RDR.RSDataReader(outpath+'/halos',lastsnap,AllParticles=True)
        haloid = get_best_halo_id(outpath,cat)

        print "rsid %i mvir %3.2e" % (haloid,cat.ix[haloid]['mvir']/cat.h0)
        haloparts = cat.get_particles_from_halo(haloid)
        halopos = np.array(cat.ix[haloid][['posX','posY','posZ']])
        halorvir = float(cat.ix[haloid]['rvir']) #h^-1 kpc
        if autorflag:
            rarr = auto_rarr(halorvir)
        #print " ",haloid,len(haloparts),halopos,halorvir
        #print " ",rarr
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
        f = open(outpath+'/rs-halo-profile.dat','w')
        f.write(str(p03rmin)+"\n")
        for r,rho in zip(rarr[1:],rhoarr):
            f.write(str(r)+" "+str(rho)+"\n")
        f.close()
        print ""

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
    compute_all_profiles(levellist=[11,12],nrvirlist=[3])
