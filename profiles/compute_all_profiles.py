import numpy as np
from optparse import OptionParser
import profilefit
from densityprofile import densityprofile,getr200
from findhalos.haloutils import get_numsnaps,get_foldername,find_halo_paths
from findhalos.haloutils import check_last_rockstar_exists,check_last_subfind_exists

import readhalos.RSDataReader as RDR
import readhalos.readsubf as readsubf
import readsnapshots.readids as readids
#import readsnapshots.readsnap as rs
#import readsnapshots.readidsHDF5 as readids
import readsnapshots.readsnapHDF5 as rs

from sets import Set
from operator import itemgetter

def get_best_halo_id(outpath,cat):
    fileparts = get_foldername(outpath).split('_')
    hid = int(fileparts[0][1:])
    lx = int(fileparts[5][2:4])
    nv = int(fileparts[7][2])
    
    print "  get_best_halo_id: guessing based on cat['npart']"
    return int(cat[cat['npart']==np.max(cat['npart'])]['id'])

def auto_rarr(rvir,dr=1):
    #return np.arange(0,(1.5*rvir)+dr,dr)/1000. #kpc to Mpc
    return np.logspace(-5,0,100)

#def compute_all_profiles(subfind=False,subfindradius=False,rockstar_get_all=False,
#                         rarr=-1,ictype="BB",
#                         levellist=[11,12,13,14],nrvirlist=[3]):
#    if sum([subfind,subfindradius,rockstar_get_all]) > 1:
#        print "ERROR: more than one thing is true"
#        return
#    #haloidlist = [268422]#, 21047, 241932, 121869]
#    #haloidlist = [1194083,1292049,1476079,230667]
#    #haloidlist = [4847,649524,1327707]
#
#    ## Sheet 1
#    haloidlist = [1194083,1232333,1292049,1725139,230667,4847,649524,706754,1327707,1476079]
#    
#    haloidlist = ['H'+str(hid) for hid in haloidlist]
#    if rarr == -1:
#        autorflag = True
#    else:
#        autorflag = False
#    halopathlist = find_halo_paths(ictype=ictype,levellist=levellist,nrvirlist=nrvirlist,onlychecklastsnap=True,verbose=False)
#    availablepathlist = []
#    if subfind or subfindradius:
#        checkfn = check_last_subfind_exists
#    else:
#        checkfn = check_last_rockstar_exists
#    for outpath in halopathlist:
#        if checkfn(outpath):
#            hname = (get_foldername(outpath).split('_'))[0]
#            if hname in haloidlist:
#                print get_foldername(outpath)
#                availablepathlist.append(outpath)
#    for outpath in availablepathlist:
#        compute_one_profile(outpath,rarr,subfind=subfind,subfindradius=subfindradius,rockstar_get_all=rockstar_get_all)

def compute_one_profile(outpath,rarr=-1,
                        subfind=False,subfindradius=False,rockstar_get_all=False):
    print "----------------------------------"
    print "computing profile for "+get_foldername(outpath)
    lastsnap = get_numsnaps(outpath)-1
    snapstr = str(lastsnap).zfill(3)
    snapfile = outpath+'/outputs/snapdir_'+snapstr+'/snap_'+snapstr
    header = rs.snapshot_header(snapfile+'.0')
    if subfindradius:
        snapIDs = rs.read_block(snapfile,"ID  ",parttype=-1)
        snapPOS = rs.read_block(snapfile,"POS ",parttype=-1)#,doubleprec=False)
    else:
        snapIDs = rs.read_block(snapfile,"ID  ",parttype=1)
        snapPOS = rs.read_block(snapfile,"POS ",parttype=1)#,doubleprec=False)
    argsorted = np.argsort(snapIDs)
    if subfind or subfindradius:
        cat = readsubf.subfind_catalog(outpath+'/outputs',lastsnap)
        bestgroup = min(enumerate(cat.group_contamination_count[0:5]),key=itemgetter(1))[0]
        if bestgroup != 0: 
            print "WARNING: bestgroup != 0 (may have inefficient computation)"
            print "group masses:",cat.group_mass[0:(bestgroup+1)]*10**10/header.hubble
            print "group npart:",cat.group_len[0:(bestgroup+1)]
            print "group contamination:",cat.group_contamination_count[0:(bestgroup+1)]
        substart = cat.group_firstsub[bestgroup]
        subnum = cat.group_nsubs[bestgroup]
        print "R200c mass of bestgroup %i: %3.2e" % (bestgroup, cat.group_m_crit200[bestgroup] * 10**10/header.hubble)
        
        if subfind:
            halopos = cat.group_pos[bestgroup]
            halorvir= cat.group_r_crit200[bestgroup]
            subii = np.arange(substart,substart+subnum)
            submass = cat.sub_mass[subii]*10**10/header.hubble
            subpos = cat.sub_pos[subii]; subdist = np.sum((subpos-halopos)**2,axis=1)
            subvmax = cat.sub_vmax[subii]
            subnpart= cat.sub_len[subii]
            #print submass; print subdist; print subdist < halorvir**2
            sub_id_list = np.where((submass <= 1e13) & (submass >= 5e11) & (subdist < halorvir**2))[0]
            if len(sub_id_list) == 1:
                subid = sub_id_list[0]+substart
            else:
                print "%i candidates in relevant mass range, skipping..." % (len(sub_id_list))
                print submass; print subdist
                print submass[sub_id_list]
                print "Closest subhalo: ",np.min(np.sqrt(subdist))
                #print subnpart[sub_id_list]
                return
            haloparts = readids.subid_file(outpath+'/outputs',lastsnap,cat.group_offset[bestgroup],cat.sub_len[cat.group_firstsub[bestgroup]])
            #haloparts = readids.subid_file(outpath+'/outputs',lastsnap)
            haloparts = haloparts.SubIDs
            halopos = cat.sub_pos[subid,:]
            halorvir = cat.group_r_crit200[bestgroup]*1000.
            halomass = cat.sub_mass[subid]*10**10/header.hubble
            print "bestgroup %i subid %i mvir %3.2e" % (bestgroup,subid,halomass)
            print " r200c %3.2f r200m %3.2f" % (cat.group_r_crit200[bestgroup]*1000.,cat.group_r_mean200[bestgroup]*1000.)
        elif subfindradius:
            halopos = cat.group_pos[bestgroup]
            halorvir = cat.group_r_crit200[bestgroup]*1000.
            halomass = cat.group_mass[bestgroup]*10**10/header.hubble
            dr = np.sqrt(np.sum((snapPOS-halopos)**2,1))
            haloparts = snapIDs[dr < (3*halorvir)]
            #groupstart = cat.group_offset[bestgroup]; groupend = groupstart+cat.group_len[bestgroup]
            #haloparts = snapIDs[groupstart:groupend]
    else:
        cat = RDR.RSDataReader(outpath+'/halos',lastsnap,AllParticles=True,version=6)
        haloid = get_best_halo_id(outpath,cat)
        print "rsid %i mvir %3.2e" % (haloid,cat.ix[haloid]['mvir']/cat.h0)
        if rockstar_get_all:
            haloparts = cat.get_all_particles_from_halo(haloid)
        else:
            haloparts = cat.get_particles_from_halo(haloid)
        halopos = np.array(cat.ix[haloid][['posX','posY','posZ']])
        halorvir = float(cat.ix[haloid]['rvir']) #h^-1 kpc; rvir is 18pi^2
        halomass = cat.ix[haloid]['mvir']/header.hubble
        
    if rarr==-1:
        rarr = auto_rarr(halorvir)
    try:
        rhoarr,p03rmin = densityprofile(rarr,snapPOS,argsorted,header,haloparts,halopos,power03=True)
        r200c = getr200(haloparts,snapPOS,argsorted,header,halopos)
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
        r200c = getr200(goodparts,snapPOS,argsorted,header,halopos)
    if subfind:
        f = open(outpath+'/subf-halo-profile.dat','w')
    elif subfindradius:
        f = open(outpath+'/subf-halo-profile-radius.dat','w')
    elif rockstar_get_all:
        f = open(outpath+'/rs-halo-profile-allpart.dat','w')
    else:
        f = open(outpath+'/rs-halo-profile.dat','w')
    f.write(str(p03rmin)+" "+str(halorvir)+" "+str(r200c)+" "+str(halomass)+"\n")
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
    # Use the antares submission script to run this automatically
    parser = OptionParser()
    options,args = parser.parse_args()
    numargs = len(args); numtypes = numargs-1
    outpath = args[0]
    print outpath
    if numtypes==0:
        print "No arguments passed to run things"
    for whichtype in args[1:]:
        if whichtype=='0':   #rockstar no subhalos
            print "---Rockstar no subhalos---"
            compute_one_profile(outpath)
        elif whichtype=='1': #rockstar yes subhalos
            print "---Rockstar with subhalos---"
            compute_one_profile(outpath,rockstar_get_all=True)
        elif whichtype=='2': #subfind particles
            print "---Subfind---"
            compute_one_profile(outpath,subfind=True)
        elif whichtype=='3': #subfind all particles within R200
            print "---Subfind R200---"
            compute_one_profile(outpath,subfindradius=True)
        else:
            print "ERROR: invalid profile type:",whichtype
    #compute_all_profiles(levellist=levellist,nrvirlist=nrvirlist)
    #compute_all_profiles(levellist=levellist,nrvirlist=nrvirlist,rockstar_get_all=True)
    #compute_all_profiles(levellist=levellist,nrvirlist=nrvirlist,subfind=True)
    #compute_all_profiles(levellist=levellist,nrvirlist=nrvirlist,subfindradius=True)
