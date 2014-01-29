import numpy as np
import pylab as plt
from haloutils import get_numsnaps,get_foldername,find_halo_paths
import os
import readsnapshots.readsnap as rs
import readhalos.readsubf as rsf
import readhalos.readgroup as rg
import readhalos.RSDataReader as rdr

def check_last_subfind_exists(outpath):
    numsnaps = get_numsnaps(outpath)
    lastsnap = numsnaps - 1; snapstr = str(lastsnap).zfill(3)
    group_tab = os.path.exists(outpath+'/outputs/groups_'+snapstr+'/group_tab_'+snapstr+'.0')
    subhalo_tab = os.path.exists(outpath+'/outputs/groups_'+snapstr+'/subhalo_tab_'+snapstr+'.0')
    return group_tab and subhalo_tab

def check_last_rockstar_exists(outpath):
    numsnaps = get_numsnaps(outpath)
    lastsnap = numsnaps - 1; snapstr = str(lastsnap)
    halo_exists = os.path.exists(outpath+'/halos/halos_'+snapstr+'/halos_'+snapstr+'.0.bin')
    part_exists = os.path.exists(outpath+'/halos/halos_'+snapstr+'/halos_'+snapstr+'.0.particles')
    return halo_exists and part_exists

def create_summary_table(filename='/bigbang/data/AnnaGroup/caterpillar/halos/halosummary.txt',
                         nrvirlist=[3,4,5,6],levellist=[11,12,13,14],ictype="BB",
                         verbose=False):
    if os.path.exists(filename):
        print filename+' exists; overwriting'
    f = open(filename,'w')
    ## These are the data written to the summary table
    f.write("#hid ictype lx nvir lastsnap mhires groupx groupy groupz fofnpart subnum subx suby subz submass subvmax subnpart rsid\n")

    halopathlist = find_halo_paths(nrvirlist=nrvirlist,levellist=levellist,ictype=ictype,onlychecklastsnap=True,verbose=False)
    subfpathlist = []
    rspathlist = []
    for outpath in halopathlist:
        if check_last_subfind_exists(outpath):
            subfpathlist.append(outpath)
        elif verbose:
            print get_foldername(outpath)+' does not have subfind run'
        if check_last_rockstar_exists(outpath):
            rspathlist.append(outpath)
        elif verbose:
            print get_foldername(outpath)+' does not have rockstar run'
    if verbose:
        print "Total not run:",len(halopathlist)-len(subfpathlist)
        print "Number of subf paths:",len(subfpathlist)
        print "Number of rs paths:",len(rspathlist)

    for outpath in subfpathlist:
        foldername = get_foldername(outpath)
        folderparts = foldername.split('_')
        numsnaps = get_numsnaps(outpath); lastsnap = numsnaps-1
        lastsnapstr = str(lastsnap).zfill(3)
        header = rs.snapshot_header(outpath+'/outputs/snapdir_'+lastsnapstr+'/snap_'+lastsnapstr+
'.0')
        hubble = header.hubble
        mhires = header.massarr[1] * 10.**10/hubble
        groupcat = rg.group_tab(outpath+'/outputs',lastsnap)
        subcat = rsf.subfind_catalog(outpath+'/outputs',lastsnap)

        ## Pick best group by the one with the most hi-res particles
        bestgroup = np.where(groupcat.GroupMassType[:,1] == np.max(groupcat.GroupMassType[:,1]))[0]
        if len(bestgroup) != 1:
            raise RuntimeError("wrong number of groups (should be 1): "+strlen(bestgroup))
        bestgroup = bestgroup[0]
        groupcm = groupcat.GroupCM[bestgroup,:]
        groupnpart = str(subcat.group_len[bestgroup])

        #hid ictype lx nvir lastsnap mhires groupx groupy groupz fofnpart
        groupstr = " ".join([folderparts[0],folderparts[1],folderparts[5][2:4],folderparts[7][2],str(lastsnap),str(mhires),str(groupcm[0]),str(groupcm[1]),str(groupcm[2]),groupnpart])
        substart = subcat.group_firstsub[bestgroup]
        subnum = subcat.group_nsubs[bestgroup]
        
        ## OLD
        #subcm = subcat.sub_cm[substart,:]
        #subid = substart
        #submass = subcat.sub_mass[substart]*10.**10/hubble
        #subvmax = subcat.sub_vmax[substart]
        #subnpart = subcat.sub_len[substart]
        #rsid = "NA"
        #substr = " ".join([str(subid),str(subcm[0]),str(subcm[1]),str(subcm[2]),str(submass),str(subvmax),str(subnpart),str(rsid)])

        ## The most massive subhalo is probably not the halo we want, so pick all halos in the mass range
        #subii = np.arange(substart,substart+subnum)
        submass = subcat.sub_mass*10**10/hubble
        subpos = subcat.sub_pos; R = np.sqrt(np.sum((subpos - groupcm)**2,1))
        subvmax = subcat.sub_vmax
        subnpart= subcat.sub_len
        sub_id_list = np.where((submass <= 9e12) & (submass >= 7e11) & (R < 4.0))[0] 
        if len(sub_id_list)==0:
            print "ERROR: no halos satisfy selection criteria for "+foldername
        else:
            print foldername+" has "+str(len(sub_id_list))+" subhalos"

        if outpath in rspathlist:
            rscat = rdr.RSDataReader(outpath+'/halos',lastsnap)
            rsids = np.array(rscat['id'])
            rsmass= np.array(rscat['mvir'])
            rspos = np.array(rscat[['posX','posY','posZ']])
            ## TODO 

        for subid in sub_id_list:
            if outpath in rspathlist:
                dist = np.sqrt(np.sum((rspos - subpos[subid,:])**2,1))
                #print "  mindist:",np.min(dist)
                massratio = np.abs(rsmass/submass[subid]-1)
                rsselect = np.where( (dist < 1) & (massratio < .8) )[0]
                if len(rsselect) == 1:
                    rsid = str(rsids[rsselect])
                elif len(rsselect) ==0:
                    rsid = "NONE!"
                else:
                    rsid = "FOUND"+str(len(rsselect))#+'_'+str(rsids[rsselect])
            else:
                rsid = "NA"
            substr = " ".join([str(subid),str(subpos[subid,0]),str(subpos[subid,1]),str(subpos[subid,2]),
                               str(submass[subid]),str(subvmax[subid]),str(subnpart[subid]),str(rsid)])
            f.write(groupstr+" "+substr+"\n")
    f.close()
    print "Summary written to "+filename

if __name__=="__main__":
    create_summary_table(nrvirlist=[3])
