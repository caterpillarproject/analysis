import numpy as np
import pylab as plt
from haloutils import get_numsnaps,get_foldername,find_halo_paths
from haloutils import check_last_subfind_exists,check_last_rockstar_exists
import os
import readsnapshots.readsnap as rs
import readhalos.readsubf as rsf
import readhalos.readgroup as rg
import readhalos.RSDataReader as rdr
import asciitable
from optparse import OptionParser
from operator import itemgetter

def create_summary_table(filename,
                         nrvirlist=[3,4,5,6],levellist=[11,12,13,14],ictype="BB",
                         verbose=False):
    if os.path.exists(filename):
        print filename+' exists; overwriting'
    f = open(filename,'w')
    ## These are the data written to the summary table
    f.write("#hid ictype lx nvir lastsnap mhires groupx groupy groupz fofnpart subfnum subfx subfy subfz subfmass subfvmax subfnpart rsid rsx rsy rsz rsmass rsvmax rsnpart\n")

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
        #groupcat = rg.group_tab(outpath+'/outputs',lastsnap)
        subcat = rsf.subfind_catalog(outpath+'/outputs',lastsnap)

        ## Pick best group by the one with the most hi-res particles
        #bestgroup = np.where(groupcat.GroupMassType[:,1] == np.max(groupcat.GroupMassType[:,1]))[0]
        bestgroup = min(enumerate(subcat.group_contamination_count[0:3]),key=itemgetter(1))[0]
        #if len(bestgroup) != 1:
        #    raise RuntimeError("wrong number of groups (should be 1): "+strlen(bestgroup))
        #bestgroup = bestgroup[0]
        #groupcm = groupcat.GroupCM[bestgroup,:]
        groupcm = subcat.group_pos[bestgroup,:]
        groupnpart = subcat.group_len[bestgroup]

        #hid ictype lx nvir lastsnap mhires groupx groupy groupz fofnpart
        groupstr = " ".join([folderparts[0],folderparts[1],folderparts[5][2:4],folderparts[7][2],str(lastsnap)])
        groupstr += " %3.2e %4.2f %4.2f %4.2f %i" % (mhires,groupcm[0],groupcm[1],groupcm[2],groupnpart)
        substart = subcat.group_firstsub[bestgroup]
        subnum = subcat.group_nsubs[bestgroup]
        
        subii = np.arange(substart,substart+subnum)
        submass = subcat.sub_mass[subii]*10**10/hubble
        subpos = subcat.sub_pos[subii]
        subvmax = subcat.sub_vmax[subii]
        subnpart= subcat.sub_len[subii]
        sub_id_list = np.where((submass <= 1e13) & (submass >= 5e11))[0]
        if len(sub_id_list)==0:
            print "ERROR: no halos satisfy selection criteria for "+foldername
        else:
            print foldername+" has "+str(len(sub_id_list))+" subhalos"

        if outpath in rspathlist:
            rscat = rdr.RSDataReader(outpath+'/halos',lastsnap)
            rsids = np.array(rscat['id'])
            rspos = np.array(rscat[['posX','posY','posZ']])
            rsmass= np.array(rscat['mvir'])/rscat.h0
            rsvmax= np.array(rscat['vmax'])
            rsnpart=np.array(rscat['npart'])

        for subid in sub_id_list:
            substr = "%i %4.2f %4.2f %4.2f %4.3e %f %i" % (subid,subpos[subid,0],subpos[subid,1],subpos[subid,2],submass[subid],subvmax[subid],subnpart[subid])
            #substr = " ".join([str(subid),str(subpos[subid,0]),
            #                   str(subpos[subid,1]),str(subpos[subid,2]),
            #                   str(submass[subid]),str(subvmax[subid]),str(subnpart[subid])])
            if outpath in rspathlist:
                dist = np.sqrt(np.sum((rspos - subpos[subid,:])**2,1))
                massratio = np.abs(rsmass/submass[subid]-1)
                rsselect = np.where( (dist < 0.01) & (massratio < .5))[0]
                if len(rsselect) == 1:
                    rsselect = rsselect[0]
                    rsstr = " ".join([str(rsids[rsselect]),str(round(rspos[rsselect][0],3)),
                                      str(round(rspos[rsselect][1],3)),str(round(rspos[rsselect][2],3)),
                                      str(rsmass[rsselect]),str(round(rsvmax[rsselect],3)),
                                      str(rsnpart[rsselect])])
                elif len(rsselect) ==0:
                    rsstr = "-999 "*7
                else:
                    rsstr = (str(-1*len(rsselect))+" ")*7
            else:
                rsstr = "-1 "*7
            f.write(groupstr+" "+substr+" "+rsstr+"\n")
    f.close()
    print "Summary written to "+filename

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("-s","--summarytable",
                      action="store_true",dest="summaryflag",default=False,
                      help="If set, recompute the halo summary table")
    parser.add_option("-f","--filename",
                      action="store",type="string",dest="filename",
                      default='/bigbang/data/AnnaGroup/caterpillar/halos/halosummary.txt',
                      help="File name for summary. Default /bigbang/data/AnnaGroup/caterpillar/halos/halosummary.txt")
                      
    parser.add_option("-H","--halo",
                      action="store",type="int",dest="haloid",default=-1,
                      help="Halo ID to plot (default is to plot all)")
    parser.add_option("-p","--plot",
                      action="store_true",dest="plotflag",default=False,
                      help="If set, plot halo properties")

    options,args = parser.parse_args()
    if options.summaryflag:
        create_summary_table(options.filename,nrvirlist=[3])
