import numpy as np
import pylab as plt
import haloutils
import asciitable
import os
from findhalos.caterpillar_findhalos import get_zoom_id
import seaborn as sns

def find_contam_dist():
    hpathlist = haloutils.find_halo_paths(contamsuite=True,onlychecklastsnap=True,
                                          require_rockstar=True,require_subfind=True)
    pcat = haloutils.load_pcatz0()

    numpaths = len(hpathlist)
    
    data = np.zeros((numpaths,),dtype={'names':['haloid','contamtype','icsize','min2','min3','min4','min5','mvir','vmax'],
                                       'formats':['i4','a10','i8','f4','f4','f4','f4','f8','f8']})
    for i,hpath in enumerate(hpathlist):
        icsize = os.path.getsize(hpath+'/ics.0')
        data[i]['icsize'] = icsize
        rscat = haloutils.load_rscat(hpath,255)
        scat = haloutils.load_scat(hpath)
        hid = haloutils.get_parent_hid(hpath)
        contamtype = haloutils.get_foldername(hpath).split('_')[-1]
        data[i]['haloid'] = hid; data[i]['contamtype']=contamtype
        print 'H'+str(hid)+'_'+contamtype
        zoomid = get_zoom_id(hid,rscat,scat,pcat)
        
        hpos = np.array(rscat.ix[zoomid][['posX','posY','posZ']])
        for parttype in [2,3,4,5]:
            ppos = haloutils.load_partblock(hpath,255,"POS ",parttype=parttype)
            dr = np.sqrt(np.sum((ppos-hpos)**2,1))
            data[i]['min'+str(parttype)] = np.min(dr)
        data[i]['mvir'] = rscat.ix[zoomid]['mvir']
        data[i]['vmax'] = rscat.ix[zoomid]['vmax']
    asciitable.write(data,'contamdist.dat')
    
def plot_contam_dist(whichfolder,haloidlist=None):
    labelmap = {'A':1,'B':2,'C':3,'D':4,'ELLIPSOID':0,'CONVEX':-1}
    labelmap = {'A4':1,'A5':2,'B':3,'C':4,'D':5,'ELLIPSOID4':-1,'ELLIPSOID5':0,'CONVEX4':-3,'CONVEX5':-2}
    labelmap = {'A4':1,'A5':2,'B':3,'C':4,'D':5,'ELLIPSOID4':-1,'ELLIPSOID5':0,'CONVEX4':-3,'CONVEX5':-2,
                'ELLIPSOIDB':6,'ELLIPSOIDC':7}
    labelmap = {'A4':3,'A5':4,'B':5,'C':6,'D':7,'ELLIPSOID4':-1,'ELLIPSOID5':0,'CONVEX4':-3,'CONVEX5':-2,
                'ELLIPSOIDB':1,'ELLIPSOIDC':2}
    labelmap = {'CA4':-2,'CA5':-1,
                'EA4':0,'EA5':1,'EB':2,'EC':3,
                'BA4':4,'BA5':5,'BB':6,'BC':7,'BD':8}


    #data = asciitable.read('contamdist.dat')
    #hidlabel = 'haloid'
    hidlabel = 'parentid'
    whichmiddle = -1
    if whichfolder=='low':
        data = asciitable.read('/bigbang/data/AnnaGroup/caterpillar/halos/low_mass_halos/contam_zoom_index.txt',Reader=asciitable.FixedWidth)
    elif 'middle' in whichfolder:
        if 'middle' != whichfolder:
            whichmiddle=int(whichfolder[-1])
            assert whichmiddle >=1 and whichmiddle <= 8
        data = asciitable.read('/bigbang/data/AnnaGroup/caterpillar/halos/middle_mass_halos/contam_zoom_index.txt',Reader=asciitable.FixedWidth)
    elif whichfolder=='high':
        data = asciitable.read('/bigbang/data/AnnaGroup/caterpillar/halos/high_mass_halos/contam_zoom_index.txt',Reader=asciitable.FixedWidth)
    print whichfolder
    hidlist = np.unique(data[hidlabel])

    if haloidlist != None:
        #print hidlist
        #print haloidlist
        hidlist = list(set(haloidlist).intersection(hidlist))

    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    for hid in hidlist:
        thisdat = data[data[hidlabel]==hid]
        contamtypes = thisdat['ictype']
        nv = thisdat['NV']
        tmplist = []
        for i in range(len(contamtypes)):
            if contamtypes[i] in ['BA','EA','CA']:
                tmplist.append(contamtypes[i]+str(nv[i]))
            else: tmplist.append(contamtypes[i])
        contamtypes = tmplist
        #print contamtypes

        contamx = np.array([labelmap[ct] for ct in contamtypes])
        iisort = np.argsort(contamx)
        contamx = contamx[iisort]
        mindist = np.array(thisdat['min2'])[iisort]
        ax1.plot(contamx,mindist,'o-',label=haloutils.hidstr(hid))
        filesizes = np.array(thisdat['icsize'])[iisort]
        ax2.plot(contamx,filesizes,'o-',label=haloutils.hidstr(hid))
    minx = np.min(labelmap.values()); maxx = np.max(labelmap.values())

    ax1.plot([minx-0.5,maxx+0.5],[1.0,1.0],'k:')
    ax1.set_ylabel('distance to parttype2 [Mpc/h]')
    ax2.set_ylabel('size of ics in MB')
    if whichmiddle != -1:
        ax1.set_ylim([0,1.8])
        #ax2.set_ylim([100,700])
    for ax in [ax1,ax2]:
        ax.set_xlim([minx-0.5,maxx+0.5])
        ax.set_xticks(np.sort(labelmap.values()))
        ax.set_xticklabels(['CV4','CV5','EL4','EL5','ELB','ELC','A4','A5','B','C','D'])
        ax.legend(loc='best')
    return fig1,fig2

if __name__=="__main__":
    mididlists = [[1079897,1129843,1130025,1232127,1232164,1268839],
                  [1292085,1327666,1354437,1386703,1387186,1422331],
                  [1599902,1599988,1631506,1666505,1697496], #1452004,
                  [1725139,1725272,1725372,1818295,1847793,1848355],
                  [196589,231532,231557,264379,264569,388476],
                  [41301,447649,5320,581141,581180,649861],
                  [706710,743569,768257,795912,796175,831001],
                  [861036,861617,918636,94638,94687]]

    midlist = ['middle'+str(i) for i in [1,2,3,4,5,6,7,8]]
    idlists = [None,None,None]+mididlists

    for whichfolder,idlist in zip(['low','high','middle']+midlist,idlists):
        fig1,fig2 = plot_contam_dist(whichfolder,haloidlist=idlist)
        fig1.savefig('contamdist'+whichfolder+'.png',bbox_inches='tight')
        fig2.savefig('contamsize'+whichfolder+'.png',bbox_inches='tight')
        plt.close("all")
