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

    #data = asciitable.read('contamdist.dat')
    #hidlabel = 'haloid'
    hidlabel = 'parentid'
    if whichfolder=='low':
        data = asciitable.read('/bigbang/data/AnnaGroup/caterpillar/halos/low_mass_halos/contam_zoom_index.txt',Reader=asciitable.FixedWidth)
    elif whichfolder=='middle':
        data = asciitable.read('/bigbang/data/AnnaGroup/caterpillar/halos/middle_mass_halos/contam_zoom_index.txt',Reader=asciitable.FixedWidth)
    elif whichfolder=='high':
        data = asciitable.read('/bigbang/data/AnnaGroup/caterpillar/halos/high_mass_halos/contam_zoom_index.txt',Reader=asciitable.FixedWidth)
    print whichfolder
    hidlist = np.unique(data[hidlabel])
    if haloidlist != None:
        hidlist = list(set(haloidlist).intersection(hidlist))

    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    for hid in hidlist:
        thisdat = data[data[hidlabel]==hid]
        contamtypes = thisdat['contamtype']
        nv = thisdat['NV']
        tmplist = []
        for i in range(len(contamtypes)):
            if contamtypes[i] in ['A','ELLIPSOID','CONVEX']:
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
    ax2.set_ylabel('size of ics.0 in MB')
    for ax in [ax1,ax2]:
        ax.set_xlim([minx-0.5,maxx+0.5])
        ax.set_xticks(np.sort(labelmap.values()))
        ax.set_xticklabels(['CV4','CV5','EL4','EL5','A4','A5','B','C','D','ELB','ELC'])
        ax.legend(loc='best')
    return fig1,fig2

if __name__=="__main__":
    for whichfolder in ['low','high']:
        fig1,fig2 = plot_contam_dist(whichfolder)
        fig1.savefig(whichfolder+'contamdist.png',bbox_inches='tight')
        fig2.savefig(whichfolder+'contamsize.png',bbox_inches='tight')
    plt.close("all")
