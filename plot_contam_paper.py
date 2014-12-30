import numpy as np
import pylab as plt
import haloutils
import asciitable
import os
from findhalos.caterpillar_findhalos import get_zoom_id
#import seaborn.apionly as sns
import seaborn as sns
sns.set_context('notebook')
sns.set_style('ticks')

from caterpillarplot import get_haloidlist

def plot_contam_dist(whichfolderlist,haloidlist):
    maxE = 5
    labelmap = {'CA4':-2,'CA5':-1,
                'EA4':0,'EA5':1,'EX4':2,'EX5':3,'EB':4,'EC':5,
                'BA4':(maxE+1),'BA5':(maxE+2),'BB':(maxE+3),'BC':(maxE+4),'BD':(maxE+5)}

    fig,axarr = plt.subplots(2,1,figsize=(6,10))
    ax1 = axarr[0]; ax2 = axarr[1]
    
    colors = sns.color_palette('Set3',n_colors=len(whichfolderlist))
    for whichfolder,hid,color in zip(whichfolderlist,haloidlist,colors):
        assert whichfolder in ['low','middle','high']
        data = asciitable.read('/bigbang/data/AnnaGroup/caterpillar/halos/'+whichfolder+'_mass_halos/contam_zoom_index.txt',Reader=asciitable.FixedWidth)
        hidlist = np.unique(data['parentid'])
        assert hid in hidlist
        thisdat = data[data['parentid']==hid]
        contamtypes = thisdat['ictype']
        nv = thisdat['NV']
        tmplist = []
        for i in range(len(contamtypes)):
            if contamtypes[i] in ['BA','EA','CA','EX']:
                tmplist.append(contamtypes[i]+str(nv[i]))
            else: tmplist.append(contamtypes[i])
        contamtypes = tmplist
        contamx = np.array([labelmap[ct] for ct in contamtypes])
        iisort = np.argsort(contamx)
        contamx = contamx[iisort]
        mindist = np.array(thisdat['min2'])[iisort]/.6711
        ax1.plot(contamx,mindist,'o-',label=haloutils.hidstr(hid),color=color)
        filesizes = np.array(thisdat['icsize'])[iisort]
        ax2.plot(contamx,filesizes,'o-',label=haloutils.hidstr(hid),color=color)
        ax1.legend(loc='upper center',bbox_to_anchor=(.5,1.05),ncol=3,fancybox=True,fontsize='small')
    minx = np.min(labelmap.values()); maxx = np.max(labelmap.values())

    ax1.plot([minx-0.5,maxx+0.5],[1.0/.6711,1.0/.6711],'k:')
    ax1.set_ylabel('Contamination Distance [Mpc]')
    ax2.set_ylabel('Size of LX11 ics [MB]')

    #ax1.set_ylim([0,1.8])
    #ax2.set_ylim([100,700])
    for ax in [ax1,ax2]:
        ax.set_xlim([minx-0.5,maxx+0.5])
        ax.set_xticks(np.sort(labelmap.values()))
    ax1.set_xticklabels(['' for i in range(len(['CA4','CA5','EA4','EA5','EX4','EX5','EB','EC','BA4','BA5','BB','BC','BD']))])
    ax2.set_xticklabels(['CA4','CA5','EA4','EA5','EX4','EX5','EB','EC','BA4','BA5','BB','BC','BD'])
    fig.subplots_adjust(hspace=.05)
    return fig

if __name__=="__main__":
    hidlist = get_haloidlist(1)
    whichfolderlist = ['middle' for i in range(12)]
    whichfolderlist[1] = 'low'

    fig = plot_contam_dist(whichfolderlist,hidlist)
    fig.savefig('paper_contam.pdf',bbox_inches='tight')
    plt.close('all')



    #mididlists = [[1079897,1129843,1130025,1232127,1232164,1268839],
    #              [1292085,1327666,1354437,1386703,1387186,1422331],
    #              [1599902,1599988,1631506,1666505,1697496], #1452004,
    #              [1725139,1725272,1725372,1818295,1847793,1848355],
    #              [196589,231532,231557,264379,264569,388476],
    #              [41301,447649,5320,581141,581180,649861],
    #              [706710,743569,768257,795912,796175,831001],
    #              [861036,861617,918636,94638,94687]]

    #midlist = ['middle'+str(i) for i in [1,2,3,4,5,6,7,8]]
    #idlists = [None,None,None]+mididlists

    #for whichfolder,idlist in zip(['low','high','middle']+midlist,idlists):
    #    fig1,fig2 = plot_contam_dist(whichfolder,haloidlist=idlist)
    #    fig1.savefig('contamdist'+whichfolder+'.png',bbox_inches='tight')
    #    fig2.savefig('contamsize'+whichfolder+'.png',bbox_inches='tight')
    #    plt.close("all")
