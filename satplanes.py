import numpy as np
import pylab as plt
import haloutils
import cPickle as pickle
import pandas as pd

from caterpillaranalysis import PluginBase
from MTanalysis2 import ExtantDataFirstPass

import scipy.linalg as linalg

class SatellitePlanes(PluginBase):
    def __init__(self):
        super(SatellitePlanes,self).__init__()
        self.filename='satplanes.dat'
        
        self.xmin = 0; self.xmax = 1
        self.ymin = 0; self.ymax = 1
        self.xlabel = 'b/a'
        self.ylabel = 'c/a'
        self.xlog=False; self.ylog=False
        self.autofigname='satplanes'

        self.extant = ExtantDataFirstPass()
        self.maxsats = 25

    def calc_eig(self,satpos):
        satdist2 = np.sum((satpos)**2,1)

        I = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                I[i,j] = np.sum((satpos[:,i]*satpos[:,j])/satdist2)
        evals, evecs = linalg.eigh(I)
        return evals,evecs
    
    def get_n_largest_sats(self,numsats,subs,extsubs):
        vinfall = extsubs['infall_vmax'].copy()
        vinfall.sort(ascending=False)
        
        satids = vinfall[0:numsats].index
        extsats = extsubs.ix[satids]
        satrsids = np.array(extsats['rsid'])
        sats = subs.ix[satrsids]
        return sats
        
    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
        extdat = self.extant.read(hpath)

        numsnaps = haloutils.get_numsnaps(hpath)
        rscat = haloutils.load_rscat(hpath,numsnaps-1)
        zoomid = haloutils.load_zoomid(hpath)
        hpos = np.array(rscat.ix[zoomid][['posX','posY','posZ']])
        hvel = np.array(rscat.ix[zoomid][['pecVX','pecVY','pecVZ']])
        subs = self.get_rssubs(rscat,zoomid)
        subrsid = np.array(subs.index)
        
        extsubs = extdat #[extdat['rsid'].isin(subrsid)] #already filtered

        sats = self.get_n_largest_sats(self.maxsats,subs,extsubs)
        extsubs.index = extsubs['rsid']

        infall_labels = ['rsid','snap','vmax','mvir','posx','posy','posz',
                         'pecvx','pecvy','pecvz','virialratio','hostid_MT',
                         'rvir','spinbullock','rs','scale_of_last_MM',
                         'Jx','Jy','Jz','xoff']
        infall_labels = ['infall_'+l for l in infall_labels]
        for col in infall_labels:
            assert col not in sats.columns
            sats[col] = extsubs.ix[sats.index][col]
        sats.sort(columns='infall_vmax',ascending=False)
        infall_scale = pd.Series([haloutils.get_scale_snap(hpath,snap) for snap in sats['infall_snap']],index=sats.index)
        infall_z = pd.Series([haloutils.get_z_snap(hpath,snap) for snap in sats['infall_snap']],index=sats.index)
        sats['infall_scale'] = infall_scale; sats['infall_z'] = infall_z

        evallist = []; eveclist = []
        for numsats in range(1,self.maxsats+1):
            thissats = sats.iloc[0:numsats]
            satpos = np.array(thissats[['posX','posY','posZ']])-hpos
            satvel = np.array(thissats[['pecVX','pecVY','pecVZ']])-hvel
            evals,evecs = self.calc_eig(satpos)
            evallist.append(evals); eveclist.append(evecs)

        assert numsats==self.maxsats and len(satpos)==numsats and len(satvel)==numsats
        satL = np.cross(satpos,satvel)
        sats['Lx'] = satL[:,0]
        sats['Ly'] = satL[:,1]
        sats['Lz'] = satL[:,2]
        sats['dx'] = satpos[:,0]
        sats['dy'] = satpos[:,1]
        sats['dz'] = satpos[:,2]
        sats['dvx'] = satvel[:,0]
        sats['dvy'] = satvel[:,1]
        sats['dvz'] = satvel[:,2]

        with open(self.get_outfname(hpath),'w') as f:
            pickle.dump([evallist,eveclist,sats],f)
    def _read(self,hpath):
        with open(self.get_outfname(hpath),'r') as f:
            evallist,eveclist,sats = pickle.load(f)
        return sats,evallist,eveclist
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,**kwargs):
        sats,evallist,eveclist = data
        raise NotImplementedError

def tab_c_over_a(hpath):
    if hpath==None: return None
    plug = SatellitePlanes()
    out = plug.read(hpath,recalc=True)
    if out == None: return None
    sats,evallist,eveclist = out
    evallist = np.array(evallist)
    ca_arr = evallist[:,0]/evallist[:,2]
    ba_arr = evallist[:,1]/evallist[:,2]
    data = tuple(ca_arr)
    names = ['ca'+str(i+1) for i in range(len(ca_arr))]
    formats = [np.float for i in range(len(ca_arr))]
    return data,names,formats

def tab_rotation(hpath):
    if hpath==None: return None
    plug = SatellitePlanes()
    out = plug.read(hpath)
    if out == None: return None
    sats,evallist,eveclist = out

    spos = np.array(sats[['dx','dy','dz']])
    svel = np.array(sats[['dvx','dvy','dvz']])
    sangmom = np.cross(spos,svel)

    output = []
    for i in range(len(sats)):
        evals = evallist[i]
        V = eveclist[i] #V[:,j] = (j+1)th eigenvector
        #Vinv = linalg.inv(V) #Change of basis matrix
        #newpos = np.dot(Vinv,spos.T).T
        #newvel = np.dot(Vinv,svel.T).T
        proj = np.dot(sangmom[0:(i+1)],V[:,0])
        numzero = np.sum(np.abs(proj) < 1e-12)
        if numzero > 0:
            print "Warning: {0}/{1} satellites have 0 angular momentum ".format(numzero,i+1)
        num_coherent = max(np.sum(proj>0),i+1-np.sum(proj>0))
        output.append(num_coherent)
    data = tuple(output)
    names = ['f'+str(i+1) for i in range(len(sats))]
    formats = [np.int for i in range(len(sats))]
    return data,names,formats

if __name__=="__main__":
    plug = SatellitePlanes()
    #df = haloutils.tabulate(tab_c_over_a,numprocs=1)
    #fig,ax = plt.subplots(figsize=(8,8))
    #hids = df.index
    #for hid in hids:
    #    ca_arr = np.array(df.ix[hid])
    #    ax.plot(np.arange(plug.maxsats)+1,ca_arr,color='gray',label=haloutils.hidstr(hid))
    #ax.set_ylim((0,1))
    #ax.plot([0,25],[0.18,0.18],'r:')
    #ax.plot([11,11],[0,1],'k:')
    #plt.savefig('halos_c_a.png',bbox_inches='tight')

    numsats = np.arange(plug.maxsats)+1
    df = haloutils.tabulate(tab_rotation,numprocs=1)
    fig,ax = plt.subplots(figsize=(8,8))
    #plt.show()
    #dummy = raw_input()
    hids = df.index
    ax.plot(numsats,numsats*0.5,'k-')
    ax.plot(numsats,numsats*0.6,'k-')
    ax.plot(numsats,numsats*0.7,'k-')
    ax.plot(numsats,numsats*0.8,'k-')
    ax.plot(numsats,numsats*0.9,'k-')
    ax.plot(numsats,numsats,'k-')
    ax.set_xlabel('num satellites')
    ax.set_ylabel('num coherently rotating satellites')
    ax.plot([11,11],[0,25],'k:')
    ax.plot([0,25],[8,8],'k:')
    for hid in hids:
        numcoherent = np.array(df.ix[hid])
        #fracarr = np.array(df.ix[hid]).astype(np.float)/(np.arange(plug.maxsats)+1)
        l, = ax.plot(numsats,numcoherent,color='r',alpha=.9,label=haloutils.hidstr(hid))

        #ax.set_title(haloutils.hidstr(hid))
        #plt.draw()
        #dummy = raw_input()
        #l.remove()
    plt.show()
