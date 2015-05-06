import matplotlib; matplotlib.use('Agg')
import numpy as np
import pylab as plt
import os,subprocess,sys,time
import cPickle as pickle
import pandas as pd

import haloutils
from caterpillaranalysis import PluginBase
from MTanalysis2 import ExtantDataFirstPass

import numpy.random as random
from scipy import linalg, stats, integrate, interpolate, optimize
import functools,itertools
from multiprocessing import Pool
import rotations
from astropy.coordinates import SkyCoord
from astropy import units as u

def angular_dists(pos):
    """
    Given positions, returns an array of angular distances in radians
    """
    N = len(pos)
    x = pos[:,0]; y = pos[:,1]; z = pos[:,2]
    c = SkyCoord(x=x,y=y,z=z,frame='galactocentric')
    angdists = np.zeros(N*(N-1)/2)
    count = 0
    for i in range(N):
        for j in range(i+1,N):
            angdists[count] = c[i].separation(c[j]).radian
            count += 1
    return angdists
def angular_correlation(pos,nbins=20):
    dists = angular_dists(pos)
    cosd = np.cos(dists)
    bins = np.linspace(-1,1,nbins+1)
    h,x = np.histogram(cosd,bins=bins)
    return h,x

class DefaultRadialDistn(stats.rv_continuous):
    """
    Zentner et al. 2005 radial distr.
    """
    def __init__(self,a=0.,b=1.):
        super(DefaultRadialDistn,self).__init__(a=a,b=b)
        self.norm = np.log((4*self.b)**3 + 1.0)/3.
    def _pdf(self,y):
        return (4*y)**2 / (self.norm*(1+(4*y)**3))
    def _cdf(self,y):
        return np.log((4*y)**3 + 1.0)/(3.*self.norm)
    def _stats(self):
        return np.nan,np.nan,np.nan,np.nan

def fit_plane(pos):
    x = pos[:,0]; y = pos[:,1]; z = pos[:,2]
    X = np.vstack([x,y]).T
    beta = linalg.inv(np.dot(X.T,X)).dot(X.T).dot(z)
    a,b = beta
    v1 = np.array([1.,0.,a])
    v2 = np.array([0.,1.,b])
    n = np.cross(v1,v2)
    n = n/linalg.norm(n)
    return n
def Drms(pos,n):
    return np.sqrt(np.sum(np.dot(n,pos)**2)/float(len(pos)))
def calc_cosw(pos,n):
    n = n/linalg.norm(n)
    d = np.sqrt(np.sum(pos**2,1))
    pos = pos/d[:,np.newaxis]
    cosw = np.abs(np.dot(pos,n))
    return cosw
def get_cosw_cdf(cosw,nbins=100):
    bins = np.linspace(0,1,nbins+1)
    h,x = np.histogram(cosw,bins=bins)
    return h.cumsum()/float(np.sum(h))
class IsotropicSatellites(object):
    def __init__(self,distn=None):
        if distn==None:
            self.distn = DefaultRadialDistn()
        else:
            self.distn = distn

    def create_random_sats(self,Nsats,N):
        return map(self.draw_random_points,itertools.repeat(Nsats,N))

    def plane_c_over_a(self,pos,satplug):
        evals,evecs = satplug.calc_eig(pos)
        c_over_a = np.sqrt(evals[0]/evals[2])
        return c_over_a
    def calc_c_over_a(self,satlist):
        satplug = SatellitePlanes()
        mapfn = functools.partial(self.plane_c_over_a,satplug=satplug)
        return np.array(map(mapfn,satlist))

    def plane_cosw_cdf(self,pos):
        n = fit_plane(pos)
        cosw = calc_cosw(pos,n)
        return get_cosw_cdf(cosw)
    def calc_cosw_cdf(self,satlist):
        cdflist = map(self.plane_cosw_cdf,satlist)
        return np.array(cdflist)

    def draw_random_points(self,N):
        """
        Draw points isotropically distributed on a sphere,
        with a random distance drawn from self.distn
        IMPORTANT: need to multiply by Rvir to get actual coordinates
        """
        points = self.uniform_sphere_points(N)
        radii = self.distn.rvs(size=N)
        return points*radii[:,np.newaxis]
    def uniform_sphere_points(self,N,d=3):
        randnorm = random.randn(N,d)
        r = np.sqrt(np.sum(randnorm**2,1))
        return (randnorm.T / r).T
    def uniform_3sphere_angles(self,N):
        randunif = random.rand(N,2)
        randunif[:,0] = np.arccos(2*randunif[:,0]-1)-np.pi/2.
        randunif[:,1] = 2*np.pi*randunif[:,1]
        return randunif

def plot_isotropic_Nsats(maxNsats=25,ax=None):
    iso = IsotropicSatellites()
    N = 10**5
    Nsatlist = np.arange(maxNsats)+1
    sigma_lo = []
    medians = []
    sigma_hi = []
    for Nsats in Nsatlist:
        Nsats = int(Nsats)
        filename = 'isosats_c_a_'+str(Nsats)+'.npy'
        if os.path.exists(filename):
            c_a = np.load('isosats_c_a_'+str(Nsats)+'.npy')
        else:
            print "Computing Nsats",Nsats,'...'
            start = time.time()
            satlist = iso.create_random_sats(Nsats,int(N))
            print "  Created random sats. Time={0:.1f}".format(time.time()-start)
            c_a = iso.calc_c_over_a(satlist)
            np.save(filename,c_a)
            print "Done! Time={0:.1f}".format(time.time()-start)
        lo,med,hi = np.percentile(c_a,[17,50,83])
        sigma_lo.append(lo); medians.append(med); sigma_hi.append(hi)
    if ax==None:
        fig,ax = plt.subplots(figsize=(8,8))
    ax.plot(Nsatlist,np.array(sigma_lo),'k:')
    ax.plot(Nsatlist,np.array(medians),'k--')
    ax.plot(Nsatlist,np.array(sigma_hi),'k:')
    return ax

def plot_isotropy(Nsatlist = [6,11,15,25]):
    iso = IsotropicSatellites()
    N = 10000.
    avgs = []
    for Nsats in Nsatlist:
        satlist = iso.create_random_sats(Nsats,int(N))
        cdfs = iso.calc_cosw_cdf(satlist)
        c_a  = iso.calc_c_over_a(satlist)
        np.save('isosats_cdf_'+str(Nsats)+'.npy',cdfs)
        np.save('isosats_c_a_'+str(Nsats)+'.npy',c_a)
        avg = np.sum(cdfs,0)/float(N)
        avgs.append(avg)
    #cdfs = []
    #c_as = []
    #for Nsats in Nsatlist:
    #    cdfs.append(np.load('isosats_cdf_'+str(Nsats)+'.npy'))
    #    #c_as.append(np.load('isosats_c_a_'+str(Nsats)+'.npy'))
    fig,ax = plt.subplots(figsize=(8,8))
    ax.set_xlabel(r'$|\cos{\omega}|$')
    ax.set_ylabel(r'$f_{\rm sat}(<|\cos{\omega}|)$')
    x_cosw = np.linspace(0,1,len(avgs[0]))
    for i,Nsats in enumerate(Nsatlist):
        ax.plot(x_cosw,avgs[i],label=r'$N_{{\rm sat}}={0}$'.format(Nsats))
    ax.legend(loc='lower right')
    fig.savefig('isotropic_cosw_cdf.png',bbox_inches='tight')

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
    out = plug.read(hpath)
    if out == None: return None
    sats,evallist,eveclist = out
    evallist = np.array(evallist)
    ca_arr = np.sqrt(evallist[:,0]/evallist[:,2])
    ba_arr = np.sqrt(evallist[:,1]/evallist[:,2])
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

def plot_ang_mom(Ldisk = [0,0,1]):
    plug = SatellitePlanes()
    hids = haloutils.cid2hid.values()
    Ldisk = np.array(Ldisk) #arbitrary for now
    rotmat = rotations.rotate_to_z(Ldisk)
    output = []
    for hid in hids:
        hpath = haloutils.get_hpath_lx(hid,14)
        if hpath==None: continue
        out = plug.read(hpath)
        if out == None: continue
        print haloutils.hidstr(hid)
        sats,evallist,eveclist = out
        spos = np.array(sats[['dx','dy','dz']])
        svel = np.array(sats[['dvx','dvy','dvz']])
        Lsat = np.array(sats[['Lx','Ly','Lz']])
        Ltot = np.linalg.norm(Lsat,axis=1)
        #r_spos = rotmat.dot(spos.T).T
        #r_svel = rotmat.dot(svel.T).T
        r_Lsat = rotmat.dot(Lsat.T).T
        rLx = r_Lsat[:,0]; rLy = r_Lsat[:,1]; rLz = r_Lsat[:,2]
        theta = np.arccos(rLz/Ltot)
        phi = np.arccos(rLx/np.sqrt(rLx**2 + rLy**2))
        phi[rLy < 0] = 2*np.pi - phi[rLy < 0]

        theta = theta - np.pi/2.
        phi = phi-np.pi

        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111,projection='mollweide')
        ax.scatter(phi,theta)
        fig.savefig('5-1/mollweideL_'+haloutils.hidstr(hid)+'.png',bbox_inches='tight')
        plt.close('all')
        h,x = angular_correlation(Lsat)
        output.append([hid,h,x,r_Lsat])
    with open('Lmom.p','w') as f:
        pickle.dump(output,f)

def plot_mollweide(sats,host,fig=None,subplots=111,useA=False):
    if fig==None:
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='mollweide')
    else:
        ax = fig.add_subplot(subplots,projection='mollweide')

    if useA:
        halo_z = np.array(host[['A2[x]','A2[y]','A2[z]']])
    else:
        halo_z = np.array(host[['Jx','Jy','Jz']])
    rotmat = rotations.rotate_to_z(halo_z)
    satsL   = np.array(sats[['Lx','Ly','Lz']])

    r_satsL = rotmat.dot(satsL.T).T
    sats_theta,sats_phi = rotations.xyz2thetaphi(r_satsL)
    sats_theta -= np.pi/2; sats_phi -= np.pi

    #halo_z = rotmat.dot(halo_z)
    #halo_theta,halo_phi = rotations.xyz2thetaphi(halo_z[np.newaxis,:])
    #halo_theta -= np.pi/2; halo_phi -= np.pi

    #ax.plot(halo_phi,halo_theta,'o',color='red',mec=None,markersize=16)
    ax.plot(sats_phi,sats_theta,'o',color='blue',mec=None,markersize=16)
    return ax

if __name__=="__main__":
    #plot_isotropy()
    #plot_isotropic_Nsats()
    #plot_ang_mom()
    hids = haloutils.cid2hid.values()
    plug = SatellitePlanes()
    for hid in hids:
        print haloutils.hidstr(hid)
        hpath = haloutils.get_hpath_lx(hid,14)
        if hpath==None: continue
        out = plug.read(hpath)
        if out==None: continue
        sats,evals,evecs = out
        rscat = haloutils.load_rscat(hpath,haloutils.get_numsnaps(hpath)-1)
        zoomid= haloutils.load_zoomid(hpath)
        host = rscat.ix[zoomid]
        fig = plt.figure()
        plot_mollweide(sats,host,fig=fig,useA=False)
        fig.savefig('5-5/mollweide/LrotJ_'+haloutils.hidstr(hid)+'.png')
        fig = plt.figure()
        plot_mollweide(sats,host,fig=fig,useA=True)
        fig.savefig('5-5/mollweide/LrotA_'+haloutils.hidstr(hid)+'.png')
        plt.close('all')

def simple_plots():
    plug = SatellitePlanes()
    df = haloutils.tabulate(tab_c_over_a,numprocs=1)
    fig,ax = plt.subplots(figsize=(8,8))
    hids = df.index
    for hid in hids:
        ca_arr = np.array(df.ix[hid])
        ax.plot(np.arange(plug.maxsats)+1,ca_arr,color='gray',label=haloutils.hidstr(hid))
    ax.set_ylim((0,1))
    ax.plot([0,25],[0.18,0.18],'r:')
    ax.plot([11,11],[0,1],'k:')
    ax.set_xlabel(r'$\rm{num\ satellites}$')
    ax.set_ylabel(r'$c/a$')
    plt.savefig('halos_c_a.png',bbox_inches='tight')

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
