import matplotlib; matplotlib.use('Agg')
import numpy as np
import pylab as plt
import os,subprocess,sys,time
import cPickle as pickle
import pandas as pd
import glob

import haloutils
from rotations import xyz2thetaphi

from caterpillaranalysis import PluginBase,MassAccrPlugin
from scipy import special

lmin = 0
lmax = 10

def a_lm(l,m,theta,phi):
    #Note special.sph_harm documentation has opposite definitions of theta/phi
    return 1./(np.pi * 4.0) * np.sum(np.conj(special.sph_harm(m,l,phi,theta)))
def C_l(l,theta,phi,rvir):
    alm_arr = np.zeros(2*l+1,dtype=complex)
    for i in range(2*l+1):
        m = i-l
        alm_arr[i] = a_lm(l,m,theta,phi)
    sumsquares = np.sum(alm_arr * np.conj(alm_arr))
    return sumsquares/(4*np.pi*(2*l+1)*rvir**2)

def filter_points(hpos,hvel,rvir,ppos,pvel,rmin=1.0,rmax=1.2):
    """
    This function selects points in some distance bin that are infalling.
    Neglects hubble expansion.

    @param hpos: host halo pos
    @param hvel: host halo vel
    @param rvir: host halo virial radius (same units as hpos)
    @param ppos: particle pos (same units as hpos)
    @param pvel: particle vel (same units as hvel)
    @param rmin: minimum distance from hpos in units of rvir (default 1.0)
    @param rmax: maximum distance from hpos in units of rvir (default 1.2)

    @return: indices of ppos and pvel which satisfy the conditions
    """

    ppos -= hpos
    pvel -= hvel
    dist = np.sqrt(np.sum(ppos**2,1))
    ii_distcut = (dist >= rmin*rvir) & (dist <= rmax*rvir)
    ii_infall = np.sum(ppos*pvel,1) < 0

    return ii_distcut & ii_infall

minsize=200
maxsize=500

cmap = plt.get_cmap('cubehelix')
fig,ax = plt.subplots(figsize=(8,8))
hids = haloutils.cid2hid.values()
for hid in hids:
    print haloutils.hidstr(hid)
    hpath = haloutils.get_hpath_lx(hid,12)
    icsize = 0
    for icpath in glob.glob(hpath+'/ics.*'):
        icsize += os.path.getsize(icpath)
    icsize /= 1.e6 #in MB
    normicsize = (icsize-minsize)/(maxsize-minsize)
    with open('5-5/'+haloutils.hidstr(hid)+'_Clarr.p','r') as f:
        Cl_arr,mb = pickle.load(f)
    times = haloutils.get_t_snap(hpath,mb['snap'])
    ratio_2_1 = Cl_arr[2,mb['snap']]/Cl_arr[1,mb['snap']]
    ratio_2_1[ratio_2_1==0] = 10**-10
    ax.plot(times,ratio_2_1,color=cmap(normicsize),label=haloutils.hidstr(hid))
ax.set_yscale('log')
plt.show()

def old():
    Cl_arr = np.zeros((lmax-lmin+1,256),dtype=complex)
    mtplug = MassAccrPlugin()
    mb = mtplug.read(hpath)
    h0 = .6711
    all_hpos = np.array(mb[['x','y','z']]).view(float).reshape((len(mb),3))
    all_hvel = np.array(mb[['vx','vy','vz']]).view(float).reshape((len(mb),3))
    all_rvir_Mpch = h0*mb['rvir']/1000.
    start = time.time()
    for j,snap in enumerate(mb['snap']):
        if snap % 10 == 0:
            print "Snap {0}, Time {1:.1f}".format(snap,time.time()-start)
        hpos = all_hpos[j,:]
        hvel = all_hvel[j,:]
        rvir_Mpch = all_rvir_Mpch[j]
        
        ppos = haloutils.load_partblock(hpath,snap,"POS ",parttype=1)
        pvel = haloutils.load_partblock(hpath,snap,"VEL ",parttype=1)

        ii_cut = filter_points(hpos,hvel,rvir_Mpch,ppos,pvel)
        theta,phi = xyz2thetaphi(ppos[ii_cut,:])
        
        for i,l in enumerate(np.arange(lmin,lmax+1)):
            Cl_arr[i,snap] = C_l(l,theta,phi,rvir_Mpch)
    a = Cl_arr.astype(float)
    a = a/a[0,:]
    a = a[:,~np.isnan(a[0,:])]
    plt.figure()
    plt.imshow(np.log10(a),interpolation='none',aspect='auto')
    plt.colorbar()
    plt.xlabel('offset snap')
    plt.ylabel(r'$\ell$')
    plt.ylim(plt.gca().get_ylim()[::-1])
    fig = plt.gcf()
    with open('5-5/'+haloutils.hidstr(hid)+'_Clarr.p','w') as f:
        pickle.dump([Cl_arr,mb],f)
    fig.savefig('5-5/'+haloutils.hidstr(hid)+'_powspec.png',bbox_inches='tight')
