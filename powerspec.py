import matplotlib; matplotlib.use('Agg')
import numpy as np
import pylab as plt
import os,subprocess,sys,time
import cPickle as pickle
import pandas as pd

import haloutils

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

def normalize_phi(phi):
    allgood = False
    while not allgood:
        ii1 = phi>2*np.pi
        ii2 = phi<0
        if np.sum(ii1) > 0:
            phi[ii1] = phi[ii1]-2*np.pi
        elif np.sum(ii2) > 0:
            phi[ii2] = phi[ii2]-2*np.pi
        else:
            allgood=True
    return phi

def xyz2thetaphi(pos):
    r = np.sqrt(np.sum(pos**2,1))
    x = pos[:,0]; y = pos[:,1]; z = pos[:,2]
    rxy = np.sqrt(x**2 + y**2)
    theta = np.arccos(z/r)
    phi = np.arccos(x/rxy)
    phi[y<0] = 2*np.pi - phi[y<0]
    phi = normalize_phi(phi)
    return theta,phi

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

if __name__=="__main__":
    hpath = haloutils.get_hpath_lx(1725139,12)
    Cl_arr = np.zeros((lmax-lmin+1,256),dtype=complex)
    mtplug = MassAccrPlugin()
    mb = mtplug.read(hpath)
    h0 = .6711
    all_hpos = np.array(mb[['x','y','z']]).view(float).reshape((len(mb),3))
    all_hvel = np.array(mb[['vx','vy','vz']]).view(float).reshape((len(mb),3))
    all_rvir_Mpch = h0*mb['rvir']/1000.
    start = time.time()
    for j,snap in enumerate(mb['snap']):
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
