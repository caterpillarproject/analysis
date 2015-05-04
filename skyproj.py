import numpy as np
import astropy.coordinates as coord
import astropy.units as u

import pylab as plt
from matplotlib.widgets import Slider

import rotations

def is_linearly_independent(vecs,tol=None):
    vecs = np.vstack(vecs)
    N,d = vecs.shape
    if N > d: return False
    rank = np.linalg.matrix_rank(vecs.T,tol=tol)
    return rank == N

def gramschmidt(X, row_vecs=True, norm = True):
    if not row_vecs:
        X = X.T
    Y = X[0:1,:].copy()
    for i in range(1, X.shape[0]):
        proj = np.diag((X[i,:].dot(Y.T)/np.linalg.norm(Y,axis=1)**2).flat).dot(Y)
        Y = np.vstack((Y, X[i,:] - proj.sum(0)))
    if norm:
        Y = np.diag(1/np.linalg.norm(Y,axis=1)).dot(Y)
    if row_vecs:
        return Y
    else:
        return Y.T

def rotate_to_galactocentric(x_sun,Ldisk):
    """
    Returns a matrix that rotates coordinates such that Ldisk is (0, 0, -L) 
    and x_sun is at (-Rsun, 0, zsun) 

    @param x_sun: the coordinate where the sun is relative to the GC
    @param Ldisk: direction of disk rotation (opposite the direction of galactic north)
    @return Rmat: Rmat.dot(pos) is in coordinates where x_sun points in -x direction and L points in -z direction
    """
    assert is_linearly_independent([x_sun,Ldisk])
    Lhat = Ldisk/np.linalg.norm(Ldisk)
    zhat = -1.0 * Lhat
    R_L_to_z = rotations.rotate_to_z(zhat)

    rotsunx,rotsuny,rotsunz = R_L_to_z.dot(x_sun)
    theta_z = np.arccos(rotsunx/np.sqrt(rotsunx**2 + rotsuny**2))
    if rotsuny > 0: theta_z = -theta_z
    R_sun_to_negx = rotations.rotmat_z(theta_z+np.pi)

    Rmat = np.dot(R_sun_to_negx,R_L_to_z)
    return Rmat

def skyproj(spos, svel, Rsun, Ldisk, Vsun=220., numphi = 100):
    Ldisk = np.array(Ldisk); assert len(Ldisk)==3
    N,d = spos.shape; assert N,d == svel.shape; assert d==3

    # Get 3 linearly independent vectors
    v1 = Ldisk
    found_i = -1
    for i in range(N):
        v2 = spos[i,:]
        if is_linearly_independent([v1,v2]):
            found_i = i
            break
    if found_i == -1:
        raise RuntimeError("No linearly independent vectors!")
    v2 = spos[found_i,:]
    v3 = np.cross(v1,v2)
    X = np.array([v1,v2,v3]) #rows are vectors
    B = gramschmidt(X) #orthonormal basis
    B1 = B[0,:]; B2 = B[1,:]; B3 = B[2,:]

    phiarr = np.linspace(0,2*np.pi,num=numphi,endpoint=False)
    output = []
    for phi in phiarr:
        x_sun = Rsun * (B2 * np.cos(phi) + B3 * np.sin(phi))

        L_hat = B1
        R_hat = x_sun/np.linalg.norm(x_sun)
        phi_hat = np.cross(L_hat,R_hat)
        v_sun = Vsun * phi_hat #km/s

        this_spos = spos - x_sun
        this_svel = svel - v_sun
        output.append([x_sun,this_spos,this_svel])
        #this_sdist = np.linalg.norm(this_spos,axis=1)
        #this_spos_normed = (this_spos.T / this_sdist).T
        #this_vrad = np.sum(this_svel*this_spos_normed,1)

    return phiarr,output

if __name__=="__main__":
    import haloutils
    from satplanes import SatellitePlanes
    hpath = haloutils.get_hpath_lx(1725139,14)
    plug = SatellitePlanes()
    sats,evallist,eveclist = plug.read(hpath)

    spos = np.array(sats[['dx','dy','dz']])
    svel = np.array(sats[['dvx','dvy','dvz']])
    Ldisk = np.cross(spos[0,:],svel[0,:])
    Rsun = .67 * 8.3/1000. #Mpc/h
    phiarr,output = skyproj(spos, svel, Rsun, Ldisk)

    #    # z ~ -Ldisk, x_sun ~ -x
    #    Rmat = rotate_to_galactocentric(x_sun,Ldisk)
    #    rot_x_sun = Rmat.dot(x_sun)
    #    rot_spos  = Rmat.dot(spos.T).T
#
#        x_sun= x_sun*1000./.67 * u.kpc
#        spos = spos*1000./.67 * u.kpc
    gcframe = coord.Galactocentric(z_sun=0*u.kpc)

    coordlist = []
    for i,(x_sun,spos,svel) in enumerate(output):
        Rmat = rotate_to_galactocentric(x_sun,Ldisk)
        rot_spos = Rmat.dot(spos.T).T * (1000./.67) * u.kpc
        rot_svel = Rmat.dot(svel.T).T * u.km/u.s
        satcoords = coord.SkyCoord(x=rot_spos[:,0],y=rot_spos[:,1],z=rot_spos[:,2],frame=gcframe)
        coordlist.append(satcoords)
        RA = satcoords.icrs.ra.wrap_at(180*u.degree)
        DEC= satcoords.icrs.dec
        L = satcoords.galactic.l.wrap_at(180*u.degree)
        B = satcoords.galactic.b

        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111,projection='mollweide')
        ax.scatter(RA.radian,DEC.radian)
        fig.savefig('4-30/mollweideRADEC_H1725139_'+str(i).zfill(3)+'.png',bbox_inches='tight')
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111,projection='mollweide')
        ax.scatter(L.radian,B.radian)
        fig.savefig('4-30/mollweideGALAC_H1725139_'+str(i).zfill(3)+'.png',bbox_inches='tight')
        plt.close('all')
        
