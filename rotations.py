import numpy as np

def rotmat_x(t):
    return np.array([[1,0,0],[0,np.cos(t),-np.sin(t)],[0,np.sin(t),np.cos(t)]])
def rotmat_y(t):
    return np.array([[np.cos(t),0,np.sin(t)],[0,1,0],[-np.sin(t),0,np.cos(t)]])
def rotmat_z(t):
    return np.array([[np.cos(t),-np.sin(t),0],[np.sin(t),np.cos(t),0],[0,0,1]])

def rotate_to_z(vec):
    """
    Creates a rotation matrix M where dot(M,vec) \propto [0,0,1]
    """
    assert len(vec)==3; vec = np.array(vec)
    u,v,w = vec/np.linalg.norm(vec)
    eps = np.finfo(float).eps
    if np.abs(u)<eps and np.abs(v) < eps:
        if w > 0: return np.eye(3)
        else: return rotmat_x(np.pi)
    #Rotate around z-axis to x-z plane
    theta_z = -np.arctan(v/u)
    if u < 0: theta_z += np.pi
    Rz = rotmat_z(theta_z)
    #Rotate around y-axis to z
    theta_y = -np.arccos(w)
    Ry = rotmat_y(theta_y)
    return np.dot(Ry,Rz)

def rotate_to_x(vec):
    """Creates a rotation matrix M where dot(M,vec) = [1,0,0]"""
    R2z = rotate_to_z(vec)
    Ry = rotmat_y(np.pi/2.)
    return np.dot(Ry,R2z)
def rotate_to_y(vec):
    """Creates a rotation matrix M where dot(M,vec) = [0,1,0]"""
    R2z = rotate_to_z(vec)
    Rx = rotmat_x(-np.pi/2.)
    return np.dot(Rx,R2z)

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
def xyz2thetaphi(pos,rotate_for_mollweide=False):
    """
    Puts [0,0,1] as pole by default.
    If rotate_for_mollweide, puts [1,0,0] at (theta,phi)=0,0 and [0,0,1] as pole
    """
    x = pos[:,0]; y = pos[:,1]; z = pos[:,2]
    rxy = np.sqrt(x**2 + y**2)
    phi = np.arccos(x/rxy)
    phi[np.isnan(phi)] = 0
    phi[y<0] = 2*np.pi - phi[y<0]
    phi = normalize_phi(phi)
    r = np.sqrt(np.sum(pos**2,1))
    theta = np.arccos(z/r)
    if rotate_for_mollweide:
        phi[phi > np.pi] = phi[phi > np.pi]-2*np.pi
        theta = np.pi/2. - theta
        #rxz = np.sqrt(x**2 + z**2)
        #theta = np.arccos(np.abs(x)/rxz)
        #theta[rxz==0] = 0
        #theta[z<0] = -theta[z<0]
    return theta,phi

