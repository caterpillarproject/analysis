import numpy as np
import annika_shapes as shapes
from scipy import integrate
import haloutils

def distance(posA, posB,boxsize=100.):
    dist = abs(posA-posB)
    tmp = dist > boxsize/2.0
    dist[tmp] = boxsize-dist[tmp]
    if dist.shape == (3,):
        return np.sqrt(np.sum(dist**2))
    else:
        return np.sqrt(np.sum(dist**2,axis=1))


# need to load all DM particles from main host at z=0
# then also load all star particles from halo at z=0.
# will need to weigh particles according to their mass.
def getShape(hpath, radius=None):
    snap_z0 = haloutils.get_numsnaps(hpath)-1
    cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
    hostID = haloutils.load_zoomid(hpath)
    center = np.array(cat.ix[hostID][['posX','posY','posZ']])
    
    if radius==None:
        radius = 1.25*float(cat.ix[hostID]['rvir'])/cat.h0
    
    posns = haloutils.load_partblock(hpath, snap_z0, "POS ")
    dist = distance(center,posns)*cat.scale/cat.h0*1000
    mask = dist<(radius)
    pos = posns[mask]
    pos2 = (pos-center)*cat.scale/cat.h0*1000 # in kpc physical
    arr_in = pos2
    ratios, evecs = shapes.axis(arr_in,radius,shell=False,axes_out=True) # this takes a long time. best to write data to file
    return ratios, evecs

# ratios = 2 element list c,b.  c<b<a. and a = 1
# eigenvectors are the normal unit vector axes associated with c,b and a.
# let x' be the axis associated with a, the semi-major axis.
# enforce that x' cross y' = z'. to determine which axis is y'
def getAngle(vec1, vec2):
    angle = np.arccos(np.dot(vec1,vec2))*180/np.pi
    return np.min([angle, 180.-angle])

