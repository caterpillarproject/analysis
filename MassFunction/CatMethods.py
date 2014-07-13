import numpy as np
import os
import readsubf

# compute distance from posA to posB.
# posA can be an array. boxsize must be in same units as positions.
def distance(posA, posB,boxsize=25.):
    dist = abs(posA-posB)
    tmp = dist > boxsize/2.0
    dist[tmp] = boxsize-dist[tmp]
    return np.sqrt(np.sum(dist**2,axis=1))

def row_magnitude(matrix):
    """
    Find magnitude of each row.
    @ param matrix: m x n matrix.
    @ return: m x 1 column vector of magnitudes.
    """
    return np.sqrt(sum((matrix**2).T))[:,np.newaxis]

### handles row vectors with 0 magnitude correctly
def row_norm(matrix):
    magnitude = row_magnitude(matrix)
    matrix = matrix.astype('d') #make sure array is all floats
    return np.nan_to_num(matrix/magnitude)

def row_dot(a,b):
    return sum((a*b).T)[:,np.newaxis]

def getMidpoints(bins):
    """
    Given a range of cut-off values, find the midpoints.
    @return: array of length len(bins)-1
    """
    spacing = bins[1:]-bins[:-1]
    return bins[:-1]+spacing/2.0

def getaqparams(halo):
    if halo == "Aq-A":
        snap = 1023
        swapval = False
        long_ids = True
    elif halo == "Aq-F":
        snap = 111
        swapval = False
        long_ids = False
    elif halo == "Aq-E":
        snap = 127
        swapval = True
        long_ids = False
    else:
        snap = 127
        swapval = False
        long_ids = True

    return snap,swapval,long_ids

def getAQcatalogs():
    halolist = ["Aq-A","Aq-B","Aq-C","Aq-D","Aq-E","Aq-F"]
    catlist = []
    for haloi in halolist:
        maxsnap,swapval,long_ids = getaqparams(haloi)
        path = "/bigbang/data/AnnaGroup/aquarius/data/" + haloi + "/2/"
        if os.path.isdir(path+'groups_'+str(maxsnap)):
            catlist.append(readsubf.subfind_catalog(path, maxsnap,long_ids = long_ids, swap = swapval))
    return catlist 
