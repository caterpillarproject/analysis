import pickle
from astropy.io import ascii
import haloutils

infilename='geometries.dat'
outfilename='geometries.p'

def bestgeom2geom(bestgeom):
    assert len(bestgeom)==3
    return bestgeom[0:2]+'_'+bestgeom[-1]

if __name__=="__main__":
    tab = ascii.read(infilename)
    tab = tab[tab['status']==1]
    geom = dict([(haloutils.hidstr(hid),bestgeom2geom(bestgeom)) for (hid,bestgeom) in zip(tab['hid'],tab['bestgeom'])])
    with open(outfilename,'w') as f:
        pickle.dump(geom,f)
