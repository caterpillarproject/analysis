import numpy as np
import readsnapshots.readsnap as rs
from scipy import interpolate

def densityprofile(rarr,filepath,header,haloparts,halopos,verbose=False,power03=False):
    """
    @param rarr: radius array in Mpc
    @param filepath: filepath to read_block (see readsnap.py)
    @param haloparts: list of halo particle IDs
    @param halopos: halo position in Mpc
    @param verbose: if True, prints things out
    @param power03: if True, also returns power03 convergence radius

    @return rhoarr: at every r in rarr, returns rho(r) in 10^10Msun/Mpc^3/Mpc
    @return p03rmin: if power03 is set to True, returns the power03 convergence radius (in Mpc)
    """
    rhoarr = np.zeros((len(rarr),5))
    if power03:
        Narr = np.zeros(len(rarr))
        marr = np.zeros(len(rarr))

    for parttype in [1,2,3,4,5]:
        snapIDs=rs.read_block(filepath,"ID  ",parttype=parttype,doubleprec=False)
        snapPOS=rs.read_block(filepath,"POS ",parttype=parttype,doubleprec=False)
        mask = np.in1d(snapIDs,haloparts,assume_unique=True) #for large numbers of particles, this is slow
        dx=snapPOS[mask][:,0]-halopos[0]
        dy=snapPOS[mask][:,1]-halopos[1]
        dz=snapPOS[mask][:,2]-halopos[2]
        dr=np.sort(np.sqrt(dx**2 + dy**2 + dz**2))
        
        if verbose:
            print "  Particle type",parttype
            print "  Number particles:",len(mask),len(dr)
            if len(dr) != 0:
                print "  dr range=",dr.min(),dr.max()
                if dr.max()>rarr[-1]:
                    nout = np.sum(dr > rarr[-1])
                    print "  densityprofile warning:",nout,"particles lie outside max(rarr)"

        if len(dr) != 0:
            h_r, x_r = np.histogram(dr, bins=np.concatenate(([0],rarr)))
            #Use Phillip's method of splining the cumulative array and taking the derivative
            m_lt_r = np.cumsum(h_r)*header.massarr[parttype]
            tck = interpolate.splrep(rarr,m_lt_r)
            rhoarr[:,parttype-1] = interpolate.splev(rarr,tck,der=1)
            if power03:
                Narr = Narr + np.cumsum(h_r)
                marr = marr + np.cumsum(h_r)*header.massarr[parttype]
    rhoarr = np.sum(rhoarr,axis=1)/(4*np.pi*rarr**2)
    if power03:
        rhobar = marr/(4*np.pi/3 * rarr**3) #10^10Msun/Mpc^3
        rhocrit = 14 #10^10Msun/Mpc^3
        p03 = np.sqrt(200)/8.0 * Narr/np.log(Narr) / np.sqrt(rhobar/rhocrit)
        p03rmin = rarr[np.min(np.where(np.logical_and(p03>=1,np.isfinite(p03)))[0])]
        return rhoarr,p03rmin
    else:
        return rhoarr
