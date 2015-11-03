import numpy as np
import haloutils
#import readsnapshots.readsnapHDF5_greg as rsg
from caterpillaranalysis import *
from scipy import interpolate
#from scipy.integrate import quad
import matplotlib.pyplot as plt
import grifflib as glib
from caterpillaranalysis import ProfilePlugin

def getMidpoints(bins):
    """
    Given a range of cut-off values, find the midpoints.
    @return: array of length len(bins)-1
    """
    spacing = bins[1:]-bins[:-1]
    return bins[:-1]+spacing/2.0

def distance(posA, posB,boxsize=100*1000.):
    dist = abs(posA-posB)
    tmp = dist > boxsize/2.0
    dist[tmp] = boxsize-dist[tmp]
    if dist.shape == (3,):
        return np.sqrt(np.sum(dist**2))
    else:
        return np.sqrt(np.sum(dist**2,axis=1))

# function for estimating particle boundedness.
# computes total energy assuming spherical symmetry and only
# the potential of the halo in question
# dr in megaparsecs, physical. Must not put in co-moving to this!
def PotentialE(dr,cat):
    from scipy.integrate import quad
    from scipy import interpolate
    G = 1.326*10**11 # in km^3/s^2/Msun
    mpc_to_km = 3.086*10**19

    rarr = 10**np.linspace(np.log10(min(dr))-.01, np.log10(max(dr))+.01,70)
    h_r, x_r = np.histogram(dr, bins=np.concatenate(([0],rarr)))
    m_lt_r = np.cumsum(h_r)*cat.particle_mass/cat.h0
    tck = interpolate.splrep(rarr,m_lt_r) # gives mass in Msun
    def Ufunc(x):
        return interpolate.splev(x,tck)/(x**2)

    # do it even faster by using an interpolative function
    # for computing potential energy
    # pick 60 - 100 data points
    # compute potential for all, then use an interpolation scheme
    U = np.zeros(len(rarr))
    for i in range(len(rarr)):
        r = rarr[i]
        if r > max(dr)+.05:
            print 'warning - particle outside of halo. likely inaccurate PE'
            U[i] = -G*m_lt_r[-1]/(r*mpc_to_km)
        else:
            tmp = -G*m_lt_r[-1]/(max(dr)*mpc_to_km)
            U[i] = tmp+G*quad(Ufunc,max(dr),r)[0]/mpc_to_km
    tck2 = interpolate.splrep(rarr,U)
    return interpolate.splev(dr,tck2), tck


# Find constant offset value for fitting a curve to scatter data
from fitting import *
def Fit(x_axis,y_axis,tck):
    """
    x_axis, y_axis in normal units.
    """
    def func(x, a):
        return interpolate.splev(x,tck)+a
    avar=0.0
    p0 = np.array([avar])

    popt, punc, rchi2, dof = general_fit(func, x_axis, y_axis, p0)
    print popt[0],'+-','{:.3f}'.format(punc[0]), 'offset'
    return popt[0]

def fitTriaxialNFW(hpath,a,b,c,ahat,bhat,chat):
    snap_z0 = haloutils.get_numsnaps(hpath)-1
    cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
    hostID = int(cat[0:1]['id'])
    center = np.array(cat.ix[hostID][['posX','posY','posZ']])
    particles = cat.get_all_particles_from_halo(hostID)
    particles=np.sort(particles)
    pos = haloutils.load_partblock(hpath, snap_z0, "POS ")

    pos = (pos-center)*cat.scale/cat.h0*1000 # in kpc physical
    # rotate pos2
    rotate_matrix = np.array([ahat,bhat,chat]).T
    pos2 = np.dot(pos,rotate_matrix)
    # after rotation, it should be the case that the eigenvalues
    # should align so that ahat is [1,0,0] with eigenvalue 1
    # bhat should be [0,1,0] with eigenvalue b
    # chat should be [0,0,1] with eigenvalue c

    R = np.sqrt(np.sum((pos2/np.array([a,b,c]))**2, axis=1))
    #R = np.sqrt((pos2[:,0]/a)**2+(pos2[:,1]/b)**2+(pos2[:,2]/c)**2)
    # Now plot profile based on R array

    maxr = np.max(R)
    minr = 0.1
    print maxr, minr, 'max and min r'
    binwidth = 0.04
    nbins = np.ceil((np.log10(maxr)-np.log10(minr))/binwidth)
    rarr = 10**np.linspace(np.log10(minr), np.log10(minr)+nbins*binwidth,nbins+1)
    h_r, x_r = np.histogram(R, bins=np.concatenate(([0],rarr)))
    m_lt_r = np.cumsum(h_r)*cat.particle_mass/cat.h0
    tck = interpolate.splrep(rarr,m_lt_r)
    return tck


def get_shape_values(hpath,radius=None,recalc=False,save=True):
    filename = hpath+'/'+'analysis'+'/shapevalues.dat'
    import os
    if os.path.exists(filename) and not recalc:
        values = np.fromfile(filename)
        a,b,c = values[0:3]
        ahat = values[3:6]
        bhat = values[6:9]
        chat = values[9:12]
        return [a,b,c], [ahat,bhat,chat]
    else:
        # must compute values and write to file
        import Shapes
        
        ratios, evecs = Shapes.getShape(hpath, radius=radius)
        ratios = np.append([1],ratios)
        # Now select what is a,b,c. Assume ahat = x-axis, make bhat = y-axis, chat=z-axis
        # ahat cross bhat = chat. a should always be 1.
        print ratios, evecs, 'first'
        if distance(np.cross(evecs[0],evecs[2]), evecs[1]) < 1e-3:
            # need some tolerance 
            tmp = evecs[2]
            evecs[2]=evecs[1]
            evecs[1]=tmp
            tmp = ratios[2]
            ratios[2] = ratios[1]
            ratios[1] = tmp
        print ratios,evecs, 'after'
        if save:
            f= open(filename,'wb')
            ratios.tofile(f)
            np.array(evecs).tofile(f)
            f.close()
        return ratios, evecs


# parametric potential forms, smooth, static.
# they found r_s and rho_0 effectively via least-squares minimization
# of the difference between analytic and exact potential.
# This is the best that can be done, and not what you would do to infer
# NFW parameters from measurements of a stream.
# Question is: how useful is this?


# currently computes NFW potential by matching the potentials to each other at the scale radius.
# could try to use 4*pi*G*rho_0* R_s^2, which is the correct term for a spherical NFW profile. Need rho_0 values
# For triaxial NFW, the equation for Phi isn't actually fully correct. Hard to estimate the constant in front. Could use current method, could minimize the error over entire halo

def cent(posn, pot):
    x = posn[:,0]; y = posn[:,1]; z = posn[:,2]
    m = pot #np.array([1]*len(x))
    xc=0
    yc=0
    zc=0
    fac=0.975
    nmin=1000

    r=np.sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
    rmax=r.max()
    nall=r.shape[0]
    n=nall
    while (n>nmin):
        ind=r<rmax*fac
        if (n>nmin):
            x=x[ind]
            y=y[ind]
            z=z[ind]
            m=m[ind]
            xc=(m*x).sum()/m.sum()
            yc=(m*y).sum()/m.sum()
            zc=(m*z).sum()/m.sum()
            r=np.sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
            rmax=r.max()
            n=r.shape[0]
    return np.array([xc,yc,zc])


class GravPlugin(PluginBase):
    def __init__(self):
        super(GravPlugin,self).__init__()
        self.filename='gravdata.dat'
        self.xmin=0;     self.xmax=256
        self.ymin=10**8; self.ymax=5*10**11
        self.xlog= False; self.ylog = False
        self.xlabel='' ; self.ylabel='' 
    
    def _analyze(self,hpath):
        if not haloutils.check_last_rockstar_exists(hpath):
            raise IOError("No rockstar")
        lx = haloutils.get_zoom_params(hpath)[1]
        if lx==13:
            nbins = 45
        if lx==14:
            nbins = 100

        offset_plot = False
        snap_z0 = haloutils.get_numsnaps(hpath)-1
        cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID = haloutils.load_zoomid(hpath)
        hosthalo = cat.ix[hostID]
        
        # get all particles within rvir*sqrt(3) kpc of halo center
        # this ensures a box of width rvir is fully sampled
        # also, larger radius makes offset estimate a little better
        rvir = hosthalo['rvir']/cat.h0
        posns = haloutils.load_partblock(hpath, snap_z0, "POS ")
        pot = haloutils.load_partblock(hpath, snap_z0, "POT ")
        hostpos1 = np.array(hosthalo[['posX','posY','posZ']])
        # use a different hostpos
        dist1 = distance(hostpos1,posns)*cat.scale/cat.h0*1000
        mask1 = dist1<(hosthalo['rs']/cat.h0)
        posns_tmp = posns[mask1]
        pot_tmp = pot[mask1]
        aa = np.argmin(pot_tmp)
        hostpos = np.array(posns_tmp[aa],dtype=np.float64)
        hostpos+=np.min(dist1)/20000.
        print distance(hostpos,hostpos1)*1000/cat.h0, 'offset of potential center to RS center. [kpc]'
        # back to normal code
        dist = distance(hostpos,posns)*cat.scale/cat.h0*1000 # type is float32 for some reason.
        mask = dist<(rvir*np.sqrt(3))
        dr = dist[mask]/1000.
        posns = posns[mask]
        pot = pot[mask]
        posns = (posns-hostpos)*cat.scale/cat.h0*1000

        # Find the offset
        U, tck_mass = PotentialE(dr,cat)
        # write this to file if possible

        argsorted = np.argsort(U)
        tck = interpolate.splrep(np.arange(len(U)), U[argsorted])
        offset = -Fit(np.arange(len(U)), pot[argsorted], tck)
        if offset_plot==True:
            plt.plot(np.arange(len(U)), U[argsorted])
            plt.scatter(np.arange(len(U)), pot[argsorted]+offset)
            plt.ylabel('Potential Energy')
            plt.xlabel('boundedness order')
            plt.show()

        # Get triaxial NFW fit
        r_s = hosthalo['rs']/cat.h0
        #ratios, evecs = get_shape_values(hpath,radius=r_s,recalc=False)
        #a,b,c = ratios; ahat, bhat, chat = evecs
        # the above two lines will get triaxial fit, but they seem to be failing
        # to be accurate
        a=b=c=1.0
        ahat=[1,0,0]; bhat=[0,1,0]; chat=[0,0,1]

        # rotate positions, use square mask on posns now
        rotate_matrix = np.array([ahat,bhat,chat]).T
        posns_r = np.dot(posns,rotate_matrix)

        # plot a x-y plane projection now
        zcut = np.abs(posns_r[:,2])<5
        pos= posns_r[zcut]
        po = pot[zcut]+offset
        po = np.abs(po)
        counts,xedges,yedges=np.histogram2d(pos[:,0],pos[:,1],bins=nbins)
        total,xedges2,yedges2=np.histogram2d(pos[:,0],pos[:,1],bins=nbins,normed=False,weights=po)
        matrix = total/counts

        # now again, plot the grav potential based on this
        grid = np.meshgrid(getMidpoints(xedges), getMidpoints(yedges))
        xdata = grid[0]
        ydata = grid[1][::-1]
        R = np.sqrt((xdata/a)**2+(ydata/b)**2)

        r_s = hosthalo['rs']/cat.h0
        R_s = np.sqrt(r_s**2/3.*(1+1/b**2+1/c**2))
        G = 1.326*10**11 # in km^3/s^2/Msun
        kpc_to_km = 3.086*10**16
        V_s = np.sqrt(G*interpolate.splev(r_s/1000.,tck_mass)/(r_s*kpc_to_km)) # this doesn't work so well

        # estimate phi at r_s from gadget POT data
        mask1 = dr<(r_s/1000*1.05)
        mask2 = dr>(r_s/1000*.95)
        mean_phi_rs = np.mean(pot[mask1*mask2])+offset
        v_s = np.sqrt(np.abs(mean_phi_rs)/np.log(2))

        # estimate v_s from rho_0
        profile= ProfilePlugin()
        r,mltr,p03r,rvir,r200c,pNFW,pEIN = profile.read(hpath)
        rho_0 = pNFW[1] # in Msun/kpc^3
        G = 4.52*10**-39 # in kpc^3/s^2/Msun
        v_s2 = np.sqrt(4*np.pi*G*rho_0*R_s**2*kpc_to_km**2) # in km^2/s^2
        
        print v_s, v_s2, 'constant from data, NFW fit'

        def compute_phi(R,R_s,V_s):
            return V_s**2 * R_s/R*np.log(1+R/R_s) # returns positive values

        def forceAspect(ax,aspect=1):
            im = ax.get_images()
            extent =  im[0].get_extent()
            ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
        smoothPhi = compute_phi(R,R_s,v_s)
        # write out smoothPhi, matrix, xedges, yedges
        g = open(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename,'wb')
        np.array(xedges).tofile(g)
        np.array(yedges).tofile(g)
        np.array(matrix).tofile(g)
        np.array(smoothPhi).tofile(g)
        g.close()

    def _read(self,hpath):
        lx = haloutils.get_zoom_params(hpath)[1]
        if lx==13:
            nbins = 45
        if lx==14:
            nbins = 100

        data = np.fromfile(hpath+'/'+self.OUTPUTFOLDERNAME+'/'+self.filename)
        xedges = data[0:nbins+1]
        yedges = data[nbins+1:(nbins+1)*2]
        matrix = data[(nbins+1)*2:(nbins+1)*2+nbins**2]
        smoothPhi = data[(nbins+1)*2+nbins**2:(nbins+1)*2+2*nbins**2]
        return xedges,yedges,np.reshape(matrix,(nbins,nbins)),np.reshape(smoothPhi,(nbins,nbins))


## This code segment used to run analyze
lx = 14
figdir = '/bigbang/data/AnnaGroup/GregFigs/lx'+str(lx)+'/'
halo_paths = haloutils.find_halo_paths(levellist=[lx],require_mergertree=True,require_subfind=False,verbose=False) 
offset_plot = False
GP = GravPlugin()

for hpath,i  in zip(halo_paths, range(len(halo_paths))):
    GP.analyze(hpath,recalc=True)
    print 'done with', i, hpath


# uncomment to make plots
"""
for hpath in halo_paths:
    xedges,yedges,matrix,smoothPhi = GP.read(hpath)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    fig, (ax1,ax2) = plt.subplots(1,2, figsize=(15,4.5),sharey=False)
    vmin= 125 ; vmax = 310; limit = 100 #limit = rvir
    vmax = np.nanmax(np.sqrt(matrix))+10
    xb = np.where(np.abs(xedges)<100)[0]; yb=np.where(np.abs(yedges)<100)[0]
    vmin = np.nanmin(np.sqrt(matrix[xb[0]:xb[-1]+1,yb[0]:yb[-1]+1]))-10   # min of potential in the -100 to 100 range

    im1=ax1.imshow(np.sqrt(matrix), extent=extent,vmin=vmin,vmax=vmax,aspect='equal',cmap='cubehelix')
    ax1.set_xlabel('x [kpc]')
    ax1.set_ylabel('y [kpc]')
    ax1.set_ylim((-limit,limit))
    ax1.set_xlim((-limit,limit))
    ax1.text(-limit*.9,limit*.8, 'True Data $\sqrt{|\phi|}$ (km/s)',color='white',fontsize=14)
    
    im2=ax2.imshow(np.sqrt(smoothPhi), extent=extent,vmin=vmin,vmax=vmax,aspect='equal',cmap='cubehelix')
    ax2.set_xlabel('x [kpc]')
    ax2.set_ylim((-limit,limit))
    ax2.set_xlim((-limit,limit))
    ax2.text(-limit*.9,limit*.8, 'NFW Fit $\sqrt{|\phi|}$ (km/s)',color='white',fontsize=14)
    plt.setp(ax2.get_yticklabels(), visible=False)
    fig.tight_layout()

    fig.subplots_adjust(right=0.605,left=0.02)
    cbar_ax = fig.add_axes([.59,.10,0.02,0.87])
    fig.colorbar(im1,cax=cbar_ax)

    percent_diff = 100*(np.sqrt(smoothPhi)-np.sqrt(matrix))/np.sqrt(matrix)
    ax3 = fig.add_axes([.67,.10,.28,.87])
    im3 = ax3.imshow(percent_diff,extent=extent,aspect='equal',cmap='cubehelix',vmin=-10,vmax=10)
    ax3.set_xlabel('x [kpc]')
    ax3.set_ylim((-limit,limit))
    ax3.set_xlim((-limit,limit))
    ax3.text(-limit*.9,limit*.8, 'Residual (%)',color='black',fontsize=14)
    cbar_ax2 = fig.add_axes([.96,.10,.02,.87])
    fig.colorbar(im3,cax=cbar_ax2)

    # for plotting circles
    snap_z0 = haloutils.get_numsnaps(hpath)-1
    cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)    
    hostID = haloutils.load_zoomid(hpath)
    hosthalo = cat.ix[hostID]
    hostpos = np.array(hosthalo[['posX','posY','posZ']])
    subs = cat.get_subhalos_within_halo(hostID)
    maskv = subs['vmax']> 25
    maskz = np.abs((np.array(subs['posZ'])-hostpos[2])*cat.scale*1000/cat.h0) < 10
    maskv = maskv*maskz
    x_arr=(np.array(subs[maskv]['posX'])-hostpos[0])*cat.scale*1000/cat.h0
    y_arr=(np.array(subs[maskv]['posY'])-hostpos[1])*cat.scale*1000/cat.h0
    r_arr=np.array(subs[maskv]['rvir']/cat.h0) # in kpc already
    for x,y,rvir in zip(x_arr,y_arr,r_arr):
        xcirc,ycirc = glib.drawcircle(x,y,rvir)
        ax1.plot(xcirc,ycirc,'k-')
        ax3.plot(xcirc,ycirc,'k-')


    plt.savefig(figdir+haloutils.hpath_name(hpath)+'_grav_potential', bbox_inches='tight')
    print 'done with', haloutils.hpath_name(hpath)
"""







"""
for hpath in halo_paths:
    #hpath = halo_paths[0]
    snap_z0 = haloutils.get_numsnaps(hpath)-1
    cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
    hostID = haloutils.load_zoomid(hpath)
    hosthalo = cat.ix[hostID]

    # get all particles within rvir*sqrt(3) kpc of halo center
    # this ensures a box of width rvir is fully sampled
    # also, larger radius makes offset estimate a little better
    rvir = hosthalo['rvir']/cat.h0
    posns = haloutils.load_partblock(hpath, snap_z0, "POS ")
    pot = haloutils.load_partblock(hpath, snap_z0, "POT ")
    hostpos = np.array(hosthalo[['posX','posY','posZ']])
    dist = distance(hostpos,posns)*cat.scale/cat.h0*1000
    mask = dist<(rvir*np.sqrt(3))
    dr = dist[mask]/1000.
    posns = posns[mask]
    pot = pot[mask]
    posns = (posns-hostpos)*cat.scale/cat.h0*1000

    # Find the offset
    U, tck_mass = PotentialE(dr,cat)
    # write this to file if possible

    argsorted = np.argsort(U)
    tck = interpolate.splrep(np.arange(len(U)), U[argsorted])
    offset = -Fit(np.arange(len(U)), pot[argsorted], tck)
    if offset_plot==True:
        plt.plot(np.arange(len(U)), U[argsorted])
        plt.scatter(np.arange(len(U)), pot[argsorted]+offset)
        plt.ylabel('Potential Energy')
        plt.xlabel('boundedness order')
        plt.show()

    # Get triaxial NFW fit
    r_s = hosthalo['rs']/cat.h0
    #ratios, evecs = get_shape_values(hpath,radius=r_s,recalc=False)
    #a,b,c = ratios; ahat, bhat, chat = evecs
    # the above two lines will get triaxial fit, but they seem to be failing
    # to be accurate
    a=b=c=1.0
    ahat=[1,0,0]; bhat=[0,1,0]; chat=[0,0,1]


    # rotate positions, use square mask on posns now
    rotate_matrix = np.array([ahat,bhat,chat]).T
    posns_r = np.dot(posns,rotate_matrix)


    # plot a x-y plane projection now
    zcut = np.abs(posns_r[:,2])<5
    pos= posns_r[zcut]
    po = pot[zcut]+offset
    po = np.abs(po)
    counts,xedges,yedges=np.histogram2d(pos[:,0],pos[:,1],bins=45)
    total,xedges2,yedges2=np.histogram2d(pos[:,0],pos[:,1],bins=45,normed=False,weights=po)
    matrix = total/counts

    # now again, plot the grav potential based on this
    grid = np.meshgrid(getMidpoints(xedges), getMidpoints(yedges))
    xdata = grid[0]
    ydata = grid[1][::-1]
    R = np.sqrt((xdata/a)**2+(ydata/b)**2)

    r_s = hosthalo['rs']/cat.h0
    R_s = np.sqrt(r_s**2/3.*(1+1/b**2+1/c**2))
    G = 1.326*10**11 # in km^3/s^2/Msun
    kpc_to_km = 3.086*10**16
    V_s = np.sqrt(G*interpolate.splev(r_s/1000.,tck_mass)/(r_s*kpc_to_km)) # this doesn't work so well

    # estimate phi at r_s from gadget POT data
    mask1 = dr<(r_s/1000*1.05)
    mask2 = dr>(r_s/1000*.95)
    mean_phi_rs = np.mean(pot[mask1*mask2])+offset
    v_s2 = np.sqrt(np.abs(mean_phi_rs)/np.log(2))

    # estimate v_s from rho_0
    profile= ProfilePlugin()
    r,mltr,p03r,rvir,r200c,pNFW,pEIN = profile.read(hpath)
    rho_0 = pNFW[1] # in Msun/kpc^3
    G = 4.52*10**-39 # in kpc^3/s^2/Msun
    v_s = np.sqrt(4*np.pi*G*rho_0*R_s**2*kpc_to_km**2) # in km^2/s^2

    print v_s, v_s2, 'constant from NFW fit, and from data'

    def compute_phi(R,R_s, V_s):
        return V_s**2 * R_s/R*np.log(1+R/R_s) # returns positive values

    def forceAspect(ax,aspect=1):
        im = ax.get_images()
        extent =  im[0].get_extent()
        ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

    smoothPhi = compute_phi(R,R_s,v_s)
"""





#im3 = plt.imshow(percent_diff,extent=extent,aspect='equal')
#plt.ylim((-limit,limit))
#plt.xlim((-limit,limit))
#plt.colorbar()
#plt.show()

"""
from matplotlib import gridspec
fig = plt.figure()
gs = gridspec.GridSpec(1,3, width_ratios=[9,9,1],height_ratios=[1,1,1])
ax1 = plt.subplot(gs[:,0])
ax2 = plt.subplot(gs[:,1])
ax3 = plt.subplot(gs[:,2])

fig.set_figheight(6.5)
fig.set_figwidth(fig.get_figheight()*(2+1/(19.)))

extent = [-200,200,-200,200]
vmin= 100 ; vmax = 300
im1=ax1.imshow(np.sqrt(matrix), extent=extent,vmin=vmin,vmax=vmax,aspect='equal')
ax1.set_xlabel('x [kpc]')
ax1.set_ylabel('y [kpc]')
#ax1.set_ylim((-200,200))
#ax1.set_xlim((-200,200))

im2=ax2.imshow(np.sqrt(smoothPhi), extent=extent,vmin=vmin,vmax=vmax,aspect='equal')
ax2.set_xlabel('x [kpc]')
#ax2.set_ylim((-200,200))
#ax2.set_xlim((-200,200))
cbar1 = plt.colorbar(im1,cax=ax3,orientation='vertical')

gs.tight_layout(fig)
print fig.get_figwidth()
print fig.get_figheight()
plt.savefig(figdir+'smooth')
plt.show()
"""







