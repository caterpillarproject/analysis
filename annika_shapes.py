import numpy as np
import scipy.linalg as linalg

def axis(arr_in, rad, shell=False, axes_out=False, fix_volume=True, quiet=False):
    """Compute axis ratios iteratively.

    WORK IN PROGRESS -- converting from gadget_profile.pro.  
    Needs to be checked.

    May want to add capability to compute axis ratios at multiple
    different radii in one function call (w/ a loop?) 

    See Dubinski & Carlberg (1991) for details

    INPUTS:
    arr_in: array of particle positions, assumed to be centered
    rad: radius to compute axis ratios.  If computing axis ratios in
    shells, rad should be a two-element list / array, with rad[0] =
    inner radius of shell and rad[1] = outer radius of shell.
    Otherwise, rad should be a real number equal to the radius within
    which to compute the axis ratios


    OPTIONAL INPUTS:
    shell=False: by default, compute cumulative axis ratio
    (i.e. for all particles within radius r).  If shell=True, compute
    axis ratio in an ellipsoidal shell instead.
    axes_out=False:  if True, also return principal axes (in ascending order)
    fix_volume=True: keep the volume of the ellipsoid constant while iterating.
    If false, keep the semi-major axis equal to the initial (spherical) search
    radius. This will result in a smaller effective volume.
    quiet=False: if set to true, suppress information that is printed
    on each iteration
    """
    def calc_inertia(arr_in, axrat_in):
	"""calculate the modified moment of inertia tensor and get its
	eigenvalues and eigenvalues"""
	tensor=np.zeros([3,3])
	# given initial values for primary axes, compute ellipsoidal
	# radius of each particle
	rp2=(arr_in[:,0]**2 +
             arr_in[:,1]**2/axrat_in[0]**2 + 
             arr_in[:,2]**2/axrat_in[1]**2)

	# construct the moment of inerial tensor to be diagonalized:
	for i in range(3):
	    for j in range(3):
		tensor[i,j]=(arr_in[:,i]*arr_in[:,j]/rp2).sum()
	
	evecs=linalg.eig(tensor)
    
	# return a tuple with eigenvalues and eigenvectors; 
	# eigenvalues=evecs[0], eigenvectors=evecs[1]
	return evecs

    cnt=0
    # initial guess for principal axes:
    evs0=np.array([[1e0,0e0,0e0],[0e0, 1e0, 0e0], [0e0, 0e0, 1e0]])
    # initial guess for axis ratios:
    axes=np.array([1e0,1e0])
    avdiff=1e0

    rad=np.asarray(rad)
    while (avdiff > 0.01) & (cnt < 100):
	axtemp=axes.copy()
	# compute ellipsoidal radius of each particle:
	dist2=(arr_in[:,0]**2 +
	       arr_in[:,1]**2/axes[0]**2 +
	       arr_in[:,2]**2/axes[1]**2)**0.5

	if shell == True:
	    # find locations of particles between rad[0] and rad[1]
            if not fix_volume:
                r_ell=rad
            else:
                r_ell=(rad**3/axes[0]/axes[1])**(1e0/3e0)
	    locs=((dist2 < r_ell[1]) & 
                  (dist2 > r_ell[0])).nonzero()[0]
	else:
	    # find locations of particles with r < rad
	    locs=(dist2 < rad).nonzero()
	    # current version keeps volume same as original spherical
	    # volume rather than keeping principal axis same as
	    # initial radius
	    # compute ellipsoidal radius = (a*b*c)^(1/3)
            if not fix_volume:
                r_ell=rad
            else:
                r_ell=(rad**3/axes[0]/axes[1])**(1e0/3e0)
	# get eigenvectors and eigenvalues.  Note: the inertia tensor
	# is a real symmetric matrix, so the eigenvalues should be real.
	axrat=calc_inertia(arr_in[locs], axtemp)
	if abs(np.imag(axrat[0])).max() > 0.:
	    #raise ValueError('Error: eigenvalues are not all real!')
            #Annika edit so that this doesn't totally crap out if I have to loop over it.
            print 'Error: eigenvalues are not all real!'
            if axes_out:
                return [-1., -1.], [evs0[:,inds[0]], evs0[:,inds[1]], evs0[:,inds[2]]]
            else:
                return [-1., -1.]
            

	evals=np.real(axrat[0])
	evecs=axrat[1]
	# get axis ratios from eigenvalues:
	axes=np.sqrt([evals[1]/evals[0], 
                      evals[2]/evals[0]])
	# rotate particles (and previous eigenvector array) into new
	# basis:
	arr_in=np.dot(arr_in, evecs)
	evs0=np.dot(evs0, evecs)
	
	# compute difference between current and previous iteration
	# for axis ratios and diagonal of eigenvector matrix
	avd0=abs((axtemp[0]-axes[0])/(axtemp[0]+axes[0]))
	avd1=abs((axtemp[1]-axes[1])/(axtemp[1]+axes[1]))

	# used to check how close most recent rotation is to zero
	# rotation (i.e. how close evecs is to the identity matrix)
	avd2=3e0-abs(evecs[0,0])-abs(evecs[1,1]) - abs(evecs[2,2])
	if quiet == False:
	    # print axis ratios relative to major axis
	    print 'axis ratios: ', (np.sort(evals/evals.max()))**0.5
            print 'deviations from previous iteration: ' + \
		'%.*e, %.*e, %.*e' % (4, avd0, 4, avd1, 4, avd2)
	    print 'number of particles in shell / sphere: ', np.size(locs)
	avdiff=max(avd0, avd1, avd2)
	cnt+=1

    # normalize eigenvalues to largest eigenvalue
    evals=evals/evals.max()
    inds=evals.argsort()
    final_axrats=evals[inds[:2]]**0.5 
    print 'number of iterations: ', cnt
    print 'number of particles in shell / sphere: ', np.size(locs), '\n'
    print 'normalized eigenvalues (sorted in ascending order) ' + \
	'and corresponding unit eigenvectors are:'
    print '%.*f' % (5, evals[inds[0]]), evs0[:,inds[0]]
    print '%.*f' % (5, evals[inds[1]]), evs0[:,inds[1]]
    print '%.*f' % (5, evals[inds[2]]), evs0[:,inds[2]]
    print 'axis ratios: ', final_axrats
    if cnt == 100:
        print 'Failed to converge after 100 iterations'
    if axes_out:
        return final_axrats, [evs0[:,inds[0]], evs0[:,inds[1]], evs0[:,inds[2]]]
    else:
        return final_axrats


# output is a, b, c, then ahat, bhat, chat. - GD a=1.0

"""Annika edit for 2d shapes"""
"""Annika Peter 12/2/2011"""
def axis2D(arr_in, rad, shell=False, axes_out=False, fix_volume=True, quiet=False):
    """Compute axis ratios iteratively.

    WORK IN PROGRESS -- converting from gadget_profile.pro.  
    Needs to be checked.

    May want to add capability to compute axis ratios at multiple
    different radii in one function call (w/ a loop?) 

    See Dubinski & Carlberg (1991) for details

    INPUTS:
    arr_in: array of particle positions, assumed to be centered
    rad: radius to compute axis ratios.  If computing axis ratios in
    shells, rad should be a two-element list / array, with rad[0] =
    inner radius of shell and rad[1] = outer radius of shell.
    Otherwise, rad should be a real number equal to the radius within
    which to compute the axis ratios


    OPTIONAL INPUTS:
    shell=False: by default, compute cumulative axis ratio
    (i.e. for all particles within radius r).  If shell=True, compute
    axis ratio in an ellipsoidal shell instead.
    axes_out=False:  if True, also return principal axes (in ascending order)
    fix_volume=True: keep the volume of the ellipsoid constant while iterating.
    If false, keep the semi-major axis equal to the initial (spherical) search
    radius. This will result in a smaller effective volume.
    quiet=False: if set to true, suppress information that is printed
    on each iteration
    """
    def calc_inertia(arr_in, axrat_in):
	"""calculate the modified moment of inertia tensor and get its
	eigenvalues and eigenvalues"""
	tensor=np.zeros([2,2])
	# given initial values for primary axes, compute ellipsoidal
	# radius of each particle
	rp2=(arr_in[:,0]**2 +
             arr_in[:,1]**2/axrat_in[0]**2)

	# construct the moment of inerial tensor to be diagonalized:
	for i in range(2):
	    for j in range(2):
		tensor[i,j]=(arr_in[:,i]*arr_in[:,j]/rp2).sum()
	
	evecs=linalg.eig(tensor)
    
	# return a tuple with eigenvalues and eigenvectors; 
	# eigenvalues=evecs[0], eigenvectors=evecs[1]
	return evecs

    cnt=0
    # initial guess for principal axes:
    evs0=np.array([[1e0,0e0],[0e0, 1e0]])
    # initial guess for axis ratios:
    axes=np.array([1e0])
    avdiff=1e0

    rad=np.asarray(rad)
    while (avdiff > 0.01) & (cnt < 100):
	axtemp=axes.copy()
	# compute ellipsoidal radius of each particle:
	dist2=(arr_in[:,0]**2 +
	       arr_in[:,1]**2/axes[0]**2)**0.5

	if shell == True:
	    # find locations of particles between rad[0] and rad[1]
            if not fix_volume:
                r_ell=rad
            else:
                r_ell=(rad**2/axes[0])**(1e0/2e0)
	    locs=((dist2 < r_ell[1]) & 
                  (dist2 > r_ell[0])).nonzero()[0]
	else:
	    # find locations of particles with r < rad
	    locs=(dist2 < rad).nonzero()
	    # current version keeps volume same as original spherical
	    # volume rather than keeping principal axis same as
	    # initial radius
	    # compute ellipsoidal radius = (a*b*c)^(1/3)
            if not fix_volume:
                r_ell=rad
            else:
                r_ell=(rad**2/axes[0])**(1e0/2e0)
	# get eigenvectors and eigenvalues.  Note: the inertia tensor
	# is a real symmetric matrix, so the eigenvalues should be real.
	axrat=calc_inertia(arr_in[locs], axtemp)
	if abs(np.imag(axrat[0])).max() > 0.:
	    #raise ValueError('Error: eigenvalues are not all real!')
            #Annika edit so that this doesn't totally crap out if I have to loop over it.
            print 'Error: eigenvalues are not all real!'
            if axes_out:
                return [-1.], [evs0[:,inds[0]], evs0[:,inds[1]]]
            else:
                return [-1.]
            

	evals=np.real(axrat[0])
	evecs=axrat[1]
	# get axis ratios from eigenvalues:
	axes=np.sqrt([evals[1]/evals[0]])
	# rotate particles (and previous eigenvector array) into new
	# basis:
	arr_in=np.dot(arr_in, evecs)
	evs0=np.dot(evs0, evecs)
	
	# compute difference between current and previous iteration
	# for axis ratios and diagonal of eigenvector matrix
	avd0=abs((axtemp[0]-axes[0])/(axtemp[0]+axes[0]))

	# used to check how close most recent rotation is to zero
	# rotation (i.e. how close evecs is to the identity matrix)
	avd2=2e0-abs(evecs[0,0])-abs(evecs[1,1])
	if quiet == False:
	    # print axis ratios relative to major axis
	    print 'axis ratios: ', (np.sort(evals/evals.max()))**0.5
            print 'deviations from previous iteration: ' + \
		'%.*e, %.*e' % (4, avd0, 4, avd2)
	    print 'number of particles in shell / sphere: ', np.size(locs)
	avdiff=max(avd0, avd2)
	cnt+=1

    # normalize eigenvalues to largest eigenvalue
    evals=evals/evals.max()
    inds=evals.argsort()
    final_axrats=evals[inds[:2]]**0.5
    print 'number of iterations: ', cnt
    print 'number of particles in shell / sphere: ', np.size(locs), '\n'
    print 'normalized eigenvalues (sorted in ascending order) ' + \
	'and corresponding unit eigenvectors are:'
    print '%.*f' % (5, evals[inds[0]]), evs0[:,inds[0]]
    print '%.*f' % (5, evals[inds[1]]), evs0[:,inds[1]]
    print 'axis ratios: ', final_axrats
    if cnt == 100:
        print 'Failed to converge after 100 iterations'
    if axes_out:
        return final_axrats, [evs0[:,inds[0]], evs0[:,inds[1]]]
    else:
        return final_axrats

