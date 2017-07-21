import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import haloutils as htils
import DwarfMethods as dm
import abundance_matching as am


"""
hpaths = dm.get_hpaths(field=False, lx=14)
vmax_ach = 9.48535156 
vmax_filt = 23.54248047
frac_to_tag = .02
print frac_to_tag, 'frac tagged'
frac_all = []
dist_all=[]
for hpath in hpaths:
    snap_z0 = htils.get_numsnaps(hpath)-1
    cat_z0 = htils.load_rscat(hpath,snap_z0,rmaxcut=False)
    hostID = htils.load_zoomid(hpath)
    dataE = dm.get_extant_data(hpath,field=False)

    # now only get rsids of halos that survive reionization #        
    mask_pre = np.array(dataE['vmax_9'] > vmax_ach)   # use redshift 9 reionization
    mask_post = np.array(dataE['peak_vmax'] >= vmax_filt)
    dataE = dataE[mask_pre & mask_post]

    frac_left = []; dists = []
    # want it to be the same ordering
    for line in dataE.index:
        print line
        snap = int(dataE.ix[line]['max_mass_snap'])
        cat = htils.load_rscat(hpath,snap,rmaxcut=False)
        maxmass_rsid = int(dataE.ix[line]['max_mass_rsid'])
        iPids = cat.get_all_particles_from_halo(maxmass_rsid)
        star_pids = iPids[0:int(np.round(len(iPids)*frac_to_tag))]
    
        rsid_z0 = dataE.ix[line]['rsid']
        zPids = cat_z0.get_all_particles_from_halo(rsid_z0)
        frac = np.sum(np.in1d(star_pids, zPids, assume_unique=True)) / float(len(star_pids))
        frac_left.append(frac)
        dist = dm.distance(np.array(cat_z0.ix[rsid_z0][['posX','posY','posZ']]), np.array(cat_z0.ix[hostID][['posX','posY','posZ']]))*1000/cat_z0.h0
        dists.append(dist)
        print frac, dist
    print frac_left, 'frac left of cat', htils.hpath_catnum(hpath)
    frac_all.append(frac_left)
    dist_all.append(dists)

np.save('frac_remaining_all2', np.array(frac_all))
np.save('dist_all2', np.array(dist_all))

"""





# plot the fraction stripped vs peak mass
hpaths = dm.get_hpaths(field=False, lx=14)
vmax_ach = 9.48535156 
vmax_filt = 23.54248047

peak_mass = []
dists = []
for hpath in hpaths:
    snap_z0 = htils.get_numsnaps(hpath)-1
    cat_z0 = htils.load_rscat(hpath,snap_z0,rmaxcut=False)
    hostID = htils.load_zoomid(hpath)
    dataE = dm.get_extant_data(hpath,field=False)

    # now only get rsids of halos that survive reionization #        
    mask_pre = np.array(dataE['vmax_9'] > vmax_ach)
    mask_post = np.array(dataE['peak_vmax'] >= vmax_filt)
    dataE = dataE[mask_pre & mask_post]

    peak_mass = np.append(peak_mass, np.array(dataE['max_mass']/cat_z0.h0))
    rsids = dataE['rsid']
    dist = dm.distance(np.array(cat_z0.ix[rsids][['posX','posY','posZ']]), np.array(cat_z0.ix[hostID][['posX','posY','posZ']]))*1000/cat_z0.h0
    dists = np.append(dists, dist)
    print htils.hpath_catnum(hpath)

peak_mass = np.array(peak_mass)
dists = np.array(dists)

#frac_left = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/frac_remaining_all5.npy')
#dists1 = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/dist_all5.npy')

frac_left = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/frac_remaining_all2.npy')
dists1 = np.load('/bigbang/data/AnnaGroup/greg_backup/Dropbox/DwarfsOfDwarfs/code/dist_all2.npy')


frac_left_flat = np.array([item for row in frac_left for item in row])
dists_flat =  np.array([item for row in dists1 for item in row])



nbins = 6
bins = np.logspace(8.3,11,nbins+1)
mean_peak=[]; mean_frac = []
frac_scatter=[]
for i in range(nbins):
    mask = (peak_mass < bins[i+1])&(peak_mass > bins[i])
    mean_peak.append(np.median(peak_mass[mask]))
    mean_frac.append(np.mean(frac_left_flat[mask]))
    frac_scatter.append(np.std(frac_left_flat[mask]))
mean_frac=np.array(mean_frac); frac_scatter = np.array(frac_scatter)
#plt.scatter(peak_mass, frac_left_flat)
plt.plot(mean_peak, mean_frac)
plt.fill_between(mean_peak, mean_frac+frac_scatter, mean_frac-frac_scatter, alpha=0.25)
plt.xscale('log')
plt.ylabel('Star Fraction Remaining')
plt.xlabel('Peak Mass')
plt.ylim((0,1))
plt.savefig('LMCPlots/strip_v_peakmass2%')
plt.close()



nbins = 12
bins = np.linspace(0,300,nbins+1)
mean_dist=[]; mean_frac = []
frac_scatter=[]
for i in range(nbins):
    mask = (dists_flat < bins[i+1])&(dists_flat > bins[i])&(peak_mass > 3*10**9)
    mean_dist.append(np.median(dists_flat[mask]))
    mean_frac.append(np.mean(frac_left_flat[mask]))
    frac_scatter.append(np.std(frac_left_flat[mask]))

mean_dist2=[]; mean_frac2 = []
frac_scatter2=[]
for i in range(nbins):
    mask = (dists_flat < bins[i+1])&(dists_flat > bins[i])&(peak_mass < 10**9)
    mean_dist2.append(np.median(dists_flat[mask]))
    mean_frac2.append(np.mean(frac_left_flat[mask]))
    frac_scatter2.append(np.std(frac_left_flat[mask]))


mean_frac=np.array(mean_frac); frac_scatter = np.array(frac_scatter)
#plt.scatter(dists_flat, frac_left_flat)
plt.plot(mean_dist, mean_frac)
plt.fill_between(mean_dist, mean_frac+frac_scatter, mean_frac-frac_scatter, alpha=0.25)
plt.plot(mean_dist2, mean_frac2)
plt.xlabel('Distance [kpc]')
plt.ylabel('Star Fraction Remaining')
#plt.ylim((0,.3))
plt.xlim((0,310))
plt.grid()
plt.savefig('LMCPlots/strip_v_dist2%')
plt.close()


"""

## Want to compare the fraction stripped of subhalos like the largest 11 in the MW
## to the subhalos within 100 kpc of MW, 50 kpc of LMC, and ~50 kpc of SMC

galx = dm.load_nearby_gx()
mass_cut  = galx['mstar']>2.2*10**5
argsort = np.argsort(galx[mass_cut]['mstar'])[::-1]

galx[mass_cut][argsort]
galx[mass_cut]['dist_GC']
np.mean(galx[mass_cut]['dist_GC'][1:])  # exclude Canis Major which is not in Shea's list

near_LMC = galx['dist_LMC'] < 50
masssort = np.argsort(galx[near_LMC]['mstar'])
galx[near_LMC][masssort]
galx[near_LMC]['dist_GC'][masssort]
np.mean(galx[near_LMC]['dist_GC'])
"""
