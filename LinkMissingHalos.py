import MTanalysis2 as MTA
import haloutils
import numpy as np

#hid = 'H1631506'
#subRSID = 148923

hid = 'H581180'
subRSID = 166360

# to find the correct row
#hid = 'H1387186'
#subRSID = 119955

"""
Make those bar plots of blue red yellow in lx 13. is it resolution dependent?
Try linking all missing halos using below code.
Make a re-plot of original diagrams with the red bars turned blue.
Also record how many red bars turned blue, and the total number 
of red bars out of all 24 halos.

If linking is good, create a "Fixed" catalog for the halos that are actually
extent, but in the destroyed catalog. Should be able to append this to the extant data system.

use 

halo_paths = haloutils.get_paper_paths_lx(13)

if hpath:

"""


def get_link_to_row(hid, subRSID):
    hpath = haloutils.get_hpath_lx(hid,14)

    #AE = MTA.AllExtantData()
    AD = MTA.AllDestroyedData()
    TM = MTA.TagMass()

    dataD = AD.read(hpath)
    #dataE = AE.read(hpath)
    idsE,massE,idsD,massD = TM.read(hpath)
 
    snap_z0 = haloutils.get_numsnaps(hpath)-1
    cat = haloutils.load_rscat(hpath,snap_z0,rmaxcut=False)
    hostID = haloutils.load_zoomid(hpath)
    hosthalo = cat.ix[hostID]
    subs = cat.get_all_subhalos_within_halo(hostID)

    pids = cat.get_particles_from_halo(subRSID)
    inhaloD = np.in1d(pids, idsD)
    print np.sum(inhaloD)/float(len(pids)), 'fraction of particles that are stars'
    star_ids = pids[inhaloD]
    mask = dataD['infall_mvir'] > subs.ix[subRSID]['mgrav']

    for i in xrange(len(dataD[mask])):
        row = np.where(mask)[0][i]
        stars = MTA.getStars(dataD, idsD, row)
        common = np.in1d(star_ids, stars,assume_unique=True)
        n = np.sum(common)
        if n>len(star_ids)/2.:
            fixed_row = row
            print n, len(star_ids), 'number of stars in original tagging, number of stars in z=0 halo'
            print n/float(len(star_ids)), len(star_ids), 'if the first number is close to 1, and the second number is reasonably large, confidence on the match is high'
    print dataD['infall_snap'][fixed_row], 'infall snap'
    print dataD['backsnap'][fixed_row], 'back snap'
    print fixed_row, 'row in dataD'
    return fixed_row

get_link_to_row(hid, subRSID)
