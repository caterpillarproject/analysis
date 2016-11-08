import haloutils

#hpaths = [haloutils.get_hpath_lx(hid,14) for hid in [447649, 95289, 1354437, 5320, 581141, 1130025, 1725139, 581180, 1232164, 1387186, 1292085, 94687, 1599988, 388476]]

hpaths = haloutils.get_paper_paths_lx(14)

for hpath in hpaths[16:]:
    print 'For hpath: ', hpath
    print haloutils.hpath_name(hpath)
    try:
        zoomid = haloutils.load_zoomid(hpath)
    except:
        print 'no rockstar catalogue for it'
        continue
    snap = haloutils.get_lastsnap(hpath)
    cat = haloutils.load_rscat(hpath,snap,rmaxcut=False)
    print zoomid, int(cat[0:1]['id']), 'should match'
    try:
        print cat.ix[zoomid]['mgrav'], 'zoomid mgrav'
    except:
        print 'zoomid not in catalog'
    
