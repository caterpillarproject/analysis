import numpy as np
import haloutils
import time,sys
from numpy.lib.recfunctions import append_fields
from caterpillaranalysis import PluginBase
from scipy import interpolate

class MinihaloFinderPlugin(PluginBase):
    def __init__(self,verbose=False,Tvir=2000):
        super(MinihaloFinderPlugin,self).__init__()
        self.filename = 'minihalo_array.npy'
        if Tvir != 2000:
            self.filename = 'minihalo_array_{}.npy'.format(int(Tvir))
        self.verbose = verbose
        self.Tvir = Tvir
        self.use_all_trees = False

    def mcrit(self,T,z):
        h = 0.6711
        omega_m = 0.3125
        M = 1e8/h * (T/(1.+z))**1.5 * (0.6*10/1.22/1.98e4)**1.5 * (18*3.14*3.14/178/omega_m)**0.5 #in solar masses
        return M

    def is_above_threshold(self,scale,mass):
        z = 1./scale-1.0
        mcrit = self.mcrit(self.Tvir,z)
        return mass > mcrit
    def is_above_threshold_and_not_phantom(self,scale,mass,phantom):
        threshold = self.is_above_threshold(scale,mass)
        not_phantom = phantom == 0
        return threshold * not_phantom
    
    def search_one_tree(self,mt,h0=0.6711):
        start = time.time()
        # This step is really slow but not bad enough to warrant fixing yet
        mt.data = append_fields(mt.data,['base_rsid','row'],data=[np.zeros(len(mt.data),dtype=int)+mt.rockstar_id,np.arange(len(mt.data),dtype=int)])
        #if self.verbose: print "Time to append: {:.1f}".format(time.time()-start)
    
        threshold = self.is_above_threshold_and_not_phantom(mt['scale'],mt['mvir']/h0,mt['phantom'])
        minihalo_candidates = mt.data[threshold]
        is_minihalo = np.zeros(len(minihalo_candidates),dtype=bool)
        for i,mh in enumerate(minihalo_candidates):
            dfid_base = mh['dfid']; dfid_last = mh['lastprog_dfid']
            mask = np.logical_and(mt.data['dfid']>=dfid_base, mt.data['dfid']<=dfid_last)
            if np.sum(threshold[mask]) == 1:
                is_minihalo[i] = True
        minihalos = minihalo_candidates[is_minihalo]
        if self.verbose: print "Time to find {} minihalos in tree {}: {:.1f}".format(len(minihalos),mt.rockstar_id,time.time()-start)
        return minihalos

    def _analyze(self,hpath):
        import sys
        if self.use_all_trees:
            print 'loading mtc'
            sys.stdout.flush()
            mtc = haloutils.load_mtc(hpath,indexbyrsid=True)
            print 'loaded all of mtc'
            sys.stdout.flush()
        else:
            mtc = haloutils.load_zoom_mtc(hpath,indexbyrsid=True)
        all_minihalos = []
        start = time.time()
        for base_rsid,mt in mtc.Trees.iteritems():
            all_minihalos.append(self.search_one_tree(mt))
            print 'finished a tree'
            sys.stdout.flush()
        print "Total time: {:.1f}".format(time.time()-start)
        all_minihalos = np.concatenate(all_minihalos)
        #assert not all_minihalos.masked,"If masked there is probably some data corruption!"
        all_minihalos = np.array(all_minihalos)
        np.save(self.get_outfname(hpath),all_minihalos)
    
    def _read(self,hpath):
        return np.load(self.get_outfname(hpath))

class AllMinihaloFinderPlugin(MinihaloFinderPlugin):
    """
    Identical to MinihaloFinderPlugin,
    except uses load_mtc instead of load_zoom_mtc
    (and has a different output file)
    """
    def __init__(self,verbose=False,Tvir=2000):
        super(AllMinihaloFinderPlugin,self).__init__(verbose=verbose,Tvir=Tvir)
        self.filename = 'all_minihalo_array.npy'
        if Tvir != 2000:
            self.filename = 'all_minihalo_array_{}.npy'.format(int(Tvir))
        self.verbose = verbose
        self.Tvir = Tvir
        self.use_all_trees = True
        print "initialized AllMinihaloFinderPlugin"

class AllFirstGalFinderPlugin(MinihaloFinderPlugin):
    """
    Identical to MinihaloFinderPlugin,
    except uses load_mtc instead of load_zoom_mtc,
    has a 10^4K cutoff by default,
    and has a different output file
    """
    def __init__(self,verbose=False,Tvir=10000.):
        super(AllFirstGalFinderPlugin,self).__init__(Tvir=Tvir)
        self.filename = 'all_firstgal_array.npy'
        if Tvir != 10000:
            self.filename = 'all_firstgal_array_{}.npy'.format(int(Tvir))
        self.verbose = verbose
        self.Tvir = Tvir
        self.use_all_trees = True

class ModifiedRicottiMinihaloFinderPlugin(MinihaloFinderPlugin):
    """
    Use a constant minimum mass until it intersects the Ricotti LW model
    """
    def __init__(self,logMmin,verbose=False,Tvir=2000):
        super(ModifiedRicottiMinihaloFinderPlugin,self).__init__(verbose=verbose,Tvir=Tvir)
        assert logMmin in [6.0, 6.5, 7.0, 7.5, 8.0]
        self.logMmin = logMmin
        self.lw_data_dir = '/bigbang/data/bgriffen/crosby2013'
        
        ## Fit Ricotti
        zr,Mcrit = np.loadtxt(self.lw_data_dir+'/ricotti.txt',delimiter=',', unpack=True)
        # Add some things to extend interpolation to low and high z (just flat)
        zr = np.array([127]+list(zr)+[0])
        Mcrit = np.array([Mcrit[0]]+list(Mcrit)+[Mcrit[-1]])
        self.logffit = interpolate.interp1d(zr,np.log10(Mcrit))
        self.ffit = lambda z: 10**self.logffit(z)
        
        self.filename = "modified_ricotti_minihalo_array_{:.1f}.npy".format(logMmin)
        self.verbose = verbose
        self.use_all_trees = False
        
    def is_above_threshold(self, scale, mass):
        ## The MT returns masked arrays because of asciitable, but nothing should be masked.
        ## interp1d does not work if it is a masked array
        ## Here I simply np.ma.filled to get rid of it.
        z = 1./np.ma.filled(scale)-1.0
        ## Just to be explicit: don't use mcrit anywhere for this one!
        #mcrit = self.mcrit(self.Tvir,z)
        mricotti = self.ffit(z)
        return (np.log10(mass) > self.logMmin) & (mass > mricotti)

class LWMinihaloFinderPlugin(MinihaloFinderPlugin):
    def __init__(self,verbose=False,lwimf='kroupa',Tvir=2000):
        super(LWMinihaloFinderPlugin,self).__init__(verbose=verbose,Tvir=Tvir)
        assert lwimf in ['kroupa','chabrier','salpeter']
        print 'using the lw imf model', lwimf
        self.lwimf = lwimf
        self.lw_data_dir = '/bigbang/data/bgriffen/crosby2013'
        ffit,x = self.get_imf_func(self.lwimf)
        self.ffit = ffit
        self.mlwmin = 7e4

        self.Tvir = Tvir
        if Tvir != 2000:
            self.filename = 'lw_minihalo_array_{}_{}.npy'.format(self.lwimf,self.Tvir)
        else:
            self.filename = 'lw_minihalo_array_'+self.lwimf+'.npy'
        self.verbose = verbose
        self.use_all_trees = False

    def get_poly(self,data):
        import numpy.polynomial.polynomial as poly
        x = data[:,0]
        y = data[:,1]
        mask = (x < 26)
        coefs = poly.polyfit(x[mask], y[mask], 4)
        x_new = np.linspace(x[mask][0], x[mask][-1], num=len(x[mask])*11)
        ffit = poly.Polynomial(coefs)
        return ffit,x_new
    def get_imf_func(self,imf):
        data = np.loadtxt(self.lw_data_dir+'/'+imf+'.txt',delimiter=',')
        ffit,x = self.get_poly(data)
        return ffit,x
    def is_above_threshold(self,scale,mass):
        z = 1./scale-1.0
        mcrit = self.mcrit(self.Tvir,z)
        mlw = self.ffit(z)
        return (mass > mcrit) & (mass > self.mlwmin) & (mass > mlw)

class AllLWMinihaloFinderPlugin(LWMinihaloFinderPlugin):
    def __init__(self,verbose=False,lwimf='kroupa',Tvir=2000):
        super(AllLWMinihaloFinderPlugin,self).__init__(verbose=verbose,lwimf=lwimf,Tvir=Tvir)
        self.filename = 'all_'+self.filename
        self.use_all_trees = True

if __name__=="__main__":
    if len(sys.argv) == 3:
        hid = int(sys.argv[1])
        lx = int(sys.argv[2])
        assert lx==14 or lx==15
        plug = MinihaloFinderPlugin(verbose=True)
        hpath = haloutils.get_hpath_lx(hid,lx)
        MHs = plug.read(hpath,recalc=True)
    elif len(sys.argv) == 4:
        if sys.argv[3] == "All":
            hid = int(sys.argv[1])
            lx = int(sys.argv[2])
            assert lx==14 or lx==15
            plug = AllMinihaloFinderPlugin(verbose=True)
            hpath = haloutils.get_hpath_lx(hid,lx)
            if lx==15 and hid ==1387186:
                hpath = "/bigbang/data/AnnaGroup/caterpillar/halos/H1387186/H1387186_EB_Z127_P7_LN7_LX15_O4_NV4"
            print 'loaded hpath in All Mini halos'
            MHs = plug.read(hpath,recalc=True) #read
        elif sys.argv[3] == "AllFirstGal":
            hid = int(sys.argv[1])
            lx = int(sys.argv[2])
            assert lx==14 or lx ==15
            plug = AllFirstGalFinderPlugin(verbose=True)
            hpath = haloutils.get_hpath_lx(hid,lx)
            if lx==15 and hid ==1387186:
                hpath = "/bigbang/data/AnnaGroup/caterpillar/halos/H1387186/H1387186_EB_Z127_P7_LN7_LX15_O4_NV4"

            MHs = plug.read(hpath,recalc=True)
        elif sys.argv[3] == "LW":
            hid = int(sys.argv[1])
            lx = int(sys.argv[2])
            assert lx==14 or lx ==15
            plug = LWMinihaloFinderPlugin(verbose=True)
            hpath = haloutils.get_hpath_lx(hid,lx)
            if lx==15 and hid ==1387186:
                hpath = "/bigbang/data/AnnaGroup/caterpillar/halos/H1387186/H1387186_EB_Z127_P7_LN7_LX15_O4_NV4"

            MHs = plug.read(hpath,recalc=True)
        elif sys.argv[3] == "LWsalpeter":
            hid = int(sys.argv[1])
            lx = int(sys.argv[2])
            assert lx==14 or lx==15
            plug = LWMinihaloFinderPlugin(verbose=True,lwimf='salpeter')
            hpath = haloutils.get_hpath_lx(hid,lx)
            if lx==15 and hid ==1387186:
                hpath = "/bigbang/data/AnnaGroup/caterpillar/halos/H1387186/H1387186_EB_Z127_P7_LN7_LX15_O4_NV4"

            MHs = plug.read(hpath,recalc=True)
        elif sys.argv[3] == "AllLW":
            hid = int(sys.argv[1])
            lx = int(sys.argv[2])
            assert lx==14 or lx==15
            plug = AllLWMinihaloFinderPlugin(verbose=True)
            hpath = haloutils.get_hpath_lx(hid,lx)
            if lx==15 and hid ==1387186:
                hpath = "/bigbang/data/AnnaGroup/caterpillar/halos/H1387186/H1387186_EB_Z127_P7_LN7_LX15_O4_NV4"

            MHs = plug.read(hpath,recalc=True)
        elif sys.argv[3] == "AllLWchabrier":
            hid = int(sys.argv[1])
            lx = int(sys.argv[2])
            assert lx==14 or lx == 15
            plug = AllLWMinihaloFinderPlugin(verbose=True,lwimf='chabrier')
            hpath = haloutils.get_hpath_lx(hid,lx)
            if lx==15 and hid ==1387186:
                hpath = "/bigbang/data/AnnaGroup/caterpillar/halos/H1387186/H1387186_EB_Z127_P7_LN7_LX15_O4_NV4"

            MHs = plug.read(hpath,recalc=True)
        elif sys.argv[3] == "AllLWsalpeter":
            hid = int(sys.argv[1])
            lx = int(sys.argv[2])
            assert lx==14 or lx==15
            plug = AllLWMinihaloFinderPlugin(verbose=True,lwimf='salpeter')
            hpath = haloutils.get_hpath_lx(hid,lx)
            if lx==15 and hid ==1387186:
                hpath = "/bigbang/data/AnnaGroup/caterpillar/halos/H1387186/H1387186_EB_Z127_P7_LN7_LX15_O4_NV4"

            MHs = plug.read(hpath,recalc=True)
        elif sys.argv[3] == "AllLWkroupa":
            hid = int(sys.argv[1])
            lx = int(sys.argv[2])
            assert lx==14 or lx==15
            plug = AllLWMinihaloFinderPlugin(verbose=True,lwimf='kroupa')
            hpath = haloutils.get_hpath_lx(hid,lx)
            if lx==15 and hid ==1387186:
                hpath = "/bigbang/data/AnnaGroup/caterpillar/halos/H1387186/H1387186_EB_Z127_P7_LN7_LX15_O4_NV4"

            MHs = plug.read(hpath,recalc=True)

        elif sys.argv[3].startswith("RF"):
            logMmin = float(sys.argv[3][2:])
            assert logMmin in [6.0, 6.5, 7.0, 7.5, 8.0]
            hid = int(sys.argv[1])
            lx = int(sys.argv[2])
            assert lx==14 or lx==15
            plug = ModifiedRicottiMinihaloFinderPlugin(logMmin,verbose=True)
            hpath = haloutils.get_hpath_lx(hid,lx)
            if lx==15 and hid ==1387186:
                hpath = "/bigbang/data/AnnaGroup/caterpillar/halos/H1387186/H1387186_EB_Z127_P7_LN7_LX15_O4_NV4"
            MHs = plug.read(hpath,recalc=True)

        else:
            raise ValueError("Only accepts \"LW\", \"AllLW\", \"All\", \"LWsalpeter\", \"AllFirstGal\", or \"RFXX\"")
