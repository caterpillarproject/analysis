import numpy as np
import pylab as plt
import os,subprocess,sys,time
import haloutils
import csv

from caterpillaranalysis import PluginBase,MassAccrPlugin
import brendanlib.conversions as bconversions
from scipy import optimize

class FormationTimePlugin(PluginBase):
    def __init__(self):
        super(FormationTimePlugin,self).__init__()
        self.filename='formationtime.csv'
        self.mbplug = MassAccrPlugin()
        self.qtynames = ['a_lmm','a_half','conc','a_exp','a_betagam','a_gmb','a_mmm']

        self.xmin = 0; self.xmax = 1
        self.ymin = 10**6;  self.ymax = 10**13
        self.xlabel = r'$\rm{scale\ factor}$'
        self.ylabel = r'$M\ (M_\odot)$'
        self.n_xmin = 0; self.n_xmax = 1
        self.n_ymin = 10**-6.5;  self.n_ymax = 10**0.5
        self.n_xlabel = r'$\rm{scale\ factor}$'
        self.n_ylabel = r'$M/M(a=1)$'
        self.xlog = False; self.ylog = True

    def _analyze(self,hpath):
        mb = self.mbplug.read(hpath)
        scale= mb['scale']
        mvir = mb['mvir']
        msam = mb['sam_mvir'] #not used, not exactly sure what it is yet but could be useful
        
        output = {}
        output['a_lmm'] = mb[-1]['scale_of_last_MM']

        a_half = self.calc_a_half(scale,mvir)
        output['a_half'] = a_half
        
        rscat = haloutils.load_rscat(hpath,haloutils.get_lastsnap(hpath))
        zoomid= haloutils.load_zoomid(hpath)
        conc = rscat.ix[zoomid]['rvir']/rscat.ix[zoomid]['rs']
        output['conc'] = conc

        alpha = self.fit_exp(scale,mvir)
        output['a_exp'] = 1./(1.+np.log(2)/alpha)
        output['a_exp_alpha'] = alpha

        beta,gamma = self.fit_exppow(scale,mvir)
        output['exppow_gamma'] = gamma
        output['exppow_beta'] = beta
        output['a_gmb'] = 1./(1.+np.log(2)/(gamma-beta))

        z_betagam = self.solve_for_z_betagam(beta,gamma,z0=(1./a_half - 1.0))
        a_betagam = 1./(1.+z_betagam)
        output['a_betagam'] = a_betagam

        a_mmm,a_mmm_mass = self.calc_max_mass_merger(scale,mvir)
        output['a_mmm'] = a_mmm
        output['a_mmm_mass'] = a_mmm_mass
        
        with open(self.get_outfname(hpath),'w') as f:
            writer = csv.writer(f)
            for k,v in output.items():
                writer.writerow([k,v])

    def calc_a_half(self,scale,mvir):
        """ Smallest scale factor where main branch reaches half of a=1 mass """
        half_m = mvir[-1]/2.0
        return np.min(scale[mvir>half_m])
    def calc_max_mass_merger(self,scale,mvir):
        dM = np.concatenate([[mvir[0]],np.diff(mvir)])
        imax = np.argmax(dM)
        return scale[imax],dM[imax]

    def exp_fn(self,z,M0,alpha):
        return M0*np.exp(-alpha*z)
    def exppow_fn(self,z,M0,beta,gamma):
        return M0*(1.+z)**beta * np.exp(-gamma*z)
    def fit_exp(self,scale,mvir,p0=(1.0)):
        z = 1./scale - 1.0
        logM = np.log10(mvir)
        logM0 = np.log10(mvir[-1])
        loge = np.log10(np.exp(1.0))
        def logexp(z,alpha):
            return logM0 - loge*alpha*z
        popt,pcov = optimize.curve_fit(logexp,z,logM,p0=p0)
        return popt[0]
    def fit_exppow(self,scale,mvir,p0=(-2.0,0.50)):
        z = 1./scale - 1.0
        logM = np.log10(mvir)
        logM0 = np.log10(mvir[-1])
        loge = np.log10(np.exp(1.0))
        def logexppow(z,beta,gamma):
            return logM0 + beta*np.log10(1.+z) - loge*gamma*z
        popt,pcov = optimize.curve_fit(logexppow,z,logM,p0=p0)
        return popt[0],popt[1]
    def solve_for_z_betagam(self,beta,gamma,z0 = 1.0):
        def betagam(z):
            return (1.+z)**beta * np.exp(-gamma*z) - 0.5
        z_betagam = optimize.fsolve(betagam,z0)
        return z_betagam[0]

    def _read(self,hpath):
        with open(self.get_outfname(hpath),'r') as f:
            reader = csv.reader(f)
            output = dict([k,float(v)] for k,v in reader)
        return output
    def _plot(self,hpath,data,ax,lx=None,labelon=False,normtohost=False,
              plottime=False,**kwargs):
        #plot mass accr M(a)
        mb = self.mbplug.read(hpath)
        scale = mb['scale']; mvir = mb['mvir']
        if plottime:
            xval = bconversions.GetTime(scale)
        else:
            xval = scale

        if lx != None: color = self.colordict[lx]
        else:
            if 'color' in kwargs: color = kwargs.pop('color')
            else: color = None
        ax.plot(xval,mvir,color=color,**kwargs)

        #label all the different formation times
        raise NotImplementedError

    def plot_exp(self,hpath,ax,**kwargs):
        pass
    def plot_exppow(self,hpath,ax,**kwargs):
        pass

def formation_tabfn(hpath):
    if hpath==None: return None
    fplug = FormationTimePlugin()
    qtys = fplug.qtynames
    ddict = fplug.read(hpath)
    data = [ddict[qty] for qty in qtys]
    names = qtys
    formats = [np.float for qty in qtys]
    return data,names,formats

if __name__=="__main__":
    tab = haloutils.tabulate(formation_tabfn,exclude_hids=[388476])
    gtab = haloutils.tabulate(formation_tabfn,exclude_hids=[388476,94687,1599988])
    #tab = haloutils.tabulate(formation_tabfn,lx=13)

    import seaborn as sns
    g = sns.PairGrid(tab,diag_sharey=False)
    g.map_lower(plt.scatter)
    g.map_upper(plt.scatter)
    g.map_diag(plt.hist)
    g.savefig('7-13/formation_time.png')

    g = sns.PairGrid(gtab,diag_sharey=False)
    g.map_lower(plt.scatter)
    g.map_upper(plt.scatter)
    g.map_diag(plt.hist)
    g.savefig('7-13/formation_time2.png')

def run_all(recalc=False):
    plug = FormationTimePlugin()
    hids = haloutils.cid2hid.values()
    for hid in hids:
        for lx in [11,12,13,14]:
            hpath = haloutils.get_hpath_lx(hid,lx)
            if hpath == None: continue
            d = plug.read(hpath,recalc=recalc)
