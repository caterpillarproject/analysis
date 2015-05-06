import numpy as np
import pylab as plt
import pandas as pd
#import seaborn as sns
#sns.set_style('ticks')
#sns.set_context('poster')

import haloutils

def tab_conc_mass(hpath):
    if hpath==None: return None
    if not haloutils.check_last_rockstar_exists(hpath): return None
    lastsnap = haloutils.get_numsnaps(hpath)-1
    rscat = haloutils.load_rscat(hpath,lastsnap)
    zoomid= haloutils.load_zoomid(hpath)
    host = rscat.ix[zoomid]
    mvir = host['mvir']
    rvir = host['rvir']
    mass = host['mgrav']
    m200c = host['altm2']
    conc = host['rvir']/host['rs']
    conc_klypin = host['rvir']/host['rs_klypin']
    Xoff = host['Xoff']; Voff = host['Voff']; TU = host['T/|U|']
    data = conc,mass,conc_klypin,mvir,m200c,Xoff,Voff,TU,rvir
    names = ['conc','mass','conc_klypin','mvir_all','m200c','Xoff','Voff','T/U','rvir']
    formats = [np.float for i in range(len(data))]
    return data,names,formats

def dutton_maccio_2014(logM,use200=False):
    """
    @param logM: log10(M) in units of h^-1 Msun 
    @param use200: if True, fit for M200c instead of Mvir
    """
    if use200:
        return 0.905 - 0.101 * (logM - 12)
    return 1.025 - 0.097*(logM-12)
    

if __name__=="__main__":
    lx = 14
    #df = haloutils.tabulate(tab_conc_mass,lx=lx,numprocs=4,exclude_hids=[94687],savefile='concmass'+str(lx)+'.csv')
    df = pd.read_csv('concmass'+str(lx)+'.csv',index_col=0)
    fig,ax = plt.subplots(figsize=(8,8))
    sc = ax.scatter(np.log10(df['mass']),np.log10(df['conc']),c=df['Xoff']/df['rvir'],s=36)
    logM = np.linspace(11.7,12.2)
    ax.plot(logM,dutton_maccio_2014(logM),'k')
    ax.plot(logM,dutton_maccio_2014(logM)+.11,'k:')
    ax.plot(logM,dutton_maccio_2014(logM)-.11,'k:')
    cb = plt.colorbar(sc)
    cb.set_label('xoff')
    ax.set_xlim((11.7,12.2))
    ax.set_xlabel(r'$\log M_{\rm vir} [M_\odot/h]$')
    ax.set_ylabel(r'$\log R_{\rm vir}/r_s$')
    plt.show()
