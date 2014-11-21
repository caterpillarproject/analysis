import numpy as np
import pylab as plt
import os,subprocess
import haloutils
from analysisplugin import AnalysisPluginBase
import pynbody as pnb
import mergertrees.MTCatalogue as MTC
from brendanlib.grifflib import makecolormap
import seaborn.apionly as sns
import gc,sys

#def get_open_fds():
#    '''
#    return the number of open file descriptors for current process
#    
#    .. warning: will only work on UNIX-like os-es.
#    '''
#    
#    pid = os.getpid()
#    procs = subprocess.check_output( 
#        [ "lsof", '-w', '-Ff', "-p", str( pid ) ] )
#    
#    nprocs = len( 
#        filter( 
#            lambda s: s and s[ 0 ] == 'f' and s[1: ].isdigit(),
#            procs.split( '\n' ) )
#        )
#    return nprocs

def make_movie(hid,lx):
    mycm = makecolormap()
    hpath = haloutils.get_hpath_lx(hid,lx)
    #if not os.path.exists(hpath+'/halos/trees/tree.bin'):
    #    print "No tree! Converting..."
    #    MTC.convertmt(hpath+'/trees',version=4,verbose=True)
    mtc = MTC.MTCatalogue(hpath+'/halos/trees',version=4,numHosts=1)
    mt = mtc[0]
    zoomid = haloutils.load_zoomid(hpath)
    if mt.rockstar_id != zoomid:
        raise ValueError()
    mb = mt.getMainBranch()

    outdir = 'visualizations/movies/'+haloutils.hidstr(hid)+"_LX"+str(lx)
    try:
        os.makedirs(outdir)
    except OSError as e:
        if e[0]==17: pass #File exists
        else: raise e

    for row in mb:
        snap = row['snap']
        #if snap >= 223: continue
        hpos = np.array([row['posX'],row['posY'],row['posZ']])
        snapstr = str(snap).zfill(3)
        snappath = hpath+'/outputs/snapdir_'+snapstr+'/snap_'+snapstr
        print snapstr
        
        imname = outdir+'/'+snapstr+'.png'
        sim = pnb.load(snappath)
        sim['pos'] -= hpos
        width = 0.6 #comoving Mpc/h
        im = pnb.plot.image(sim,qty='rho',units="Msol kpc^-2",
                            width=width,cmap=mycm,
                            vmin=2*10**6,vmax=7*10**8,show_cbar=False,
                            noplot=True)
                            #show_cbar=True)
        plt.xlabel(''); plt.ylabel('')
        plt.xticks([]); plt.yticks([])
        #plt.savefig(imname,bbox_inches='tight')
        plt.close('all')
        sim.close()
        sys.stdout.flush()
        np.save(outdir+'/im'+snapstr,im)
        gc.collect()

if __name__=="__main__":
    #print 11
    #make_movie(95289,11)
    print "LX12"
    make_movie(95289,12)
    #print "LX13"
    #make_movie(95289,13)
    #print "LX14"
    #make_movie(95289,14)
