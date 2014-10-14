import numpy as np
import pylab as plt
import haloutils
import asciitable

if __name__=="__main__":
    vmaxbins = np.logspace(-1,3,101)
    varr = vmaxbins[1:]
    
    snap = 83
    hid = 1327707
    lxlist = [14,13,12,11]

    colordict = {14:'m',13:'g',12:'r',11:'b'}
    fig,axarr = plt.subplots(3,2,sharex=True,sharey=True,figsize=(8,11))
    axlist = axarr.reshape(-1)

    snaplist = [82,126,144,169,204,255]

    for lx in lxlist:
        print "LX =",lx
        hpath = haloutils.get_hpath(hid,'bb',lx,4)
        zoomid = haloutils.load_zoomid(hpath)
        mtc = haloutils.load_mtc(hpath,haloids=[zoomid])
        mb = mtc[0].getMainBranch()

        for snap,ax in zip(snaplist,axlist):
            rscat = haloutils.load_rscat(hpath,snap)
            snapid = mb[mb['snap']==snap]['origid'][0]
            print "z=%3.2f" % (1./rscat.scale-1),zoomid,snapid

            subs = rscat.get_all_subhalos_within_halo(snapid)
            svmax = np.array(subs['vmax'])
            #srmax = np.array(subs['rvmax'])
            h,x = np.histogram(svmax,bins=vmaxbins)
            Nvmax = np.cumsum(h[::-1])[::-1]
            
            ii = varr > np.min(svmax)
            ax.plot(varr[ii],Nvmax[ii],color=colordict[lx])
            ax.text(30,10**4,'z=%3.2f' % (1./rscat.scale - 1))

    for ax in axlist:
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlim((.3,100))
        ax.set_ylim((1,10**5))
    fig.savefig("zNvmax.png",bbox_inches='tight')
    plt.show()
