import haloutils
import numpy as np
import pylab as plt
import readsnapshots.readsnapHDF5_greg as rsnap


if __name__=="__main__":
    
    vmax_scale = 10
    vmin_scale = 0.
  
    fig,axarr = plt.subplots(4,4,figsize=(15,13))
    fig.subplots_adjust(wspace=0.05)

    print "LX | vmax |  rs   | np  |    mp    |  min_vir |   20*mp    | eta [pc/h]"

    for figi,lx in enumerate([11,12,13,14]):

        nv=4;haloid=1327707;
        hpath = haloutils.get_hpath(haloid,"BB",lx,nv)
        cat = haloutils.load_rscat(hpath,255)
        header = rsnap.snapshot_header(hpath+"/outputs/snapdir_255/snap_255")
        mp = header.massarr[1]*1e10/header.hubble
        lxd = float(lx)
        soft = 1e3*100./(2**lxd)/40. #kpc
        
        cat_hosts = cat.get_hosts()

        rmax   = np.array(cat_hosts['rvmax'])
        vmax   = np.array(cat_hosts['vmax'])
        rs     = np.array(cat_hosts['rs'])
        tnpart = np.array(cat_hosts['total_npart'])
        mvir   = np.array(cat_hosts['mvir']/cat.h0)
        
        print "%2i | %3.2f | %3.2f | %3i | %3.2e | %3.2e |  %3.2e  | %3.2f" % (lx,np.min(vmax),np.min(rs),np.min(tnpart),mp,np.min(mvir),20.*mp,soft)
       
        
        x = mvir; xlabel = 'mvir'
        y = vmax; ylabel = 'vmax'
        
        fig_num = figi+1
        ax1 = axarr[0,figi]
        #ax1 = fig.add_subplot(3,4,fig_num)
        ax1.plot(x,y,'bo',markeredgecolor='b')
        ax1.set_xlabel(xlabel)

        if fig_num == 1:
            ax1.set_ylabel(ylabel)

        mass_range_high = np.logspace(0,np.log10(mp*20.),20)
        mass_range_low  = np.logspace(0,np.log10(mp*10.),20)
        mass_low_lim_mp = np.logspace(0,np.log10(mp),20)
        
        ax1.fill_between(mass_range_high, 0, vmax_scale, facecolor='green', alpha=0.25)
        ax1.fill_between(mass_range_low,  0, vmax_scale, facecolor='red',   alpha=0.25)
        ax1.fill_between(mass_low_lim_mp, 0, vmax_scale, facecolor='black', alpha=1.0)
        ax1.set_xlim([mp/10.,5e7])
        ax1.set_ylim([vmin_scale,vmax_scale])
        ax1.set_title("LX:"+str(lx))
        
        x = tnpart; xlabel = 'total # particles'
        ax2 = axarr[1,figi]
        #ax2 = fig.add_subplot(3,4,4+fig_num)
        ax2.plot(x,y,'bo',markeredgecolor='b')
        np_range = np.arange(21)
        ax2.fill_between(np_range, 0, vmax_scale, facecolor='green', alpha=0.5)
        ax2.set_xlabel(xlabel)

        if fig_num == 1:
            ax2.set_ylabel(ylabel)

        ax2.set_ylim([vmin_scale,vmax_scale])
        ax2.set_xlim([10,250])
        
        x = rs; xlabel = 'rs [kpc]'
        ax3 = axarr[2,figi]
        #ax3 = fig.add_subplot(3,4,8+fig_num)
        ax3.plot(x,y,'bo',markeredgecolor='b')
        rs_range = [0,soft]
        ax3.fill_between(rs_range, 0, vmax_scale, facecolor='green', alpha=0.25)
        ax3.set_xlabel(xlabel)
        ax3.set_xlim([0.,20])
        ax3.set_ylim([vmin_scale,vmax_scale])

        ax4 = axarr[3,figi]
        ax4.plot(rmax,vmax,'bo',markeredgecolor='b')
        ax4.fill_between(rs_range, 0, vmax_scale, facecolor='green', alpha=0.25)
        ax4.set_xlabel('rmax [kpc]')
        ax4.set_xlim([0.,20])
        ax4.set_ylim([vmin_scale,vmax_scale])

        if fig_num == 1:
            ax4.set_ylabel(ylabel)

        if fig_num > 1:
            ax1.set_yticks([])
            ax2.set_yticks([])
            ax3.set_yticks([])
            ax4.set_yticks([])
            
    fig.savefig("halofunkyness.png",bbox_inches='tight')
    #plt.show()

