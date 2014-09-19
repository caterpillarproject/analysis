import readhalos.RSDataReader as RDR
import readhalos.readsubf as rsf
import pylab as plt
import numpy as np

if __name__=="__main__":
    vmaxbins = np.logspace(-1,3,101)
    varr = vmaxbins[1:]
    levellist = [2,3,4,5]
    snaplist = [126,511,1023,127]
    zoomidlist = [33387,157193,23676,2855] #halos
    #zoomidlist = [33387,156883,23628,2848] #halossoft
    fig,axarr = plt.subplots(2,3)
    colorlist = ['purple','green','red','blue']
    epslist = [.0658,.1205,.3425,.6849] #kpc, from Springel et al 2008
    for level,snap,zoomid,col,eps in zip(levellist,snaplist,zoomidlist,colorlist,epslist):
        if level==2: cat = RDR.RSDataReader("/bigbang/data/AnnaGroup/aquarius/Aq-A/"+str(level)+"/BACKhalos",snap,version=7)
        else: 
            cat = RDR.RSDataReader("/bigbang/data/AnnaGroup/aquarius/Aq-A/"+str(level)+"/halos",snap,version=7)
        subs = cat.get_all_subhalos_from_halo(zoomid)
        svmax = np.array(subs['vmax'])
        srmax = np.array(subs['rvmax'])
        svmaxp= svmax * np.sqrt(1+(eps/srmax)**2)
        h,x = np.histogram(svmax,bins=vmaxbins)
        Nvmax = np.cumsum(h[::-1])[::-1]
        h,x = np.histogram(svmaxp,bins=vmaxbins)
        Nvmaxp = np.cumsum(h[::-1])[::-1]
        
        if level != 2:
            ii = varr >= np.min(svmax)
            ax1 = axarr[0,0]; ax2 = axarr[1,0]
            ax1.plot(varr[ii],Nvmax[ii],color=col)
            ax1.set_xlabel('rs subs vmax (soft)'); ax1.set_ylabel('N(>vmax)')
            ax2.plot(varr[ii],Nvmaxp[ii],color=col)
            ax2.set_xlabel('rs subs vmax (soft,corr)'); ax1.set_ylabel('N(>vmax)')

        hostpos = np.array(cat.ix[zoomid][['posX','posY','posZ']])
        halopos = np.array(cat[['posX','posY','posZ']])
        halodr  = np.sqrt(np.sum((halopos-hostpos)**2,1))
        ii = (halodr < 0.4)
        subs = cat.data[ii]
        svmax = np.array(subs['vmax'])
        srmax = np.array(subs['rvmax'])
        svmaxp= svmax * np.sqrt(1+(eps/srmax)**2)
        h,x = np.histogram(svmax,bins=vmaxbins)
        Nvmax = np.cumsum(h[::-1])[::-1]
        h,x = np.histogram(svmaxp,bins=vmaxbins)
        Nvmaxp = np.cumsum(h[::-1])[::-1]
        
        if level != 2:
            ii = varr >= np.min(svmax)
            ax1 = axarr[0,1]; ax2 = axarr[1,1]
            ax1.plot(varr[ii],Nvmax[ii],color=col)
            ax1.set_xlabel('rs <400kpc vmax (soft)'); ax1.set_ylabel('N(>vmax)')
            ax2.plot(varr[ii],Nvmaxp[ii],color=col)
            ax2.set_xlabel('rs <400kpc vmax (soft,corr)'); ax1.set_ylabel('N(>vmax)')

        if level==2: subsnap = 1023
        else: subsnap = snap
        scat = rsf.subfind_catalog("/bigbang/data/AnnaGroup/aquarius/Aq-A/"+str(level),subsnap)
        svmax = scat.sub_vmax[0:scat.group_nsubs[0]]
        srmax = scat.sub_vmaxrad[0:scat.group_nsubs[0]]
        svmaxp= svmax * np.sqrt(1+(eps/1000./srmax)**2)
        h,x = np.histogram(svmax,bins=vmaxbins)
        sNvmax = np.cumsum(h[::-1])[::-1]
        h,x = np.histogram(svmaxp,bins=vmaxbins)
        sNvmaxp = np.cumsum(h[::-1])[::-1]
        
        ii = varr >= np.min(svmax)
        ax1 = axarr[0,2]; ax2 = axarr[1,2]
        ax1.plot(varr[ii],sNvmax[ii],color=col)
        ax1.set_xlabel('subf vmax'); ax1.set_ylabel('N(>vmax)')
        ax2.plot(varr[ii],sNvmaxp[ii],color=col)
        ax2.set_xlabel('subf vmax (corr)'); ax1.set_ylabel('N(>vmax)')
        #print 'level',level,'mass',cat.ix[zoomid]['mvir'],'totnp',cat.ix[zoomid]['total_npart'],'minvmax',np.min(svmax)
    for ax in [axarr[0,0],axarr[0,1],axarr[0,2],axarr[1,0],axarr[1,1],axarr[1,2]]:
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlim((.1,100)); ax.set_ylim((1,10**5))
    plt.show()
