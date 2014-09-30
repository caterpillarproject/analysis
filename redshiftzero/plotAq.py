import readhalos.RSDataReader as RDR
import readhalos.readsubf as rsf
import pylab as plt
import numpy as np

def plot_AqA():
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
        ax2.set_xlabel('subf vmax (corr)'); ax2.set_ylabel('N(>vmax)')
        #print 'level',level,'mass',cat.ix[zoomid]['mvir'],'totnp',cat.ix[zoomid]['total_npart'],'minvmax',np.min(svmax)
    for ax in [axarr[0,0],axarr[0,1],axarr[0,2],axarr[1,0],axarr[1,1],axarr[1,2]]:
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlim((.1,100)); ax.set_ylim((1,10**5))
    plt.show()


if __name__=="__main__":
    #plotAqA()

    vmaxbins = np.logspace(-1,3,101)
    varr = vmaxbins[1:]
    namelist = ['A','B','C','D','E','F']
    allreslist = [[2,3,4,5],
                  [2,4,5],
                  [2,4],
                  [2,4,5],
                  [2,4],
                  [2,3]]
    allsnaplist = [[1023,511,1023,127],
                   [127,127,200],
                   [127,127],
                   [127,127,200],
                   [127,127],
                   [111,110]]

    colordict = {2:'m',3:'g',4:'r',5:'b'}
    epsdict = {2:0.0658,3:0.1205,4:0.3425,5:0.6849}

    fig1,axarr1 = plt.subplots(3,2,figsize=(8,12),sharex=True,sharey=True)
    fig2,axarr2 = plt.subplots(3,2,figsize=(8,12),sharex=True,sharey=True)
    for i in range(6):
        name = namelist[i]
        irow,icol = divmod(i,2)
        ax1 = axarr1[irow,icol]
        ax2 = axarr2[irow,icol]
        reslist = allreslist[i]
        snaplist= allsnaplist[i]
        for level,snap in zip(reslist,snaplist):
            print name,level,snap
            col = colordict[level]
            eps = epsdict[level]

            hpath = "/bigbang/data/AnnaGroup/aquarius/Aq-"+name+"/"+str(level)
            try:
                scat = rsf.subfind_catalog(hpath,snap)
            except:
                print "ERROR"
                continue

            svmax = scat.sub_vmax[0:scat.group_nsubs[0]]
            srmax = scat.sub_vmaxrad[0:scat.group_nsubs[0]]
            svmaxp= svmax * np.sqrt(1+(eps/1000./srmax)**2)
            h,x = np.histogram(svmax,bins=vmaxbins)
            sNvmax = np.cumsum(h[::-1])[::-1]
            h,x = np.histogram(svmaxp,bins=vmaxbins)
            sNvmaxp = np.cumsum(h[::-1])[::-1]

            ii = varr >= np.min(svmax)
            ax1.plot(varr[ii],sNvmax[ii],color=col)
            ax2.plot(varr[ii],sNvmaxp[ii],color=col)
            
            if irow==2:
                ax1.set_xlabel('subf vmax'); ax2.set_xlabel('subf vmax (corr)')
            if icol==0:
                ax1.set_ylabel('N(>vmax)'); ax2.set_ylabel('N(>vmax)')
        ax1.text(40,10**4,'Aq'+name,fontsize=14)
        ax2.text(40,10**4,'Aq'+name,fontsize=14)

        for ax in [ax1,ax2]:
            ax.set_xscale('log'); ax.set_yscale('log')
            ax.set_xlim((.3,100)); ax.set_ylim((1,10**5))

    fig1.savefig("Aq-all_Nvmax.png",bbox_inches='tight')
    fig2.savefig("Aq-all_Nvmaxp.png",bbox_inches='tight')
    plt.show()
