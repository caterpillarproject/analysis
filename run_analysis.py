import haloutils
from caterpillaranalysis import *
from caterpillarplot import *
import pylab as plt

if __name__=="__main__":
    Nvmax = NvmaxPlugin()
    SHMF = SHMFPlugin()
    Profile = ProfilePlugin()
    VelProfile = VelocityProfilePlugin()
    MassAccr = MassAccrPlugin()
    SubPhase = SubPhaseContourPlugin()
    SubProfile = SubProfilePlugin()
    SubVelProfile = SubVelocityProfilePlugin()
    #convergeplot(1,Nvmax,lw=2,figfilename='NvmaxLX_s1.png')
    #convergeplot(1,SHMF,lw=2,figfilename='SHMFLX_s1.png')
    #convergeplot(1,Profile,figfilename='rhor2LX_s1.png',recalc=True)
    #convergeplot(1,VelProfile,figfilename='velprofLX_s1.png')
    #convergeplot(1,MassAccr,figfilename='massaccrLX_s1.png')
    convergeplot(1,MassAccr)
    #convergeplot(1,SubPhase,whichlx=[13,14],figfilename='subphaseLX_s1.png')
    #convergeplot(1,SubProfile,whichlx=[12],figfilename='subprofLX12_s1.png',recalc=True)
    #convergeplot(1,SubProfile,whichlx=[13],figfilename='subprofLX13_s1.png',recalc=True)
    #convergeplot(1,SubProfile,whichlx=[14],figfilename='subprofLX14_s1.png',recalc=True)
    #convergeplot(1,SubVelProfile,whichlx=[12],figfilename='subvelprofLX12_s1.png')
    #convergeplot(1,SubVelProfile,whichlx=[13],figfilename='subvelprofLX13_s1.png')
    #convergeplot(1,SubVelProfile,whichlx=[14],figfilename='subvelprofLX14_s1.png')

    firsttwelve = get_haloidlist(1)
    #stackplot(firsttwelve,14,MassAccr,color='k',alpha=.2)
    
    testhid = firsttwelve[2]
    #haloplot(testhid,12,pluglist=[Nvmax,SHMF,Profile,MassAccr])
    #haloplot(testhid,14,pluglist=[SubPhase],recalc=True)
    #haloplot(testhid,12,pluglist=[SubProfile],recalc=True)
    plt.show()