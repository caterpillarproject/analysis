import haloutils
from caterpillaranalysis import *
from caterpillarplot import *
import pylab as plt
from optparse import OptionParser

if __name__=="__main__":
    parser = OptionParser()
    (options,args) = parser.parse_args()
    firsttwelve = get_haloidlist(1)
    if len(args) == 0:
        #plug = NvmaxPlugin()
        #plug = MassAccrPlugin()
        #stackplot(firsttwelve,14,plug,lw=2,autocolor=1,figfilename='testcolor1.png',normtohost=True)
        #stackplot(firsttwelve,14,plug,lw=2,autocolor=2,figfilename='testcolor2.png',normtohost=True)
        #stackplot(firsttwelve,14,plug,lw=2,autocolor=3,figfilename='testcolor3.png',normtohost=True)
        arg = None
    else:
        arg = int(args[0])

    myrecalc = False
    if arg == 1:
        Nvmax      = NvmaxPlugin()
        convergeplot(1,Nvmax,lw=2,figfilename='NvmaxLX_s1.png',recalc=myrecalc)
        convergeplot(1,Nvmax,lw=2,figfilename='normNvmaxLX_s1.png',recalc=False,normtohost=True)
        stackplot(firsttwelve,14,Nvmax,figfilename='stackLX14_normNvmax_s1.png',lw=2,color='k',alpha=.2,recalc=False,normtohost=True)
    if arg == 2:
        SHMF       = SHMFPlugin()
        intSHMF    = IntegrableSHMFPlugin()
        BoundSHMF  = BoundSHMFPlugin()
        convergeplot(1,SHMF,lw=2,figfilename='SHMFLX_s1.png',recalc=myrecalc)
        convergeplot(1,SHMF,lw=2,figfilename='normSHMFLX_s1.png',recalc=False,normtohost=True)
        convergeplot(1,BoundSHMF,lw=2,figfilename='BoundSHMFLX_s1.png',recalc=False)
        convergeplot(1,BoundSHMF,lw=2,figfilename='normBoundSHMFLX_s1.png',recalc=False,normtohost=True)
        convergeplot(1,intSHMF,lw=2,figfilename='intSHMFLX_s1.png',recalc=False)
        stackplot(firsttwelve,14,SHMF,figfilename='stackLX14_normSHMF_s1.png',lw=2,color='k',alpha=.2,recalc=False,normtohost=True)
        stackplot(firsttwelve,14,BoundSHMF,figfilename='stackLX14_normBoundSHMF_s1.png',lw=2,color='k',alpha=.2,recalc=False,normtohost=True)
    if arg == 3:
        Profile    = ProfilePlugin()
        VelProfile = VelocityProfilePlugin()
        convergeplot(1,Profile,figfilename='rhor2LX_s1.png',recalc=myrecalc)
        convergeplot(1,VelProfile,figfilename='velprofLX_s1.png',recalc=False)
        convergeplot(1,VelProfile,figfilename='normvelprofLX_s1.png',recalc=False,normtohost=True)
        stackplot(firsttwelve,14,VelProfile,figfilename='stackLX14_normvelprof_s1.png',color='k',alpha=.2,recalc=False,normtohost=True)
    if arg == 4:
        #MassAccr   = MassAccrPlugin()
        LinMassAccr= LinearMassAccrPlugin()
        #convergeplot(1,MassAccr,figfilename='massaccrLX_s1.png',recalc=myrecalc)
        #convergeplot(1,MassAccr,figfilename='normmassaccrLX_s1.png',recalc=False,normtohost=True)
        #stackplot(firsttwelve,14,MassAccr,figfilename='stackLX14_massaccr_s1.png',lw=2,color='k',alpha=.2,recalc=False)
        #stackplot(firsttwelve,14,MassAccr,figfilename='stackLX14_normmassaccr_s1.png',lw=2,color='k',alpha=.2,recalc=False,normtohost=True)
        convergeplot(1,LinMassAccr,figfilename='linmassaccrLX_s1.png',recalc=False)
        convergeplot(1,LinMassAccr,figfilename='normlinmassaccrLX_s1.png',recalc=False,normtohost=True)
        stackplot(firsttwelve,14,LinMassAccr,figfilename='stackLX14_linmassaccr_s1.png',lw=2,color='k',alpha=.2,recalc=False)
        stackplot(firsttwelve,14,LinMassAccr,figfilename='stackLX14_normlinmassaccr_s1.png',lw=2,color='k',alpha=.2,recalc=False,normtohost=True)
    if arg == 5:
        SubPhase   = SubPhaseContourPlugin()
        convergeplot(1,SubPhase,whichlx=[13,14],figfilename='subphaseLX_s1.png',recalc=myrecalc)
    if arg == 6:
        SubProfile     = SubProfilePlugin()
        SubVelProfile  = SubVelocityProfilePlugin()
        convergeplot(1,SubProfile,whichlx=[12],figfilename='subprofLX12_s1.png',recalc=myrecalc)
        convergeplot(1,SubProfile,whichlx=[13],figfilename='subprofLX13_s1.png',recalc=myrecalc)
        convergeplot(1,SubProfile,whichlx=[14],figfilename='subprofLX14_s1.png',recalc=myrecalc)
        convergeplot(1,SubVelProfile,whichlx=[12],figfilename='subvelprofLX12_s1.png',recalc=False)
        convergeplot(1,SubVelProfile,whichlx=[13],figfilename='subvelprofLX13_s1.png',recalc=False)
        convergeplot(1,SubVelProfile,whichlx=[14],figfilename='subvelprofLX14_s1.png',recalc=False)
        convergeplot(1,SubVelProfile,whichlx=[12],figfilename='normsubvelprofLX12_s1.png',recalc=False,normtohost=True)
        convergeplot(1,SubVelProfile,whichlx=[13],figfilename='normsubvelprofLX13_s1.png',recalc=False,normtohost=True)
        convergeplot(1,SubVelProfile,whichlx=[14],figfilename='normsubvelprofLX14_s1.png',recalc=False,normtohost=True)
    if arg == 7:
        SubRad         = SubhaloRadialPlugin()
        intSubRad      = IntegrableSubhaloRadialPlugin()
        SubRadByMass   = SubhaloRadialByMassPlugin()
        intSubRadByMass= IntegrableSubhaloRadialByMassPlugin()
        convergeplot(1,SubRad,figfilename='subradLX_s1.png',recalc=myrecalc)
        convergeplot(1,intSubRad,whichlx=[14],figfilename='intsubradLX14_s1.png',recalc=False)
        convergeplot(1,SubRadByMass,whichlx=[14],figfilename='subradbymassLX14_s1.png',recalc=False)
        convergeplot(1,intSubRadByMass,whichlx=[14],figfilename='intsubradbymassLX14_s1.png',recalc=False)
    if arg == 8:
        SubRadMass     = SubhaloRadialMassPlugin()
        SubRadMassFrac = SubhaloRadialMassFracPlugin()
        convergeplot(1,SubRadMass,figfilename='subradmassLX_s1.png',recalc=myrecalc)
        convergeplot(1,SubRadMassFrac,figfilename='subradmassfracLX_s1.png',recalc=False)
        stackplot(firsttwelve,14,SubRadMassFrac,figfilename='stackLX14_subradmassfrac_s1.png',color='k',alpha=.2,recalc=False)
    if arg == 9:
        SubRadBoundMass = SubhaloRadialBoundMassPlugin()
        SubRadBoundMassFrac = SubhaloRadialBoundMassFracPlugin()
        convergeplot(1,SubRadBoundMass,figfilename='subradboundmassLX_s1.png',recalc=myrecalc)
        convergeplot(1,SubRadBoundMassFrac,figfilename='subradboundmassfracLX_s1.png',recalc=False)
        stackplot(firsttwelve,14,SubRadBoundMassFrac,figfilename='stackLX14_subradboundmassfrac_s1.png',color='k',alpha=.2,recalc=False)

        SubRadBound95Mass = SubhaloRadialBound95MassPlugin()
        SubRadBound95MassFrac = SubhaloRadialBound95MassFracPlugin()
        convergeplot(1,SubRadBound95Mass,figfilename='subradbound95massLX_s1.png',recalc=myrecalc)
        convergeplot(1,SubRadBound95MassFrac,figfilename='subradbound95massfracLX_s1.png',recalc=False)
        stackplot(firsttwelve,14,SubRadBound95MassFrac,figfilename='stackLX14_subradbound95massfrac_s1.png',color='k',alpha=.2,recalc=False)

        SubRadBound99Mass = SubhaloRadialBound99MassPlugin()
        SubRadBound99MassFrac = SubhaloRadialBound99MassFracPlugin()
        convergeplot(1,SubRadBound99Mass,figfilename='subradbound99massLX_s1.png',recalc=myrecalc)
        convergeplot(1,SubRadBound99MassFrac,figfilename='subradbound99massfracLX_s1.png',recalc=False)
        stackplot(firsttwelve,14,SubRadBound99MassFrac,figfilename='stackLX14_subradbound99massfrac_s1.png',color='k',alpha=.2,recalc=False)

    #testhid = firsttwelve[2]
    #haloplot(testhid,12,pluglist=[Nvmax,SHMF,Profile,MassAccr],recalc=myrecalc)
    #haloplot(testhid,14,pluglist=[SubPhase],recalc=myrecalc)
    #haloplot(testhid,12,pluglist=[SubProfile],recalc=myrecalc)
    #haloplot(1130025,14,pluglist=[SubVelProfile])

    plt.close('all')
