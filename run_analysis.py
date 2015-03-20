import haloutils
from caterpillaranalysis import *
from subradplugin import *
from caterpillarplot import *
from tbtfplugin import *
from subprofileplugin import *
import pylab as plt
from optparse import OptionParser

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("--recalc",action="store_true",dest='recalc',default=False)
    parser.add_option("--eps",action="store_true",dest='eps',default=False)
    parser.add_option("--pdf",action="store_true",dest='pdf',default=False)
    parser.add_option("-s","--sheet",action="store",type=int,dest='sheet',default=1)
    parser.add_option("--numlabel",action="store_true",dest='numlabel',default=False)
    parser.add_option("--stoponerror",action="store_true",dest='stoponerror',default=False)
    (options,args) = parser.parse_args()
    firsttwelve = get_haloidlist(1)
    ### Use this section for quick testing scripts 
    if len(args) == 0:
        #plug = SubhaloRadialSubmassFracPlugin()
        #convergeplot(1,plug,figfilename='subradsubmassfracLX_s1.png',stop_on_error=True)
        #plug = SubProfileSoftPlugin()
        #convergeplot(1,plug,whichlx=[14],figfilename='subprofcorrectionLX14_s1.png',stop_on_error=True)
        plug = TBTFSoftPlugin()
        #convergeplot(1,plug,whichlx=[14],figfilename='tbtfsoftLX14_s1.png',lw=2,stop_on_error=True)
        convergeplot(1,plug,whichlx=[13],figfilename='tbtfsoftLX13_s1.png',lw=2,stop_on_error=True)
        #plug = NvinfallPlugin()
        #convergeplot(1,plug,whichlx=[13,14],figfilename='NvinfallLX1314_s1.png',stop_on_error=True)
        arg = None
    else:
        assert len(args)==1
        arg = int(args[0])

    ### Use this section for essentially finalized analysis plots
    sheet = options.sheet; ext = '_s'+str(sheet)
    myrecalc = options.recalc
    recalckwargs = {'recalc':options.recalc,
                    'stop_on_error':options.stoponerror}
    usehaloname = options.numlabel
    if options.eps: ext+='.eps'
    elif options.pdf: ext+='.pdf'
    else: ext+='.png'
    if arg == 1:
        Nvmax      = NvmaxPlugin()
        convergeplot(sheet,Nvmax,lw=2,figfilename='NvmaxLX'+ext,usehaloname=usehaloname,**recalckwargs)
        convergeplot(sheet,Nvmax,lw=2,figfilename='normNvmaxLX'+ext,normtohost=True,usehaloname=usehaloname)
        paper_stackplot(14,Nvmax,figfilename='stackLX14_normNvmax'+ext,lw=2,normtohost=True)
    if arg == 2:
        SHMF       = SHMFPlugin()
        intSHMF    = IntegrableSHMFPlugin()
        convergeplot(sheet,SHMF,lw=2,figfilename='SHMFLX'+ext,usehaloname=usehaloname,**recalckwargs)
        convergeplot(sheet,SHMF,lw=2,figfilename='normSHMFLX'+ext,normtohost=True,usehaloname=usehaloname)
        convergeplot(sheet,intSHMF,lw=2,figfilename='intSHMFLX'+ext,usehaloname=usehaloname)
        paper_stackplot(14,SHMF,figfilename='stackLX14_normSHMF'+ext,lw=2,normtohost=True)
    if arg == 3:
        Profile    = ProfilePlugin()
        VelProfile = VelocityProfilePlugin()
        convergeplot(sheet,Profile,figfilename='rhor2LX'+ext,usehaloname=usehaloname,**recalckwargs)
        convergeplot(sheet,VelProfile,figfilename='velprofLX'+ext,usehaloname=usehaloname)
        convergeplot(sheet,VelProfile,figfilename='normvelprofLX'+ext,normtohost=True,usehaloname=usehaloname)
        paper_stackplot(14,VelProfile,figfilename='stackLX14_normvelprof'+ext,normtohost=True)
    if arg == 4:
        MassAccr   = MassAccrPlugin()
        LinMassAccr= LinearMassAccrPlugin()
        convergeplot(sheet,MassAccr,figfilename='massaccrLX'+ext,usehaloname=usehaloname,**recalckwargs)
        convergeplot(sheet,MassAccr,figfilename='normmassaccrLX'+ext,normtohost=True,usehaloname=usehaloname)
        paper_stackplot(14,MassAccr,figfilename='stackLX14_massaccr'+ext,lw=2)
        paper_stackplot(14,MassAccr,figfilename='stackLX14_normmassaccr'+ext,lw=2,normtohost=True)
        convergeplot(sheet,LinMassAccr,figfilename='linmassaccrLX'+ext,usehaloname=usehaloname)
        convergeplot(sheet,LinMassAccr,figfilename='normlinmassaccrLX'+ext,normtohost=True,usehaloname=usehaloname)
        paper_stackplot(14,LinMassAccr,figfilename='stackLX14_linmassaccr'+ext,lw=2)
        paper_stackplot(14,LinMassAccr,figfilename='stackLX14_normlinmassaccr'+ext,lw=2,normtohost=True)
    if arg == 5:
        SubPhase   = SubPhaseContourPlugin()
        convergeplot(sheet,SubPhase,whichlx=[13,14],figfilename='subphaseLX'+ext,usehaloname=usehaloname,**recalckwargs)
    if arg == 6:
        SubProfile     = SubProfilePlugin()
        SubVelProfile  = SubVelocityProfilePlugin()
        TBTF = TBTFPlugin()
        for lx in [12,13,14]:
            convergeplot(sheet,SubProfile,whichlx=[lx],figfilename='subprofLX'+str(lx)+''+ext,usehaloname=usehaloname,**recalckwargs)
            convergeplot(sheet,SubVelProfile,whichlx=[lx],figfilename='subvelprofLX'+str(lx)+''+ext,usehaloname=usehaloname)
            convergeplot(sheet,SubVelProfile,whichlx=[lx],figfilename='normsubvelprofLX'+str(lx)+''+ext,normtohost=True,usehaloname=usehaloname)
        convergeplot(sheet,TBTF,whichlx=[14],figfilename='tbtfLX'+str(14)+''+ext,usehaloname=usehaloname)

    if arg == 7:
        SubRad         = SubhaloRadialPlugin()
        intSubRad      = IntegrableSubhaloRadialPlugin()
        SubRadByMass   = SubhaloRadialByMassPlugin()
        intSubRadByMass= IntegrableSubhaloRadialByMassPlugin()
        convergeplot(sheet,SubRad,figfilename='subradLX'+ext,usehaloname=usehaloname,**recalckwargs)
        convergeplot(sheet,intSubRad,whichlx=[14],figfilename='intsubradLX14'+ext,usehaloname=usehaloname)
        convergeplot(sheet,SubRadByMass,whichlx=[14],figfilename='subradbymassLX14'+ext,usehaloname=usehaloname)
        convergeplot(sheet,intSubRadByMass,whichlx=[14],figfilename='intsubradbymassLX14'+ext,usehaloname=usehaloname)

        SubRadSubmassFrac = SubhaloRadialSubmassFracPlugin()
        convergeplot(sheet,SubRadSubmassFrac,figfilename='subradsubmassfracLX'+ext,usehaloname=usehaloname,**recalckwargs)
        paper_stackplot(14,SubRadSubmassFrac,figfilename='stackLX14_subradsubmassfrac'+ext)

    #testhid = firsttwelve[2]
    #haloplot(testhid,12,pluglist=[Nvmax,SHMF,Profile,MassAccr],**recalckwargs)
    #haloplot(testhid,14,pluglist=[SubPhase],**recalckwargs)
    #haloplot(testhid,12,pluglist=[SubProfile],**recalckwargs)
    #haloplot(1130025,14,pluglist=[SubVelProfile])

    plt.close('all')
