import haloutils
from caterpillaranalysis import *
from caterpillarplot import *

if __name__=="__main__":
    Nvmax = NvmaxPlugin()
    SHMF = SHMFPlugin()
    Profile = ProfilePlugin()
    MassAccr = MassAccrPlugin()
    #convergeplot(1,Nvmax,lw=2,figfilename='NvmaxLX_s1.png')
    #convergeplot(1,SHMF,lw=2,figfilename='SHMFLX_s1.png')
    #convergeplot(1,Profile,figfilename='rhor2LX_s1.png')
    #convergeplot(1,MassAccr,figfilename='massaccrLX_s1.png')

    firsttwelve = get_haloidlist(1)
    stackplot(firsttwelve,12,MassAccr,color='k',alpha=.2)
    
    testhid = firsttwelve[0]
    #haloplot(testhid,12,pluglist=[Nvmax,SHMF,Profile,MassAccr])
