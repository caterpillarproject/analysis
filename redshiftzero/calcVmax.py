import numpy as np
import haloutils
import asciitable
from optparse import OptionParser

def calcNvmax(hpath):
    minvmax = 10**-1; maxvmax = 10**3
    vmaxbins = np.logspace(-1,3,101)
    filename = 'Nvmax.dat'

    numsnaps = haloutils.get_numsnaps(hpath)
    rscat = haloutils.load_rscat(hpath,numsnaps-1)
    zoomid = haloutils.load_zoomid(hpath)

    #subs = rscat.get_all_subhalos_from_halo(zoomid)
    soft = 1000*haloutils.load_soft(hpath)
    eps = soft/2.8
    hostpos = np.array(rscat.ix[zoomid][['posX','posY','posZ']])
    halopos = np.array(rscat[['posX','posY','posZ']])
    halodr  = np.sqrt(np.sum((halopos-hostpos)**2,1))
    ii = (halodr < 0.4)
    subs = rscat.data[ii]
    svmax = np.array(subs['vmax'])
    srmax = np.array(subs['rvmax'])
    svmaxp = svmax * np.sqrt(1+(eps/srmax)**2)

    h,x = np.histogram(svmax,bins=vmaxbins)
    Nvmax = np.cumsum(h[::-1])[::-1]
    h,x = np.histogram(svmaxp,bins=vmaxbins)
    Nvmaxp = np.cumsum(h[::-1])[::-1]

    try:
        scat = haloutils.load_scat(hpath)
        #TODO
        bestgroup = 0

        ssvmax = scat.sub_vmax[0:scat.group_nsubs[0]]
        ssrmax = scat.sub_vmaxrad[0:scat.group_nsubs[0]]
        ssvmaxp = ssvmax*np.sqrt(1+((eps/1000.)/ssrmax)**2)

        h,x = np.histogram(ssvmax,bins=vmaxbins)
        sNvmax = np.cumsum(h[::-1])[::-1]
        h,x = np.histogram(ssvmaxp,bins=vmaxbins)
        sNvmaxp = np.cumsum(h[::-1])[::-1]
    except: #LX14 doesn't have subfind yet
        ssvmax = 0
        sNvmax = np.zeros(len(Nvmax))
        sNvmaxp = np.zeros(len(Nvmax))
    
    f = open(hpath+'/'+filename,'w')
    f.write(str(np.min(svmax))+" "+str(np.min(ssvmax))+'\n')
    for v,N,sN,Np,sNp in zip(vmaxbins[1:],Nvmax,sNvmax,Nvmaxp,sNvmaxp):
        f.write(str(v)+" "+str(N)+" "+str(sN)+" "+str(Np)+" "+str(sNp)+'\n')
    f.close()

if __name__=="__main__":
    parser = OptionParser()
    options,args = parser.parse_args()
    hpath = args[0]
    calcNvmax(hpath)
    
