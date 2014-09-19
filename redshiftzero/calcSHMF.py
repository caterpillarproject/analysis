import numpy as np
import haloutils
import asciitable
from optparse import OptionParser
import MassFunctions

def calcSHMF(hpath):
    logMmax = 10.6; logMmin = 4.0
    histrange = np.arange(logMmin,logMmax,0.2)
    filename = 'SHMF.dat'

    numsnaps = haloutils.get_numsnaps(hpath)
    rscat = haloutils.load_rscat(hpath,numsnaps-1)
    zoomid = haloutils.load_zoomid(hpath)

    subs = rscat.get_all_subhalos_from_halo(zoomid)
    subM = np.array(subs['mvir'])
    x,y = MassFunctions.MassFunc_dNdM(subM,histrange)
    
    try:
        scat = haloutils.load_scat(hpath)
        #TODO
        bestgroup = 0

        ssubM = scat.sub_mass[0:scat.group_nsubs[0]]*10**10
        sx,sy = MassFunctions.MassFunc_dNdM(ssubM,histrange)
    except: #LX14 doesn't have subfind yet
        sx = np.zeros(len(x))
        sy = np.zeros(len(y))
    
    f = open(hpath+'/'+filename,'w')
    for a,b,sa,sb in zip(x,y,sx,sy):
        f.write(str(a)+' '+str(b)+' '+str(sa)+' '+str(sb)+'\n')
    f.close()

if __name__=="__main__":
    parser = OptionParser()
    options,args = parser.parse_args()
    hpath = args[0]
    calcSHMF(hpath)
    
