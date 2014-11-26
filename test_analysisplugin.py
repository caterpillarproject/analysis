import haloutils
from analysisplugin import *
from analysisreaders import *
from analysisplotters import *
from sheetplot import sheetplot,get_haloidlist
import pylab as plt

def test_calculation():
    Nvmax = NvmaxPlugin()
    SHMF = SHMFPlugin()
    Prof = ProfilePlugin()
    Proj = ProjPlugin()
    MT = MassAccrPlugin()

    haloidlist = get_haloidlist(1)
    #haloidlist = [649861,1725139]
    #haloidlist = [1130025]
    for hid in haloidlist:
        hpaths = haloutils.get_available_hpaths(hid)
        hpaths = haloutils.restrict_halopaths(hpaths,require_rockstar=True)
        for hpath in hpaths:
            print hpath
            Nvmax(hpath)
            SHMF(hpath)
            Prof(hpath)
            #Proj(hpath)
            #MT(hpath)

    #for hpath in ["/bigbang/data/AnnaGroup/caterpillar/halos/H95289/H95289_BB_Z127_P7_LN7_LX11_O4_NV4",
    #              "/bigbang/data/AnnaGroup/caterpillar/halos/H95289/H95289_BB_Z127_P7_LN7_LX12_O4_NV4",
    #              "/bigbang/data/AnnaGroup/caterpillar/halos/H95289/H95289_BB_Z127_P7_LN7_LX13_O4_NV4",
    #              "/bigbang/data/AnnaGroup/caterpillar/halos/H95289/H95289_BB_Z127_P7_LN7_LX14_O4_NV4"]:
    #for hpath in ["/bigbang/data/AnnaGroup/caterpillar/halos/H581141/H581141_EB_Z127_P7_LN7_LX11_O4_NV4",
    #              "/bigbang/data/AnnaGroup/caterpillar/halos/H581141/H581141_EB_Z127_P7_LN7_LX12_O4_NV4",
    #              "/bigbang/data/AnnaGroup/caterpillar/halos/H581141/H581141_EB_Z127_P7_LN7_LX13_O4_NV4",
    #              "/bigbang/data/AnnaGroup/caterpillar/halos/H581141/H581141_EB_Z127_P7_LN7_LX14_O4_NV4"]:
    
def plot_analysis():
    Reader =NvmaxReader()
    Plotter=NvmaxPlotter()
    #Plotterp=NvmaxpPlotter()
    #sPlotter=sNvmaxPlotter()
    #sPlotterp=sNvmaxpPlotter()

    SHMFR = SHMFReader()
    SHMFP = SHMFPlotter()
    #sSHMFP = sSHMFPlotter()

    ProfR = ProfileReader()
    ProfP = ProfilePlotter()
    ProjR11 = ProjReader(11)
    ProjP11 = ProjPlotter(11,vmin=None,vmax=None)
    ProjR12 = ProjReader(12)
    ProjP12 = ProjPlotter(12,vmin=None,vmax=None)
    ProjR13 = ProjReader(13)
    ProjP13 = ProjPlotter(13,vmin=None,vmax=None)

    MTR = MassAccrReader()
    MTP = MassAccrPlotter()
    for sheet in [1]:
        sheetplot(sheet,SHMFR,SHMFP)
        #sheetplot(sheet,SHMFR,sSHMFP)
        sheetplot(sheet,Reader,Plotter)
        #sheetplot(sheet,Reader,Plotterp)
        #sheetplot(sheet,Reader,sPlotter)
        #sheetplot(sheet,Reader,sPlotterp)
        sheetplot(sheet,ProfR,ProfP)
        #print '11'
        #sheetplot(sheet,ProjR11,ProjP11,aspect1=True)
        #print '12'
        #sheetplot(sheet,ProjR12,ProjP12,aspect1=True)
        #print '13'
        #sheetplot(sheet,ProjR13,ProjP13,aspect1=True)
        sheetplot(sheet,MTR,MTP)
    plt.close('all')

if __name__=="__main__":
    test_calculation()
    plot_analysis()
    
