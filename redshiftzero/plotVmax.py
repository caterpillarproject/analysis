import haloutils
from NvmaxPlotter import *
from SHMFPlotter import *
from sheetplot import sheetplot

if __name__=="__main__":
    Reader =NvmaxReader()
    Plotter=NvmaxPlotter()
    Plotterp=NvmaxpPlotter()
    sPlotter=sNvmaxPlotter()
    sPlotterp=sNvmaxpPlotter()

    SHMFR = SHMFReader()
    SHMFP = SHMFPlotter()
    sSHMFP = sSHMFPlotter()
    for sheet in [1]:
        #sheetplot(sheet,SHMFR,SHMFP)
        #sheetplot(sheet,SHMFR,sSHMFP)
        sheetplot(sheet,Reader,Plotter)
        sheetplot(sheet,Reader,Plotterp)
        sheetplot(sheet,Reader,sPlotter)
        sheetplot(sheet,Reader,sPlotterp)
