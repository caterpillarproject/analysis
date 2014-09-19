import haloutils
from ProfilePlotter import ProfileReader,ProfilePlotter
from sheetplot import sheetplot

if __name__=="__main__":
    #whichtype = 0 #rockstarnosubs
    whichtype = 1 #rockstarall
    #whichtype = 2 #subfind: broken right now
    #whichtype = 3 #subfindradius: broken right now
    Reader = ProfileReader(whichtype)
    Plotter= ProfilePlotter(whichtype)

    for sheet in [1]:
        sheetplot(sheet,Reader,Plotter)
