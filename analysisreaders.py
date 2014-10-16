import numpy as np
import pylab as plt
import os,asciitable
import haloutils,sheetplot
import analysisplugin as plug

class ReaderBase(object):
    def get_filename(self,hpath):
        return hpath+'/'+plug.OUTPUTFOLDERNAME+'/'+self.filename

class NvmaxReader(ReaderBase):
    def __init__(self):
        self.filename = 'Nvmax.dat'
    def __call__(self,hid):
        lxhpaths = haloutils.get_lxlist(hid,gethpaths=True)
        lxlist = haloutils.get_lxlist(hid)
        
        vlist = []; Nlist = []; sNlist = []
        Nplist = []; sNplist = []
        minvlist = []; sminvlist = []
        
        for hpath in lxhpaths:
            thisfilename = self.get_filename(hpath)
            data = asciitable.read(thisfilename,delimiter=' ',data_start=1)
            
            vlist.append(data['col1'])
            Nlist.append(data['col2'])
            sNlist.append(data['col3'])
            Nplist.append(data['col4'])
            sNplist.append(data['col5'])
            with open(thisfilename,'r') as f:
                split = f.readline().split(" ")
                minvlist.append(float(split[0])); sminvlist.append(float(split[1]))
        return hid,lxlist,vlist,Nlist,minvlist,sNlist,sminvlist,Nplist,sNplist

class SHMFReader(ReaderBase):
    def __init__(self):
        self.filename = 'SHMF.dat'
    def __call__(self,hid):
        lxhpaths = haloutils.get_lxlist(hid,gethpaths=True)
        lxlist = haloutils.get_lxlist(hid)
        
        xlist = []; ylist = []
        sxlist = []; sylist = []
        for hpath in lxhpaths:
            thisfilename = self.get_filename(hpath)
            data = asciitable.read(thisfilename,delimiter=' ')
            xlist.append(data['col1'])
            ylist.append(data['col2'])
            sxlist.append(data['col3'])
            sylist.append(data['col4'])
        return hid,lxlist,xlist,ylist,sxlist,sylist

