import numpy as np
import pylab as plt
import os,asciitable
import haloutils,sheetplot
import analysisplugin as plug

class ReaderBase(object):
    def get_filename(self,hpath):
        fname = hpath+'/'+plug.OUTPUTFOLDERNAME+'/'+self.filename
        if not os.path.exists(fname): raise IOError
        return fname

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

class ProfileReader(ReaderBase):
    def __init__(self):
        self.filename='rsprofile.dat'
    def __call__(self,hid):
        lxhpaths = haloutils.get_lxlist(hid,gethpaths=True)
        lxlist = haloutils.get_lxlist(hid)

        rlist = []; rholist = []; p03rlist = []
        rvirlist = []; r200clist = []
        for hpath in lxhpaths:
            thisfilename = self.get_filename(hpath)
            data = np.array(asciitable.read(thisfilename,delimiter=" ",data_start=1))
            rlist.append(1000.*data['col1'])
            rholist.append(data['col2'])
            f = open(thisfilename,'r')
            p03r,rvir,r200c,halomass = f.readline().split(" ")
            p03r = 1000.*float(p03r); rvir = float(rvir); r200c = float(r200c)
            f.close()
            p03rlist.append(p03r); rvirlist.append(rvir)
            r200clist.append(r200c)
        return hid,lxlist,rlist,rholist,p03rlist,rvirlist,r200clist

class ProjReader(ReaderBase):
    def __init__(self,lx,snap=255):
        self.snap=snap
        self.lx=lx
        if snap==255:
            self.filename='xyproj.npy'
        else:
            self.filename='xyproj'+str(snap).zfill(3)+'.npy'
    def __call__(self,hid):
        hpaths = haloutils.get_available_hpaths(hid)
        hpathmatches = []
        for hpath in hpaths:
            ictype,lx,nv=haloutils.get_zoom_params(hpath)
            if lx==self.lx: hpathmatches.append(hpath)
        if len(hpathmatches) != 1: 
            raise IOError("More than one LX matches for "+haloutils.hidstr(hid))
        hpath = hpathmatches[0]
        im,width = np.load(self.get_filename(hpath))
        return hid,self.lx,im,width
