import numpy as np
import matplotlib.pyplot as plt
import pynbody as pnb
from brendanlib.grifflib import makecolormap
import readhalos.readgroup as RG
import readhalos.readsubf as readsubf
from operator import itemgetter
import glob as glob
import os

halolist = [268422,121869,233776]

#halolist = [268422]
hubble = 0.6711

basepath = "/bigbang/data/AnnaGroup/caterpillar/halos/"
halolist = glob.glob(basepath + "H*")

for halo in halolist:
    halos = glob.glob(halo+"/*")
    for haloi in halos:
        numblocks = []
        path = haloi + "/outputs/"
        if os.path.isdir(path+'groups_255') and "LX14" in haloi:
            print "---------------------"
            s = readsubf.subfind_catalog(path, 255)
            contamNR = s.group_contamination_count
            contamMR = s.group_contamination_mass*10**10/hubble
            minindex = min(enumerate(s.group_contamination_count[0:3]), key=itemgetter(1))[0]
            groupmass = s.group_mass[minindex]*10**10/hubble
            groupos = s.group_pos[minindex]
            print haloi.replace(basepath,"")
            print "Group Position:",groupos
            print "Candidate FOF Mass: %0.2e" % (groupmass)
            print "Maximum Subhalo Mass %0.2e" % (max(s.sub_mass*10**10/hubble))
            print "Contatmination Num.:",contamNR[minindex]
            print "Contatmination Mass: %0.2e" % (contamMR[minindex]/0.6711)

