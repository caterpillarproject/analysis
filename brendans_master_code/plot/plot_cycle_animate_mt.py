import os
import subprocess as sub
from brendanlib.grifflib import getcaterpillarcandidates


"""
Contact: Brendan Griffen <brendan.f.griffen@gmail.com>

This script will submit N background processes set by the user.

It will then call plot_parent_mt_animation.py.

"""

# SET MAXIMUM NUMBER OF RUNS
max_nruns = 6

nrun = 0

# SET WHERE CANDIDATE LIST IS TAKEN FROM
candidatelist = 


for i,haloid in enumerate(candidatelist):
    if not os.path.isdir("/bigbang/data/bgriffen/caterpillar/mergertree/animation/H" + str(haloid)):
        if nrun < max_nruns:
        print 
            print "Running Halo:",haloid
            cmd = "nohup python plot_parent_mt_animation.py " + str(haloid) + " > H" + str(haloid) + ".out &"
            print cmd
            sub.call([cmd],shell=True)

        nrun += 1
        