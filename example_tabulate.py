import numpy as np
import os,subprocess,sys,time
import asciitable
import haloutils
import functools

import pandas as pd

def example_tabfn(hpath):
    if hpath==None: return None
    try:
        hid = haloutils.get_parent_hid(hpath)
        idstr = haloutils.hidstr(hid)
    except:
        print "  test_tabulate2 failed on "+hpath
        return None
    data = (hid,idstr)
    names = ['id','idstr']
    formats = [int,object]
    return data,names,formats

if __name__=="__main__":
    df = haloutils.tabulate(example_tabfn,numprocs=1)
