# -*- coding: utf-8 -*-
import numpy as np
from wrf import getvar,ALL_TIMES

def info(ff):
    wa = getvar(ff, "wa", timeidx=0)
    nk = np.int32(len(wa[:,0,0]))
    ny = np.int32(len(wa[0,:,0]))
    nx = np.int32(len(wa[0,0,:])); del wa
    t = getvar(ff, "times", timeidx=ALL_TIMES) #-1 isn't work
    nt = np.int32(len(t))
    return nt, nk, ny, nx, t

if __name__ == '__main__':
   info()
 
