# -*- coding: utf-8 -*-
import wrf_dim_info #local function
from time import time as cpu_time
import numpy as np
import os
#from os.path import isfile, join
import glob
import Ngl; import Nio 
import netCDF4 as nc
import xarray as xr
from wrf import getvar,ALL_TIMES,interplevel
from numba import jit, vectorize, int32, float32, types
from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor
#from joblib import Parallel, delayed

#print(Nio.__version__)
##._FillValue = 9.96920996839e+36
def interp_var(lev):
    global q_all, z
    result = interplevel(q_all, z, lev)
    return result

def delta_z():
    kk = 900.0 #m bottom
    nk = 122   #top 13.0 km
    lev = np.zeros(nk,dtype=np.float32)
    for k in range(nk):
        lev[k] = kk; #print(lev[k])
        kk += 100.0
    return lev, nk

def get_variable(ff,time):
    z = getvar(ff, "z", timeidx=time)
    qc = getvar(ff, "QCLOUD", timeidx=time) #kg kg-1
    qi = getvar(ff, "QICE", timeidx=time)
    qs = getvar(ff, "QSNOW", timeidx=time)
    temp = qc + qi + qs 
    return temp,z 

def cloud(ff):
    global lev, nk
    (t,k,y,x,tt) = wrf_dim_info.info(ff);del t,k,tt
    """ interpolation """
    tot = np.zeros(shape=(nk,y,x), dtype=np.float32)
    with ProcessPoolExecutor(max_workers=12) as executor:
         for k,output in enumerate(executor.map(interp_var, lev)):
             tot[k,:,:] = output; #print(k)
 
    tot = np.where(np.isnan(tot) , -999.9, tot) 
    """ cloud fraction """
    cldfra = np.zeros(shape=(nk,y,x), dtype=int) 
    cldfra = np.where(tot > 10.0**-5, 1, 0); del tot
    #print(cldfra.shape) 
    return cldfra
#--------------------------------------------------------------------
start_time = cpu_time()
files = glob.glob("wrfout_d03*") #list
#print(type(files))
nf = len(files); #print(nf); print(files[0:nf])
#files = np.array(files) #numpy.ndarray
#print(files.shape)
(lev,nk) = delta_z()#dz for interpolation

#a = nc.MFDataset(files)
for nid,f in enumerate(files[0:nf-1]):
    File_out = str(f)[11:24]+"_cloud_char.bin"; print("output:",File_out)
    newFile = open(File_out,"wb")
    a = nc.Dataset(f)
    (ntimes, nlev, nlat, nlon, times) = wrf_dim_info.info(a)
    print("File:",nid,f,"\nDomain:{}".format(ntimes)+"|{}".format(nlev)+\
    "|{}".format(nlat)+"|{}".format(nlon)," .......Get Dimensions!!")
    for t in range(ntimes):
        print("\tCloud variables Cal. Time:"+str(times[t].values)[0:19])
        (q_all,z) = get_variable(a,t) #global var.
        cldfra = cloud(a); #print(cldfra.shape)
        c_thick = np.zeros(shape=(nlat,nlon),dtype=np.float32) 
        c_top = np.zeros(shape=(nlat,nlon),dtype=np.float32) 
        c_base = np.zeros(shape=(nlat,nlon),dtype=np.float32) 
        for j in range(nlat):
            for i in range(nlon):
                c_thick[j,i] = max(float( sum(cldfra[:,j,i])-1 ) *100.0, 0.0)
                for k in range(1,nk):
                    if (cldfra[k,j,i] == 1 and cldfra[k-1,j,i] == 0 \
                        and c_base[j,i] == 0.0):
                        c_base[j,i] = lev[k] #get cloud base
                   #------------------------------------------------- 
                    if (k < nk-1):
                       if (cldfra[k+1,j,i] == 0 and cldfra[k,j,i] == 1):
                           c_top[j,i] = lev[k]
                    elif (k == nk-1 and cldfra[k,j,i] == 1):
                          c_top[j,i] = lev[k] 
        #------------------------------------------------------------ 
        #print(str(cldfra[:,190,1]))
        del cldfra 
        #c_thick = c_top - c_base;# print(c_thick.shape)
        print(str(c_thick[190,1]),str(c_top[190,1]),str(c_base[190,1])) 
        newFile.write(bytearray(c_thick))
        newFile.write(bytearray(c_top))
        newFile.write(bytearray(c_base))
    del a #loop t


end_time = cpu_time(); end_time = (end_time - start_time)/60.0
print("cal_cloud.py has done!\nTime elapsed: {:.2f}".format(end_time), "mins.")
Ngl.end()
#Nio.end()
#exit()
