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
#from joblib import Parallel, delayed

#print(Nio.__version__)
##._FillValue = 9.96920996839e+36
def cloud(ff,time):
    (t,k,y,x,tt) = wrf_dim_info.info(ff);del t,k,tt
    qc = getvar(ff, "QCLOUD", timeidx=time) #kg kg-1
    qi = getvar(ff, "QICE", timeidx=time)
    qs = getvar(ff, "QSNOW", timeidx=time)
    temp = qc + qi + qs ;del qc, qi, qs
    """ interpolation """
    kk = 900.0 # m, bottom 0.9 km
    nk = 122 #top 13.0 km 
    lev = np.zeros(nk,dtype=np.float32) 
    z = getvar(ff, "z", timeidx=time)
    tot = np.zeros(shape=(nk,y,x), dtype=np.float32) 
    for k in range(nk):
        lev[k] = kk; #print(lev[k])
        tot[k,:,:] = interplevel(temp, z, lev[k])    
        kk += 100.0 # kk = kk + 100.0
    tot = np.where(np.isnan(tot) ,-999.9,tot)
    """ cloud fraction """
    cldfra = np.zeros(shape=(nk,y,x), dtype=int) 
    cldfra = np.where(tot > 10.0**-5, 1, 0); del tot
    #print(cldfra.shape) 
    return cldfra, lev, nk
#--------------------------------------------------------------------
start_time = cpu_time()
files = glob.glob("wrfout_d03*") #list
#print(type(files))
nf = len(files); #print(nf); print(files[0:nf])
#files = np.array(files) #numpy.ndarray
#print(files.shape)

#a = nc.MFDataset(files)
for nid,f in enumerate(files[0:1]):#nf-1]):
    File_out = str(f)[11:24]+"_cloud_char.bin"; print("output:",File_out)
    newFile = open(File_out,"wb")
    a = nc.Dataset(f)
    (ntimes, nlev, nlat, nlon, times) = wrf_dim_info.info(a)
    print("File:",nid,f,"\nDomain:{}".format(ntimes)+"|{}".format(nlev)+\
    "|{}".format(nlat)+"|{}".format(nlon)," .......Get Dimensions!!")
    for t in range(ntimes):
        print("\tCloud variables Cal. Time:"+str(times[t].values)[0:19])
        (cldfra, lev, nk) = cloud(a,t); #print(cldfra.shape)
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
        del cldfra, lev, nk 
        #c_thick = c_top - c_base;# print(c_thick.shape)
        print(str(c_thick[190,1]),str(c_top[190,1]),str(c_base[190,1])) 
        c_thick_byte = bytearray(c_thick)
        c_top_byte = bytearray(c_top)
        c_base_byte = bytearray(c_base)
        newFile.write(c_thick_byte)
        newFile.write(c_top_byte)
        newFile.write(c_base_byte)
    del a #loop t


end_time = cpu_time(); end_time = (end_time - start_time)/60.0
print("cal_cloud.py has done!\nTime elapsed: {:.2f}".format(end_time), "mins.")
Ngl.end()
#Nio.end()
#exit()
