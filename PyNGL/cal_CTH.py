# -*- coding: utf-8 -*-
import sys, os
sys.path.append("/work3/artirain/Cloud_seeding/PyNGL")
#local function
import wrf_tools 
from time import time as cpu_time
import numpy as np
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
    kk = 300.#900.0 #m bottom
    nk = 128 #122   #top 13.0 km
    lev = np.zeros(nk,dtype=np.float32)
    for k in range(nk):
        lev[k] = kk; #print(lev[k])
        kk += 100.0
    return lev, nk

def get_variable(ff,time):
    z  = getvar(ff, "z", timeidx=time)
    qc = getvar(ff, "QCLOUD", timeidx=time) #[kg kg-1]
    qi = getvar(ff, "QICE", timeidx=time)
    qs = getvar(ff, "QSNOW", timeidx=time)
    tmp = qc + qi + qs 
    return tmp, z 

def cloud(ff):
    global lev, nk, nlat, nlon
    """ interpolation """
    tot = np.zeros(shape=(nk,nlat,nlon), dtype=np.float32)
    with ProcessPoolExecutor(max_workers=12) as executor:
         for k,output in enumerate(executor.map(interp_var, lev)):
             tot[k,:,:] = output; #print(k)
 
    tot = np.where(np.isnan(tot) , -999.9, tot) 
    """ cloud fraction """
    cldfra = np.zeros(shape=(nk,nlat,nlon), dtype=int) 
    cldfra = np.where(tot > 0.01e-3, 1, 0); del tot
    #print(cldfra.shape) 
    return cldfra
#--------------------------------------------------------------------
start_time = cpu_time()
files = glob.glob("wrfout_d03*") #list
#print(type(files))
nf    = len(files); #print(nf); print(files[0:nf])
#files = np.array(files) #numpy.ndarray
#print(files.shape)
(lev,nk) = delta_z()#dz for interpolation

#a = nc.MFDataset(files)
print("\033[92m\""+os.path.basename(__file__)+"\" is running ...... \033[0m")
for nid,f in enumerate(files[0:nf-1]):
    a = nc.Dataset(f)
    if nid == 0: wrf_tools.Reservoir_loc(a)
    (ntimes, nlev, nlat, nlon, \
     times, lat, lon, res_base) = wrf_tools.info(a, pyngl_map=False)
    print("File #{}: {}".format(nid,f))
    if ntimes < 24:
       print("\033[91mError: Data missing in time. (nt = {})\033[0m".format(ntimes)) 
       exit()
       
    File_out = str(f)[11:24]+"_cloud_char.bin"
    print("\033[43m\033[30mOutput:{}\033[0m".format(File_out))
    newFile  = open(File_out,"wb")
    print("Calculation of Cloud Base/Thickness/Top")
    for t in range(ntimes):
        print("\r\tTime({:02d}/{:02d}): {}".format(t+1,ntimes,times[t]),end="")
        (q_all, z) = get_variable(a,t) #global var.
        cldfra     = cloud(a); #print(cldfra.shape)
        c_thick = np.zeros(shape=(nlat,nlon), dtype=np.float32) 
        c_top   = np.zeros(shape=(nlat,nlon), dtype=np.float32) 
        c_base  = np.zeros(shape=(nlat,nlon), dtype=np.float32) 
        for j in range(nlat):
            for i in range(nlon):
                c_thick[j,i] = max(float( sum(cldfra[:,j,i])-1 )*1.e2, 0.)
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
        #print(str(c_thick[190,1]),str(c_top[190,1]),str(c_base[190,1])) 
        newFile.write(bytearray(c_thick))
        newFile.write(bytearray(c_top))
        newFile.write(bytearray(c_base))
    print("\n"+"-"*80)
    del a #loop t

end_time = cpu_time(); end_time = (end_time - start_time)/60.0
print("\033[92m\"{}\" has done!\nTime elapsed: {:.2f}" \
 .format(os.path.basename(__file__), end_time), "mins.\033[0m")
Ngl.end()
#Nio.end()
#exit()
