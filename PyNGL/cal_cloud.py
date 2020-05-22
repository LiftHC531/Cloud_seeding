# -*- coding: utf-8 -*-
import wrf_dim_info #local function
from time import time as cpu_time
import numpy as np
import glob
import Ngl; import Nio
import netCDF4 as nc
import xarray as xr
from wrf import getvar,ALL_TIMES,interplevel, g_cape
from numba import jit, vectorize, int32, float32
#from joblib import Parallel, delayed

#import metpy.calc as mpcalc
#print(Nio.__version__)
##._FillValue = 9.96920996839e+36

@vectorize("float32(float32,float32)", nopython=True)
def add_with_vec(a,b):
    return a + b

def variables(ff,t):
    """ air density & qc for lwp"""
    z = getvar(ff, "z", timeidx=t, units="m")
    p = getvar(ff, "pressure", timeidx=t) #hPa
    tk = getvar(ff, "tk", timeidx=t) #K
    rho = p*100.0/tk/287.0;del tk, p #kg m-3
    qc = getvar(ff, "QCLOUD", timeidx=t); qc = qc*10.0**3 #g/kg
    rhoxqc = rho*qc; del rho, qc
    """ for CTT, Max.dBZ, 500 m W, 500 m Wind vel. """
    dbz = getvar(ff, "REFL_10CM", timeidx=t) #radar reflectivity
    """
    ctt = getvar(ff, "ctt", timeidx=t, units="degC") #Cloud top T
    rh = getvar(ff, "rh", timeidx=t) #relative humidity %
    wa = getvar(ff, "wa", timeidx=t, units="m s-1") 
    ua = getvar(ff, "ua", timeidx=t, units="m s-1") 
    va = getvar(ff, "va", timeidx=t, units="m s-1") 
    w500 = interplevel(wa, z, 500.0); del wa #500 m height
    u500 = interplevel(ua, z, 500.0); del ua
    v500 = interplevel(va, z, 500.0); del va
    vel_500 = (u500**2+v500**2)**0.5; del u500, v500 """
    """ 0:MCAPE, 1:MCIN, 2:LCL, 3:LFC """
    thermo = g_cape.get_2dcape(ff,timeidx=t)
    lcl = thermo[2,:,:] ; del thermo
    #td = getvar(ff, "td", timeidx=t, units="K") #Dew point T
    #print("Get all variables!")
    return rhoxqc, z,dbz

def rainint(ff,nt):
    def calc_rain(rain_exp,nt):
        result = rain_exp; #print(result.shape)
        for t in range(1,nt):
            result[t,:,:] = rain_exp[t,:,:] - rain_exp[t-1,:,:]

        result[0,:,:] = result[1,:,:] 
        return result
 
    """ Rainfall rate """
    rain_exp = getvar(ff, "RAINNC", timeidx=ALL_TIMES) #-1 isn't work
    rain_con = getvar(ff, "RAINC", timeidx=ALL_TIMES)
    rain_exp = rain_exp + rain_con; del rain_con #total rainfall
    rain_instant = calc_rain(rain_exp,nt)
    #print("Get Rainfall rate!")
    return rain_instant

@jit(float32[:,:](int32, int32, int32, float32[:,:,:], float32[:,:,:]), nopython=True, nogil=True)
def lwp_kernel(nz,ny,nx,rhoxqc,z):
    #def lwp_kernel(nz: int,ny: int,nx: int,rhoxqc: xr.core.dataarray.DataArray,z: xr.core.dataarray.DataArray) -> np.ndarray: 
    result = np.zeros(shape=(ny,nx),dtype=np.float32) #np.empty
    #print(result.dtype)
    for j in range(ny):
        for i in range(nx):
            for k in range(1,nz):
                #result[j,i] += (rhoxqc[k-1,j,i] + rhoxqc[k,j,i])/2.0 \
                #       *(z[k,j,i]-z[k-1,j,i])
                result[j,i] = add_with_vec(result[j,i], \
                    (rhoxqc[k-1,j,i] + rhoxqc[k,j,i])/np.float32(2) \
                              *(z[k,j,i]-z[k-1,j,i]) ) 
    return result

@jit(float32[:,:](int32, int32, float32[:,:,:]), nopython=True, nogil=True)
def dbz_kernel(ny,nx,dbz):
    result = np.empty(shape=(ny,nx),dtype=np.float32) 
    for j in range(ny):
        for i in range(nx):
            result[j,i] = max(dbz[:,j,i]) 
  
    return result          
#--------------------------------------------------------------------
start_time = cpu_time()
files = glob.glob("wrfout_d03*"); del glob #list
#print(type(files))
nf = len(files); #print(nf); print(files[0:nf])
#files = np.array(files) #numpy.ndarray
#print(files.shape)

#a = nc.MFDataset(files)
for nid,f in enumerate(files[0:1]):#nf-1]):
    File_out = str(f)[11:24]+"_lwplclwdbz.bin"; print("output:",File_out)
    newFile = open(File_out,"wb")
    a = nc.Dataset(f)
    (ntimes, nlev, nlat, nlon, times) = wrf_dim_info.info(a)
    print("File:",nid,f,"\nDomain:{}".format(ntimes)+"|{}".format(nlev)+\
    "|{}".format(nlat)+"|{}".format(nlon)," .......Get Dimensions!!")
    rain_instant = rainint(a,ntimes)
    lwp = np.zeros(shape=(ntimes,nlat,nlon),dtype=np.float32)
    mdbz = np.zeros(shape=(ntimes,nlat,nlon),dtype=np.float32)
    for t in range(ntimes):
        print(str(times[t].values)[0:19])
        (rhoxqc, z, dbz) = variables(a,t); #print(rhoxqc.dtype)
        lwp[t,:,:] = lwp_kernel(nlev,nlat,nlon,np.float32(rhoxqc),np.float32(z)) 
        mdbz[t,:,:] = dbz_kernel(nlat,nlon,np.float32(dbz))

end_time = cpu_time(); end_time = (end_time - start_time)/60.0
print("cal_cloud.py has done!\nTime elapsed: {:.2f}".format(end_time), "mins.")
Ngl.end()
#Nio.end()
#exit()
