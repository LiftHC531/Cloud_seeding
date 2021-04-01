# -*- coding: utf-8 -*-
import sys, os
sys.path.append("/work3/artirain/Cloud_seeding/PyNGL")
#local function
import wrf_tools
from time import time as cpu_time
import numpy as np
import glob
import Ngl; import Nio
import netCDF4 as nc
import xarray as xr
from wrf import getvar,ALL_TIMES,interplevel, g_cape
from numba import jit, vectorize, int32, float32, types
#from joblib import Parallel, delayed

#import metpy.calc as mpcalc
#print(Nio.__version__)
##._FillValue = 9.96920996839e+36

@vectorize("float32(float32,float32,float32,float32,float32)", nopython=True)
def integral_with_vec(a,b,c,d,e):
    return a + (b+c)*(d-e)/np.float32(2) #Trapezoid method

def variables(ff,t):
    """ air density & qc for lwp"""
    z   = getvar(ff, "z", timeidx=t, units="m")
    p   = getvar(ff, "pressure", timeidx=t) #[hPa]
    tk  = getvar(ff, "tk", timeidx=t)       #[K]
    rho = p*100.0/(tk*287.0);del tk, p      #[kg m-3]
    qc  = getvar(ff, "QCLOUD", timeidx=t); qc *= 1.e3 #[g kg-1]
    rhoxqc = rho*qc; del rho, qc
    """ for CTT, Max.dBZ, 500 m W, 500 m Wind vel. """
    dbz = getvar(ff, "REFL_10CM", timeidx=t) #radar reflectivity
    ctt = getvar(ff, "ctt", timeidx=t, units="degC") #Cloud top T
    rh  = getvar(ff, "rh", timeidx=t) #relative humidity %
    wa  = getvar(ff, "wa", timeidx=t, units="m s-1") 
    ua  = getvar(ff, "ua", timeidx=t, units="m s-1") 
    va  = getvar(ff, "va", timeidx=t, units="m s-1") 
    w500 = interplevel(wa, z, 500.0); del wa #500 m height
    u500 = interplevel(ua, z, 500.0); del ua
    v500 = interplevel(va, z, 500.0); del va
    vel_500 = (u500**2+v500**2)**0.5; del u500, v500 
    """ 0:MCAPE, 1:MCIN, 2:LCL, 3:LFC """
    thermo = g_cape.get_2dcape(ff,timeidx=t)
    lcl = thermo[2,:,:] ; del thermo #m 
    #td = getvar(ff, "td", timeidx=t, units="K") #Dew point T
    #print("Get all variables!")
    return rhoxqc, z, dbz, ctt, rh, lcl, w500, vel_500

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
                result[j,i] = integral_with_vec(result[j,i], \
                            rhoxqc[k-1,j,i], rhoxqc[k,j,i],  \
                            z[k,j,i], z[k-1,j,i] ) 
    return result

@jit(float32[:,:](int32, int32, float32[:,:,:]), nopython=True, nogil=True)
def dbz_kernel(ny,nx,dbz):
    result = np.empty(shape=(ny,nx),dtype=np.float32) 
    for j in range(ny):
        for i in range(nx):
            result[j,i] = max(dbz[:,j,i]) 
  
    return result  

@jit(types.UniTuple(float32[:,:],2)\
(int32, int32, int32, float32[:,:,:], float32[:,:,:]) \
, nopython=True, nogil=True)
def RH_kernel(nz,ny,nx,rh,z):
    rh15 = np.empty(shape=(ny,nx),dtype=np.float32) 
    rh50 = np.empty(shape=(ny,nx),dtype=np.float32) 
    rh15[:,:] = 100.0; rh50[:,:] = 100.0 
    for j in range(ny):
        for i in range(nx):
            for k in range(nz):
                if (z[k,j,i] <= 1500.0):
                    if (rh15[j,i] >= rh[k,j,i]):
                        rh15[j,i] = rh[k,j,i]
                elif (z[k,j,i] > 1500.0 and z[k,j,i] <= 5000.0):
                    if (rh50[j,i] >= rh[k,j,i]):
                        rh50[j,i] = rh[k,j,i]

    return rh15, rh50        
#--------------------------------------------------------------------
start_time = cpu_time()
files = glob.glob("wrfout_d03*"); del glob #list
#print(type(files))
nf = len(files); #print(nf); print(files[0:nf])
#files = np.array(files) #numpy.ndarray
#print(files.shape)

#a = nc.MFDataset(files)
print("\033[92m\""+os.path.basename(__file__)+"\" is running ...... \033[0m")
for nid,f in enumerate(files[0:nf-1]):
    a = nc.Dataset(f)
    (ntimes, nlev, nlat, nlon, \
     times, lat, lon, res_base) = wrf_tools.info(a, pyngl_map=False)
    print("File #{}: {}".format(nid,f))
    if ntimes < 24:
       print("\033[91mError: Data missing in time. (nt = {})\033[0m".format(ntimes))
       exit()
    rain_instant = np.float32(rainint(a,ntimes))
    lwp = np.zeros(shape=(ntimes,nlat,nlon),dtype=np.float32)
    mdbz = np.zeros(shape=(ntimes,nlat,nlon),dtype=np.float32)
    rh1500 = np.zeros(shape=(ntimes,nlat,nlon),dtype=np.float32)
    rh5000 = np.zeros(shape=(ntimes,nlat,nlon),dtype=np.float32)
    lcl = np.zeros(shape=(ntimes,nlat,nlon),dtype=np.float32)
    w500 = np.zeros(shape=(ntimes,nlat,nlon),dtype=np.float32)
    ctt = np.zeros(shape=(ntimes,nlat,nlon),dtype=np.float32)
    vel_500 = np.zeros(shape=(ntimes,nlat,nlon),dtype=np.float32)
    for t in range(ntimes):
        print("\r\tTime({:02d}/{:02d}): {}".format(t+1,ntimes,times[t]),end="")
        (rhoxqc, z, dbz, ctt[t,:,:], rh, lcl[t,:,:], \
         w500[t,:,:], vel_500[t,:,:]) = variables(a,t); #print(rhoxqc.dtype)
        lwp[t,:,:] = lwp_kernel(nlev,nlat,nlon,np.float32(rhoxqc),np.float32(z)) 
        mdbz[t,:,:] = dbz_kernel(nlat,nlon,np.float32(dbz))
        (rh1500[t,:,:],rh5000[t,:,:]) = RH_kernel(nlev,nlat,nlon,np.float32(rh),np.float32(z))
    #----------------------------------------------------------------
    File_out = str(f)[11:24]+"_lwplclwdbz.bin"
    print("\n\033[43m\033[30mOutput: {}\033[0m".format(File_out))
    newFile = open(File_out,"wb")
    newFile.write(bytearray(lwp))
    newFile.write(bytearray(lcl))
    newFile.write(bytearray(w500))
    newFile.write(bytearray(mdbz))
    newFile.write(bytearray(ctt))
    newFile.write(bytearray(vel_500))
    newFile.write(bytearray(rh1500))
    newFile.write(bytearray(rh5000))
    newFile.write(bytearray(rain_instant))
    print("-"*80)

end_time = cpu_time(); end_time = (end_time - start_time)/60.0
print("\033[92m\"{}\" has done!\nTime elapsed: {:.2f}" \
 .format(os.path.basename(__file__), end_time), "mins.\033[0m")
Ngl.end()
#Nio.end()
#exit()
