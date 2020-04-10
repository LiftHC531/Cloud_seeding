load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "/home/mlhchen/big_WRFUserARW.ncl"

begin
;-------------------------------------------------------------------------
  wrf_path = "./"
  wrf_cont = "wrfout_d03_2020-03-10_000000.nc"
  wrf_Name = wrf_path + wrf_cont
  a = addfile( wrf_Name , "r")
  times    = wrf_user_getvar(a,"times",-1)
  ntimes   = dimsizes(times) ; number of times in the files
  w_noplt  = wrf_user_getvar(a,"wa",0)
  ;printVarSummary(w_noplt)
  nlev     = dimsizes(w_noplt(:,0,0))
  nlat     = dimsizes(w_noplt(0,:,0))
  nlon     = dimsizes(w_noplt(0,0,:))
  lon = wrf_user_getvar(a,"XLONG",0)
  lat = wrf_user_getvar(a,"XLAT",0)
print("Domain:"+ntimes+"|"+nlev+"|"+nlat+"|"+nlon+" .......Get Dimensions!!")
hgt = wrf_user_getvar(a,"HGT",0)  ;hgt = a->HGT(0,:,:)
z = wrf_user_getvar(a,"z",-1)
;print(z(1,:,120,123))
print(z(1,:,149,51))
print(max(hgt))
print(times)
;----------------------------------
 p = wrf_user_getvar(a,"pressure",-1) ; hPa
 ;printVarSummary(p)
 tk = wrf_user_getvar(a,"tk",-1)
 rho =  p*100.0/tk/287.0 ; kg/m3
 print(max(rho(:,:,:,:)))
 qc = wrf_user_getvar(a,"QCLOUD",-1)
 qc = qc*1.E3 ; g/kg
;------
 ctt = wrf_user_getvar(a,"ctt",-1)   
;------
 td = wrf_user_getvar(a,"td",-1)
 td = td + 273.15
;------
 dbz = wrf_user_getvar(a,"REFL_10CM",-1)
 mdbz = dbz(:,0,:,:)
;------
  wa  = wrf_user_getvar(a,"wa",-1)
  w1000 = wrf_user_intrp3d(wa,z,"h",1000.0,0.,False)
  ua  = wrf_user_getvar(a,"ua",-1)
  va  = wrf_user_getvar(a,"va",-1)
  u1000 = wrf_user_intrp3d(ua,z,"h",1000.0,0.,False)
  v1000 = wrf_user_intrp3d(va,z,"h",1000.0,0.,False)
  vel_1000 = sqrt(u1000^2.0+v1000^2.0)
;------
  qv = wrf_user_getvar(a,"QVAPOR",-1);kg/kg
  rh = relhum(tk,qv,p*100.0)
;------
  rain_exp = wrf_user_getvar(a,"RAINNC",-1)
  rain_con = wrf_user_getvar(a,"RAINC",-1)
  rain_tot = rain_exp + rain_con
  rain_instant = rain_tot
do t = 1,ntimes-1
   rain_instant(t,:,:) =  rain_tot(t,:,:) - rain_tot(t-1,:,:)
end do
rain_instant(0,:,:) = rain_instant(1,:,:)

;----------------------------------
 lwp = new((/ntimes,nlat,nlon/),float)
 lwp = 0.0
 lcl = new((/2,ntimes,nlat,nlon/),float)
 lcl = 0.0
 rh1500 = rh(:,0,:,:) 
 rh5000 = new((/ntimes,nlat,nlon/),float)
 rh5000 = 100.0
 do t = 0,ntimes-1
   do j = 0,nlat-1
     do i = 0,nlon-1
           mdbz(t,j,i) = max(dbz(t,:,j,i))
           temp = 0.0
           lcl(0,t,j,i) = lclvl(p(t,0,j,i),tk(t,0,j,i),td(t,0,j,i))
           sp = 0.0
        do k = 1,nlev-1
           temp = temp + (rho(t,k-1,j,i)*qc(t,k-1,j,i)+rho(t,k,j,i)*qc(t,k,j,i))/2.0 \
                 *(z(t,k,j,i)-z(t,k-1,j,i))
           if (p(t,k,j,i).le.lcl(0,t,j,i).and.sp.eq.0.0) then
               lcl(1,t,j,i) = z(t,k-1,j,i) + (z(t,k,j,i)-z(t,k-1,j,i))*(lcl(0,t,j,i)-p(t,k-1,j,i))/(p(t,k,j,i)-p(t,k-1,j,i))
               sp = 1.0
           end if
          ;=================
           if (z(t,k,j,i).le.1500.and.rh1500(t,j,i).ge.rh(t,k,j,i)) then
              rh1500(t,j,i) = rh(t,k,j,i)
           end if
           if (z(t,k,j,i).gt.1500.and.z(t,k,j,i).le.5000.and.rh5000(t,j,i).ge.rh(t,k,j,i)) then
              rh5000(t,j,i) = rh(t,k,j,i)
           end if
        end do
        lwp(t,j,i) = temp 
     end do
   end do
 end do
 ;print(max(lwp(:,:,:)))

;----------------------------------
 File_out = "lwplclwdbz.bin"
 system("rm "+File_out)
 setfileoption("bin","WriteByteOrder","LittleEndian")
 fbindirwrite(File_out,lwp)
 fbindirwrite(File_out,lcl(1,:,:,:))
 fbindirwrite(File_out,w1000)
 fbindirwrite(File_out,mdbz)
 fbindirwrite(File_out,ctt)
 fbindirwrite(File_out,vel_1000)
 fbindirwrite(File_out,rh1500)
 fbindirwrite(File_out,rh5000)
 fbindirwrite(File_out,rain_instant)

;------------------------------------------------------------------------
end

