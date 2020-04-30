load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin
;-------------------------------------------------------------------------
  wrf_path = "./"
  wrf_cont = "wrfout_d03_2020-04-29_00:00:00.nc"
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
delete(w_noplt)
;----------------------------------
print(lon(75,0)+", "+lat(0,75))
dx = lon(75,76)-lon(75,75)
dy = lat(76,75)-lat(75,75)
print(dx+", "+dy)
ic = 88;70;75-1;85
jc = 39;30;33-1;55
print(lon(jc,ic)+", "+lat(jc,ic))

;pt_lon = 121.1929121  ;shimen
;pt_lat =  24.8153923
;pt_lon = 120.5361886  ;zengwen 
;pt_lat =  23.2478001
;pt_lon = 120.5374418  ;nanhua
;pt_lat =  23.089879
;pt_lon = 121.0555066   ;baoshen
;pt_lat =  24.7278826
pt_lon = 120.4323562   ;kaoping
pt_lat =  22.6925145
       print("pt:lon="+pt_lon+", lat="+pt_lat)
  do j = 0,nlat-1
   do i = 0,nlon-1
      diffi = abs(lon(j,i)-pt_lon)
      diffj = abs(lat(j,i)-pt_lat)
      if ( (diffi+diffj) .le. 2.8*10^-2) then
         pt_i = i
         pt_j = j
         print("lon="+lon(j,i)+", lat="+lat(j,i)+", i="+i+", j="+j)
      end if
   end do
  end do

;print(lon(45,87)+", "+lat(45,87))


;------------------------------------------------------------------------
end

