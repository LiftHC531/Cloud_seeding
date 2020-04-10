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
print(max(hgt))
print(times)
setfileoption("bin","WriteByteOrder","LittleEndian")
File_out = "cloud_char.bin"
system("rm "+File_out)
;----------------------------------
do t = 0,ntimes-1
  print(times(t))
  qc = wrf_user_getvar(a,"QCLOUD",t)
  qi = wrf_user_getvar(a,"QICE",t)
  qs = wrf_user_getvar(a,"QSNOW",t)
  ;printVarSummary(qc)
  tot = qc + qi + qs 

  nk = 122
  lev = new((/nk/),float)
  tot_cloud = new((/nk,nlat,nlon/),float)
  z  = wrf_user_getvar(a,"z",t)    ; z on mass points
  kk = 900.0
  do k= 0,nk-1
   ;print(kk)
   lev(k) = kk
   tot_cloud(k,:,:) = wrf_user_intrp3d(tot,z,"h",kk,0.,False)
   kk = kk + 100.0
  end do
  ;lev = lev*10.0^(-3) ; m-->km
  ;printVarSummary(tot_cloud)
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cldfra = new((/nk,nlat,nlon/),integer)
  cldfra = 0
  do k = 0,nk-1
   do j = 0,nlat-1
    do i = 0,nlon-1
      if (.not.ismissing(tot_cloud(k,j,i)) .and. \
          tot_cloud(k,j,i).gt.10.0^(-5)) then 
        cldfra(k,j,i) = 1      
      end if
    end do
   end do
  end do

;print(lev(:)+"km,"+cldfra(:,90,165))
;print(lev(:)+"km,"+cldfra(:,148,35))
  c_base = new((/nlat,nlon/),float)
  c_top  = new((/nlat,nlon/),float)
  c_base = 0.0
  c_top  = 0.0
  do j = 0,nlat-1
   do i = 0,nlon-1
    do k = 1,nk-1
       sp1 = 0.0
       if(cldfra(k,j,i).eq.1.and.cldfra(k-1,j,i).eq.0.and.sp1.eq.0.0) then
         c_base(j,i) = lev(k)
         sp1 = 1.0
       end if
       if (k.lt.nk-1) then
         if(cldfra(k+1,j,i).eq.0.and.cldfra(k,j,i).eq.1) then
         c_top(j,i) = lev(k)
         end if
       else if (k.eq.nk-1.and.cldfra(k,j,i).eq.1.) then
         c_top(j,i) = lev(k)
       end if
       end if
    end do
    ;sp2 = 0.0
    ;do k = nk-2,0,-1
    ;   if(cldfra(k,j,i).eq.1.and.cldfra(k+1,j,i).eq.0.and.sp2.eq.0.0) then
    ;     c_top(j,i) = lev(k)
    ;     sp2 = 1.0
    ;   end if
    ;end do
   end do
  end do
;print(c_base(148,35))
;print(c_top(148,35))
  c_thick = new((/nlat,nlon/),float)
  c_thick = c_top - c_base


   fbindirwrite(File_out,c_thick(:,:))
   fbindirwrite(File_out,c_top(:,:))
   fbindirwrite(File_out,c_base(:,:))
delete(lev)
delete(tot_cloud)
delete(cldfra)
delete(c_thick)
delete(c_top)
delete(c_base)
end do

;------------------------------------------------------------------------
end


