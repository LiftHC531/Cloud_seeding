load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;load "/home/mlhchen/big_WRFUserARW.ncl"

begin
;-------------------------------------------------------------------------
 files = systemfunc("ls -d ./wrfout_d03_*") + ".nc"
  member = dimsizes(files)
  a = addfiles( files , "r")
lat_c = (/136,132,80,74,60/)  ; y
lon_c = (/ 93, 88,71,71,68/)  ; x
ndam = dimsizes(lat_c)

do mm = 0,member-1-1        ;not ouput last file
   date = stringtochar(files(mm))
   print(date(13:25))
   times    = wrf_user_getvar(a[mm],"times",-1)
   ntimes   = dimsizes(times) ; number of times in the files
   w_noplt  = wrf_user_getvar(a[mm],"wa",0)
   ;printVarSummary(w_noplt)
   nlev     = dimsizes(w_noplt(:,0,0))
   nlat     = dimsizes(w_noplt(0,:,0))
   nlon     = dimsizes(w_noplt(0,0,:))
   lon = wrf_user_getvar(a[mm],"XLONG",0)
   lat = wrf_user_getvar(a[mm],"XLAT",0)
print("Domain:"+ntimes+"|"+nlev+"|"+nlat+"|"+nlon+" .......Get Dimensions!!")
hgt = wrf_user_getvar(a[mm],"HGT",0)  ;hgt = a->HGT(0,:,:)
print(max(hgt))
print(times)
;----------------------------------
 UNDEF =  9.96921e+36                     ; UNDEF
rain_exp = new( (/ntimes,nlat,nlon/), float,UNDEF)
rain_con = new( (/ntimes,nlat,nlon/), float,UNDEF)
rain_tot = new( (/ntimes,nlat,nlon/), float,UNDEF)
  rain_exp = wrf_user_getvar(a[mm],"RAINNC",-1)
  rain_con = wrf_user_getvar(a[mm],"RAINC",-1)
  rain_tot = rain_exp + rain_con

  rain_int = new ( (/ndam,ntimes/), float, UNDEF)
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (date(24:25).eq."12") then
    ti = 12
    tf = 23
    day = sprinti("%2.2i", stringtointeger(date(21:22)) + 1)
else
    ti = 0
    tf = 12
    day = sprinti("%2.2i", stringtointeger(date(21:22)) )
end if
print("from: "+times(ti))
print("  to: "+times(tf))
do t = 1,ntimes-1
  do dd = 0,ndam-1
     rain_int(dd,t) =  max(rain_tot(t,lat_c(dd)-1:lat_c(dd)+1,lon_c(dd)-1:lon_c(dd)+1)) - \
                       max(rain_tot(t-1,lat_c(dd)-1:lat_c(dd)+1,lon_c(dd)-1:lon_c(dd)+1))
  end do
end do
rain_int(:,0) = rain_int(:,1)
;=============================================================================
if (True) then
wks_type = "png"
wks_type@wkWidth = 3000 ;page resolution
wks_type@wkHeight = 3000
  imagename = date(13:20)+day+"_QPF_intensity"
;=================
  wks = gsn_open_wks(wks_type,imagename)
;=================
opts = True
opts@gsnMaximize = True
;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
res = opts
res@gsnDraw = False
res@gsnFrame = False
res@vpHeightF= 0.4           ; change aspect ratio of plot
res@vpWidthF = 0.8
;res@xyLineColor = "black"
;res@trYMinF            =  0   ; Limits for Y axis.  The limits
;res@trYMaxF            =  5
res@trXMinF            = ti    ; Limits for X axis.  The limits
res@trXMaxF            = tf
res@tmXBMode   = "Explicit"
res@tmXBValues   = ispan(ti,tf,2)  ;ispan(0,20,4)
;res@tmXBLabels   = (/"0800","1200","1600","2000","0000","0400"/)
res@tmXBLabels   = (/"0800","1000","1200","1400","1600","1800","2000"/)
res@tmXBMinorValues   = ispan(0,23,1)
;res@tmYLMode   = "Explicit"
;res@tmYLValues   = ispan(1,5,2)
;res@tmYLLabels   = (/"1","3","5"/)
;res@tmYLMinorValues   = ispan(0,5,1)
res@tmBorderThicknessF = 10.0
res@tmXBMajorThicknessF = 10.0
res@tmXBMinorThicknessF = 10.0
res@tmYLMajorThicknessF = 10.0
res@tmYLMinorThicknessF = 10.0
res@tmXBLabelFont = 26
res@tmXBLabelFontHeightF = 0.017
res@tmYLLabelFont = 26
res@tmYLLabelFontHeightF = 0.017
res@tiXAxisString = "~F25~Time(LST)"
res@tiYAxisString = "~F25~Rain Intensity"
res@tiMainString = "~F25~109"+date(17:20)+day
res2 = res
res2@xyLineThicknesses = 30
  res@pmLegendDisplayMode    = "Always"            ; turn on legend

  res@pmLegendSide           = "Right"               ; Change location of
  res@pmLegendParallelPosF   = 0.8;.45    ;+ top - bottom
  res@pmLegendOrthogonalPosF = -0.3;-0.4  ;+ right - left

  res@pmLegendWidthF         = 0.10              ; Change width and
  res@pmLegendHeightF        = 0.15                ; height of legend.
  res@lgLabelFontHeightF     = 0.017                 ; change font height
  res@lgTitleOn              = True                ; turn on legend title
  res@lgTitleString          = "Reservoir"         ; create legend title
  res@lgTitleFontHeightF     = .020                ; font of legend title
  res@xyExplicitLegendLabels = (/"  Shi-Men","  Bao-Shan II","  Zeng-Wen","  Nan-Hua","  Kao-Ping"/)         ; explicit labels
  res@xyLineColors = (/"black","orange","red","blue","purple"/)
  res@xyLineThicknesses = (/30,30,30,30,30/)
  res@xyDashPatterns = (/0,0,0,0,0/)
;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
res_pn = True
;>>>>>>>>>>>>>>
plot_pn = new(20,graphic)
;dot_pn = new(ntime,graphic)
plot_pn(0) = gsn_csm_y (wks,rain_int(:,:),res)
res2@xyLineColor = "orange"
res2@xyDashPattern = 14
;plot_pn(1) = gsn_csm_y (wks,std,res2)

;overlay(plot_pn(0),plot_pn(1))
draw (plot_pn(0))
frame(wks)

  delete(rain_exp)
  delete(rain_con)
  delete(rain_tot)
  delete(rain_int)
  delete(res)
  delete(res2)
  delete(plot_pn)
end if

end do

;------------------------------------------------------------------------
end



