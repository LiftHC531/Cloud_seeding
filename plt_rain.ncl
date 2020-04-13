load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;load "/home/mlhchen/big_WRFUserARW.ncl"

begin
;-------------------------------------------------------------------------
 files = systemfunc("ls -d ../wrfout/wrfout_d03_*") + ".nc"
  member = dimsizes(files)
  a = addfiles( files , "r")
do mm = 0,member-1
  date = stringtochar(files(mm))
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
;a_lon = 121.07
;a_lat = 24.73
;do j = 0,nlat-1
; do i = 0,nlon-1
;    dis = sqrt( (lon(j,i)-a_lon)^2+(lat(j,i)-a_lat)^2 )
;    if (dis.le.10^-2) then
;       print("i="+i+", j="+j+",dis="+dis)
;    end if
; end do
;end do
;----------------------------------
 UNDEF =  9.96921e+36                     ; UNDEF
rain_exp = new( (/ntimes,nlat,nlon/), float,UNDEF)
rain_con = new( (/ntimes,nlat,nlon/), float,UNDEF)
rain_tot = new( (/ntimes,nlat,nlon/), float,UNDEF)
  rain_exp = wrf_user_getvar(a[mm],"RAINNC",-1)
  rain_con = wrf_user_getvar(a[mm],"RAINC",-1)
  rain_tot = rain_exp + rain_con


  rain_12h = new ( (/nlat,nlon/), float, UNDEF)
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ti =  0;11 
  tf =  12;23
print("from: "+times(ti))
print("to: "+times(tf))
     rain_12h(:,:) =  rain_tot(tf,:,:) - rain_tot(ti,:,:)
;=============================================================================
if (True) then
wks_type = "png"
wks_type@wkWidth = 3000 ;page resolution
wks_type@wkHeight = 3000
  imagename = date(21:33)+"_12hr_QPF"
;=================
  wks = gsn_open_wks(wks_type,imagename)
;=================
colors = (/"white","black",\
 "white","(/0.608,1,1/)","(/0.,0.812,1.0/)","(/0.039,0.596,1.0/)",\
 "(/0.039,0.396,1.0/)",\ ;---blue array---
 "(/0.188,0.6,0.039/)","(/0.196,1.0,0.0/)",\ ;---green array---
"(/0.973,1.0,0.0/)","(/1.0,0.796,0.0/)","(/1.0,0.603,0.0/)",\ ;-yellow array-
"(/0.98,0.012,0.0/)","(/0.8,0.0,0.012/)","(/0.627,0.0,0.0/)",\;--red array--
"(/0.596,0.,0.604/)","(/0.765,0.016,0.8/)",\
"(/0.973,0.02,0.953/)","(/0.996,0.796,1.0/)"/)

gsn_define_colormap(wks,colors)
opts = True
opts@NoHeaderFooter = True
tres = opts
;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
res = opts
res@cnFillOn = True
;res@UnitLabel = "mm"
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = (/ 1,2,6,10,15,20,30,40,50,70,90,110,130,150,200,300 /)
res@lbLabelBarOn = False ;Turn off individual labelbars
 res@cnFillOpacityF = 0.8

;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pltres = True
pltres@PanelPlot = True
;>>>>>>>>>>>>>>
mpres = True
  mpres@mpFillOn              = False
  mpres@mpOutlineOn           = False
mpres@gsnMaximize = True
mpres@mpGeophysicalLineColor = "Black"
mpres@mpGeophysicalLineThicknessF = 25.0
mpres@mpGridAndLimbOn = False
mpres@mpDataBaseVersion = "MediumRes"
mpres@tmXBLabelFontHeightF = 0.015;0.02 ;the font size~~~~~~~~~~~~~
mpres@tmYLLabelFontHeightF = 0.015;0.02
mpres@tmBorderThicknessF = 25.0
mpres@tmXBLabelFont = 26
mpres@tmXBMajorThicknessF = 25
mpres@tmYLLabelFont = 26
mpres@tmYLMajorThicknessF = 25
;>>>>>>>>>>>>>>
res_pn = True
res_pn@gsnMaximze = True
res_pn@txString ="~F26~                "+date(21:30)
;res_pn@gsnPanelYWhiteSpacePercent = 10 ;Add white space bar/wks plots
res_pn@gsnPanelLabelBar = True        ;Turn on common labelbar
res_pn@lbBoxSeparatorLinesOn= False
res_pn@lbTitleOn  =  True
res_pn@lbTitleString = "~F25~mm"
res_pn@lbTitlePosition  = "Top"
res_pn@lbLabelFont          = 26
res_pn@lbLabelFontHeightF= 0.013
;res_pn@lbLabelAlignment = "BoxCenters"         ; label orientation
res_pn@pmLabelBarHeightF      = 0.06
res_pn@pmLabelBarWidthF      = 0.65
res_pn@pmLabelBarOrthogonalPosF = -0.02
res_pn@pmLabelBarParallelPosF = 0.035
res_pn@lbTitleFontHeightF= 0.02               ; make title smaller
res_pn@lbTitleDirection = "Across"             ; title direction
;>>>>>>>>>>>>>>
;res@FieldTitle = "from "+times(ti)+" to "+times(tf)
res@FieldTitle = "~F26~12-hr accumulated precipitation"
contour = wrf_contour(a[mm],wks,rain_12h(:,:),res)
plot_pn =  wrf_map_overlays(a[mm],wks,(/contour/),pltres,mpres)


if (True) then
 polyres               = True          ; poly marker mods desired
 polyres@gsMarkerIndex =  12;4           ; choose circle as polymarker
 polyres@gsMarkerSizeF = 20.0           ; select size to avoid streaking
 polyres@gsMarkerThicknessF = 10 
; polyres@gsMarkerOpacityF = 0.3

dumii = new(3,graphic)
;ic = (/ 93,71,71/) 
;jc = (/136,80,74/)
xd = new((/3/),float)
yd = new((/3/),float)
xd = (/121.1929121,120.5361886,120.5374418/)
yd = (/24.8153923,23.2478001,23.089879/)
markcolor = (/"black","red","blue"/)
do rr = 0,3-1
   ;xd(rr) = lon(jc(rr),ic(rr))
   ;yd(rr) = lat(jc(rr),ic(rr))
 polyres@gsMarkerColor = markcolor(rr)   ; choose color
 dumii(rr) = gsn_add_polymarker(wks,plot_pn(0),xd(rr),yd(rr),polyres)  ; draw polymarkers
end do
delete(xd)
delete(yd)
end if

if (False) then
 polyres@gsMarkerColor = (/"blue"/)   ; choose color

ic = 60-1
jc = 25-1
  xd = lon(jc,ic)
  yd = lat(jc,ic)
dumii2 = gsn_add_polymarker(wks,plot_pn(0),xd,yd,polyres)  ; draw polymarkers

end if


;---Attach shapefile polylines to map
  dir1       = "/home/WRF/shapefile/Shimen/"
  filename1  = "COUNTY_MOI_1080617.shp";"TWN_adm2.shp"
  ;dir2       = "/home/mlhchen/shapefile/"
  ;filename2  = "TWN_water_lines_dcw.dbf"
  dir3       = "/home/WRF/shapefile/" ;"/home/mlhchen/shapefile/"
  filename3  = "CHN_adm0.shp"
  lnres                  = True
  lnres@gsLineThicknessF = 10.0
  lnres1                  = True
  lnres1@gsLineThicknessF = 3.0
  lnres1@gsLineColor = "blue"
  poly1 = gsn_add_shapefile_polylines(wks,plot_pn(0),dir1+filename1,lnres)
  ;poly2 = gsn_add_shapefile_polylines(wks,plot_pn(0),dir2+filename2,lnres1)
  poly3 = gsn_add_shapefile_polylines(wks,plot_pn(0),dir3+filename3,lnres)


end if
  gsn_panel(wks,plot_pn(0),(/1,1/),res_pn)

  delete(rain_exp)
  delete(rain_con)
  delete(rain_tot)
  delete(rain_12h)
  delete(res)
  delete(pltres)
  delete(mpres)
end do
;------------------------------------------------------------------------
end



