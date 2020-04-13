load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;load "/home/mlhchen/big_WRFUserARW.ncl"

begin
;------------------------------------------------------------------------
  wrf_path = "./"
  wrf_cont = "wrfout_d03_2020-03-10_000000.nc"
  wrf_Name = wrf_path + wrf_cont
  a = addfile( wrf_Name , "r")
  times    = wrf_user_getvar(a,"times",-1)
  ntimes   = dimsizes(times) ; number of times in the files
  w_noplt  = wrf_user_getvar(a,"wa",0)
  ;printVarSummary(w_noplt)
  nlev     = 38  ;dimsizes(w_noplt(:,0,0))
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
 ctt = wrf_user_getvar(a,"ctt",-1)
 printVarSummary(ctt)

;----------------------------------
if (True) then
wks_type = "png"
wks_type@wkWidth = 3000 ;page resolution
wks_type@wkHeight = 3000
do t = 1,1
print(t)
dtime = str_split_by_length(times(t), (/5,2,1,2,1,2,1,2,3/) )
print(dtime(1)+dtime(3)+dtime(5)+dtime(7))
imagename = dtime(1)+dtime(3)+dtime(5)+dtime(7)+"_CTT"
wks = gsn_open_wks(wks_type,imagename)

gsn_define_colormap(wks,"precip2_15lev")

opts = True
opts@NoHeaderFooter = True
opts@FieldTitle = ""
;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
res = opts
res@cnFillOn = True
tres = res
tres2 = tres
;res@UnitLabel = "dBZ"
;res@cnLevelSelectionMode = "ManualLevels"
;res@cnMinLevelValF        = 0.0
;res@cnMaxLevelValF        =70.0
;res@cnLevelSpacingF       = 5.0
res@lbLabelBarOn = False ;Turn off individual labelbars
res@lbTitleString = "1/s"
res@cnLineLabelFormat="@*+^sg"
 ;res@cnFillOpacityF = 0.7
;XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
tres2@cnFillOn = False
tres2@cnLinesOn = True        ; turn on contour lines
tres2@cnLineLabelsOn = False
tres2@cnLineThicknessF = 20
;tres2@cnLineLabelFontHeightF = 0.02
tres2@cnInfoLabelOn  = False  ; no info label
tres2@cnLevelSelectionMode = "ExplicitLevels"
tres2@cnLevels = (/1000/)
tres2@cnLineColor = "gray40"


;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pltres = True
pltres@PanelPlot = True
;>>>>>>>>>>>>>>
mpres = True
mpres@gsnMaximize = True
  mpres@mpFillOn              = False
  mpres@mpOutlineOn           = False
mpres@mpGeophysicalLineColor = "Black"
mpres@mpGeophysicalLineThicknessF = 25.0
mpres@mpGridAndLimbOn = False
mpres@mpDataBaseVersion = "MediumRes"
mpres@tmXBLabelFontHeightF = 0.015 ;the font size~~~~~~~~~~~~~
mpres@tmYLLabelFontHeightF = 0.015
mpres@tmBorderThicknessF = 25.0
mpres@tmXBLabelFont = 26
mpres@tmXBMajorThicknessF = 25
mpres@tmYLLabelFont = 26
mpres@tmYLMajorThicknessF = 25
;>>>>>>>>>>>>>>
res_pn = True
res_pn@gsnMaximze = True
res_pn@txString ="~F25~             Wind~B~850~N~"
;res_pn@gsnPanelYWhiteSpacePercent = 10 ;Add white space bar/wks plots
res_pn@gsnPanelLabelBar = True        ;Turn on common labelbar
;res_pn@lbBoxSeparatorLinesOn= False
res_pn@lbTitleOn  =  True
res_pn@lbTitleString = "~F25~~S~o~N~C"
res_pn@lbTitlePosition  = "Top"
res_pn@lbLabelFont          = 26
;res_pn@lbLabelAlignment = "BoxCenters"         ; label orientation
res_pn@pmLabelBarHeightF      = 0.06;0.07
res_pn@pmLabelBarWidthF      = 0.65
res_pn@pmLabelBarOrthogonalPosF = -0.02
res_pn@pmLabelBarParallelPosF = 0.02
res_pn@lbTitleFontHeightF= 0.015               ; make title smaller
res_pn@lbTitleDirection = "Across"             ; title direction
;>>>>>>>>>>>>>>
contour = new(61,graphic)
plot_pn = new(72,graphic)
hh = t/2.0
;res@FieldTitle = "~F26~Valid at "+dtime(5)+dtime(7)+" UTC("+hh+"-hr) "\
;                 +dtime(3)+"~S~th~N~ Jun., 2008"
res@FieldTitle = "~F26~Valid at "+dtime(5)+dtime(7)+" UTC"\
                 +dtime(3)+"~S~th~N~ Mar., 2020"
hgt@description = ""
hgt@units = ""
ctt@description = ""
ctt@units = ""
;==============
contour(1) = wrf_contour(a,wks,ctt(t,:,:),res)
;==============
 contour(2) = wrf_contour(a,wks,hgt(:,:),tres2)

plot_pn(0) = wrf_map_overlays(a,wks,(/contour(1),contour(2)/),pltres,mpres)


if (False) then
 polyres               = True          ; poly marker mods desired
 polyres@gsMarkerIndex =  12;4           ; choose circle as polymarker
 polyres@gsMarkerSizeF = 20.0           ; select size to avoid streaking
 polyres@gsMarkerThicknessF = 10;30.0
; polyres@gsMarkerOpacityF = 0.3
 polyres@gsMarkerColor = (/"purple4"/)   ; choose color

ic = 165
jc = 90
  xd = lon(jc,ic)
  yd = lat(jc,ic)
dumii = gsn_add_polymarker(wks,plot_pn(0),xd,yd,polyres)  ; draw polymarkers

end if

;---Attach shapefile polylines to map
  dir       = "/home/mlhchen/shapefile/"
  filename1  = "TWN_adm0.shp"
  ;filename2  = "TWN_water_lines_dcw.dbf"
  filename3  = "CHN_adm0.shp"
  lnres                  = True
  lnres@gsLineThicknessF = 10.0
  lnres1                  = True
  lnres1@gsLineThicknessF = 3.0
  lnres1@gsLineColor = "blue"
  poly1 = gsn_add_shapefile_polylines(wks,plot_pn(0),dir+filename1,lnres)
  ;poly2 = gsn_add_shapefile_polylines(wks,plot_pn(0),dir+filename2,lnres1)
  poly3 = gsn_add_shapefile_polylines(wks,plot_pn(0),dir+filename3,lnres)

  gsn_panel(wks,plot_pn(0),(/1,1/),res_pn)

end do ; TIME


end if
;------------------------------------------------------------------------
end

