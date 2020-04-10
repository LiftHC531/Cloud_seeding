load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;load "/home/mlhchen/big_WRFUserARW.ncl"

begin
;-------------------------------------------------------------------------
setfileoption("bin","ReadByteOrder","LittleEndian")
UNDEF = 9.96921e+36
file_a = "./2020-03-10_score.bin"
nt = 24
ndam = 3
score = new ( (/ndam,nt/), float, UNDEF)
score = fbindirread(file_a, 0, (/ndam,nt/), "float")
std = score
std = 3.0
date = stringtochar(file_a)
print(date(2:11))
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (True) then
wks_type = "png"
wks_type@wkWidth = 3000 ;page resolution
wks_type@wkHeight = 3000

imagename = date(2:11)+"score_series"
wks = gsn_open_wks(wks_type,imagename)
opts = True
opts@gsnMaximize = True
;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
res = opts
res@gsnDraw = False
res@gsnFrame = False
res@vpHeightF= 0.4           ; change aspect ratio of plot
res@vpWidthF = 0.8
;res@xyLineColor = "black"
res@trYMinF            =  0   ; Limits for Y axis.  The limits
res@trYMaxF            =  5
res@trXMinF            =  0    ; Limits for X axis.  The limits
res@trXMaxF            = 12;23 
res@tmXBMode   = "Explicit"
res@tmXBValues   = ispan(0,12,2);ispan(0,20,4)
;res@tmXBLabels   = (/"0800","1200","1600","2000","0000","0400"/)
res@tmXBLabels   = (/"0800","1000","1200","1400","1600","1800","2000"/)
res@tmXBMinorValues   = ispan(0,23,1)
res@tmYLMode   = "Explicit"
res@tmYLValues   = ispan(1,5,2)
res@tmYLLabels   = (/"1","3","5"/)
res@tmYLMinorValues   = ispan(0,5,1)
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
res@tiYAxisString = "~F25~Score"
res@tiMainString = "~F25~"+date(2:11)
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
  res@xyExplicitLegendLabels = (/"  Shi-Men","  Zeng-Wen","  Nun-Hua"/)         ; explicit labels
  res@xyLineColors = (/"black","red","blue"/)
  res@xyLineThicknesses = (/30,30,30/)
  res@xyDashPatterns = (/0,0,0/)
;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
res_pn = True
;>>>>>>>>>>>>>>
plot_pn = new(20,graphic)
;dot_pn = new(ntime,graphic)
plot_pn(0) = gsn_csm_y (wks,score(:,:),res)
res2@xyLineColor = "orange"
res2@xyDashPattern = 14
plot_pn(1) = gsn_csm_y (wks,std,res2)

overlay(plot_pn(0),plot_pn(1))
draw (plot_pn(0))
frame(wks)

end if
;------------------------------------------------------------------------
end