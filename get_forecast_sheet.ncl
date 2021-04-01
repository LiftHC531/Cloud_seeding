load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;load "/home/mlhchen/big_WRFUserARW.ncl"

;nrh   variable
;0      score
;1      LWP
;2      LCL
;3      w 1000m
;4      dBZ
;5      ctt
;6      V 1000m
;7      RH 1500 m
;8      RH 5000 m
;9      Rainfall rate
;10     thickness of cloud
;11     top of cloud
;12     base of cloud    
begin
err = NhlGetErrorObjectId()
  setvalues err
    "errLevel" : "Fatal"          ; only report Fatal errors
  end setvalues
;-------------------------------------------------------------------------
setfileoption("bin","ReadByteOrder","LittleEndian")
UNDEF = 9.96921e+36
files_a = systemfunc("ls -d ./*_score.bin")
member = dimsizes(files_a)
print(files_a)
dam = (/"Shimen","Baoshan","Zengwen","Nanhua","Feitsui"/)
nv = 1+9+1;3
ndam = dimsizes(dam) 
nt = 24
var = new ( (/nv,ndam,nt/), float, UNDEF)
;----
id = (/5,10,2,3,7,8,6/)
nid = dimsizes(id)
var_rg = new ( (/member,nv,ndam,3/), float, UNDEF) ;min & max
do mm = 0,member-1
   ;std = var
   ;std = 3.0
   date = stringtochar(files_a(mm))
   year = stringtointeger(date(2:5))
   month_int = stringtointeger(date(7:8))
   day_int = stringtointeger(date(10:11))
   ;print(date(2:14))
   if (date(13:14).eq."12") then
       ti = 12
       tf = 23
       if (month_int .eq. 1 .or. month_int .eq. 3 .or. \
           month_int .eq. 5 .or. month_int .eq. 7 .or. \
           month_int .eq. 8 .or. month_int .eq. 10 .or. \
           month_int .eq. 12) then
           if (day_int+1 .gt. 31) then
              day = "01"
              month = sprinti("%2.2i", month_int + 1)
           else
              day = sprinti("%2.2i", day_int+1)
              month = sprinti("%2.2i", month_int)
           end if
       else if (month_int .eq. 2) then
           if (day_int+1 .gt. 28) then
              day = "01"
              month = sprinti("%2.2i", month_int + 1)
           else
              day = sprinti("%2.2i", day_int+1)
              month = sprinti("%2.2i", month_int)
           end if
       else
           if (day_int+1 .gt. 30) then
              day = "01"
              month = sprinti("%2.2i", month_int + 1)
           else
              day = sprinti("%2.2i", day_int+1)
              month = sprinti("%2.2i", month_int)
           end if
       end if
       end if
   else
       ti = 0
       tf = 12
       day = sprinti("%2.2i", day_int)
       month = sprinti("%2.2i", month_int)
   end if
   nrh = 0 
   do vv = 0,nv-1
      var(vv,:,:) = fbindirread(files_a(mm), nrh, (/ndam,nt/), "float")
      nrh = nrh + 1
      if (vv.eq.2.or.vv.eq.10.or.vv.eq.12) then
         var(vv,:,:) = where(var(vv,:,:).ne.0.0,var(vv,:,:),UNDEF)
      end if
      do dd = 0,ndam-1
         var_rg(mm,vv,dd,0) = min(var(vv,dd,ti:tf))
         var_rg(mm,vv,dd,1) = max(var(vv,dd,ti:tf))
         var_rg(mm,vv,dd,2) = avg(var(vv,dd,ti:tf))
      end do
   end do
   var_rg(mm,2,:,:) = var_rg(mm,2,:,:)*1.E-3 ;base of cloud
   var_rg(mm,10,:,:) = var_rg(mm,10,:,:)*1.E-3 ; thickness of cloud
   do vv = 0,nv-1
     do dd = 0,ndam-1
      do ddd = 0,3-1
         if (ismissing(var_rg(mm,vv,dd,0))) then
             var_rg(mm,vv,dd,0) = -999
         end if
         if (ismissing(var_rg(mm,vv,dd,1))) then
             var_rg(mm,vv,dd,1) = -999
         end if
         if (ismissing(var_rg(mm,vv,dd,2))) then
             var_rg(mm,vv,dd,2) = -999
         end if
      end do
     end do
   end do
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ;sprintf("%4.1f", x)
   fname = year+"-"+month+"-"+day+"_sheet.prn"
   write_table(fname, "w", [/date(2:9)+day/], "%s ")
   do dd = 0,ndam-1
      write_table(fname, "a", [/dam(dd)+"---------"/], "%s ")
      do vv = 0,2
         write_table(fname, "a", [/var_rg(mm,id(vv),dd,2)/], "%4.1f")
      end do
      do vv = 2,nid-1
         write_table(fname, "a", [/var_rg(mm,id(vv),dd,0),"~",var_rg(mm,id(vv),dd,1)/], "%4.1f%s%4.1f")
      end do
   end do

end do

files_b = "12hr_rain.bin"
nday = 5
rain = new ( (/nday,ndam/), float, UNDEF)
rh = new ( (/nday,ndam/), float, UNDEF)
nrh = 0
do dd = 0,ndam-1
  rain(:,dd) = fbindirread(files_b, nrh, (/nday/), "float")
  rh(:,dd) = fbindirread(files_b, nrh+1, (/nday/), "float")
  nrh = nrh + 2 
end do
do dd = 0,ndam-1
   do da = 0,nday-1
      if (ismissing(rain(da,dd))) then
          rain(da,dd) = 0.0
      end if
      if (ismissing(rh(da,dd))) then
          rh(da,dd) = 0.0
      end if
   end do
end do 

   fname = "forecast_sheet.prn"
      write_table(fname, "w", [/"---------"/], "%s ")
   do dd = 0,ndam-1
      write_table(fname, "a", [/dam(dd)+"---------"/], "%s ")
      write_table(fname, "a", [/rh(0,dd),",",rh(1,dd),",",rh(2,dd),",",rh(3,dd),",",rh(4,dd)/], \
                                 "%4.1f%s%4.1f%s%4.1f%s%4.1f%s%4.1f%s%4.1f")
      write_table(fname, "a", [/rain(0,dd),",",rain(1,dd),",",rain(2,dd),",",rain(3,dd),",",rain(4,dd)/], \
                                 "%4.1f%s%4.1f%s%4.1f%s%4.1f%s%4.1f%s%4.1f")
      do vv = 0,2
         write_table(fname, "a", [/var_rg(0,id(vv),dd,2),",",var_rg(1,id(vv),dd,2),",", \
                                   var_rg(2,id(vv),dd,2),",",var_rg(3,id(vv),dd,2),",", \
                                   var_rg(4,id(vv),dd,2)/], "%4.1f%s%4.1f%s%4.1f%s%4.1f%s%4.1f%s%4.1f")
      end do
      do vv = 3,nid-1
         if (vv.ne.4.and.vv.ne.5) then
         write_table(fname, "a", [/var_rg(0,id(vv),dd,0),"~",var_rg(0,id(vv),dd,1),",", \
                                   var_rg(1,id(vv),dd,0),"~",var_rg(1,id(vv),dd,1),",", \
                                   var_rg(2,id(vv),dd,0),"~",var_rg(2,id(vv),dd,1),",", \
                                   var_rg(3,id(vv),dd,0),"~",var_rg(3,id(vv),dd,1),",", \
                                   var_rg(4,id(vv),dd,0),"~",var_rg(4,id(vv),dd,1) \
                             /], "%4.1f%s%4.1f%s%4.1f%s%4.1f%s%4.1f%s%4.1f%s%4.1f%s%4.1f%s%4.1f%s%4.1f")
         else
         write_table(fname, "a", [/var_rg(0,id(vv),dd,0),"~",var_rg(0,id(vv),dd,1),",", \
                                   var_rg(1,id(vv),dd,0),"~",var_rg(1,id(vv),dd,1),",", \
                                   var_rg(2,id(vv),dd,0),"~",var_rg(2,id(vv),dd,1),",", \
                                   var_rg(3,id(vv),dd,0),"~",var_rg(3,id(vv),dd,1),",", \
                                   var_rg(4,id(vv),dd,0),"~",var_rg(4,id(vv),dd,1) \
                             /], "%3.0f%s%3.0f%s%3.0f%s%3.0f%s%3.0f%s%3.0f%s%3.0f%s%3.0f%s%3.0f%s%3.0f")
         
         end if
      end do
   end do


end

