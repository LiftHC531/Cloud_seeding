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
;-------------------------------------------------------------------------
setfileoption("bin","ReadByteOrder","LittleEndian")
UNDEF = 9.96921e+36
files_a = systemfunc("ls -d ./*_score.bin")
member = dimsizes(files_a)
print(files_a)
dam = (/"Shimen","Zengwen","Nunhua"/)
nv = 1+9+1;3
ndam = 3
nt = 24
var = new ( (/nv,ndam,nt/), float, UNDEF)
;----
id = (/5,10,2,3,7,8,6/)
nid = dimsizes(id)
var_rg = new ( (/nv,ndam,2/), float, UNDEF) ;min & max
do mm = 0,member-1
   ;std = var
   ;std = 3.0
   date = stringtochar(files_a(mm))
   ;print(date(2:14))
   if (date(13:14).eq."12") then
       ti = 12
       tf = 23
       day = sprinti("%2.2i", stringtointeger(date(10:11)) + 1)
   else
       ti = 0
       tf = 12
       day = sprinti("%2.2i", stringtointeger(date(10:11)) )
   end if
   nrh = 0 
   do vv = 0,nv-1
      var(vv,:,:) = fbindirread(files_a(mm), nrh, (/ndam,nt/), "float")
      nrh = nrh + 1
      if (vv.eq.2.or.vv.eq.10.or.vv.eq.12) then
         var(vv,:,:) = where(var(vv,:,:).ne.0.0,var(vv,:,:),UNDEF)
      end if
      do dd = 0,ndam-1
         var_rg(vv,dd,0) = min(var(vv,dd,ti:tf))
         var_rg(vv,dd,1) = max(var(vv,dd,ti:tf))
      end do
   end do
   var_rg(2,:,:) = var_rg(2,:,:)*1.E-3 ;base of cloud
   var_rg(10,:,:) = var_rg(10,:,:)*1.E-3 ; thickness of cloud
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ;sprintf("%5.2f", x)
   fname = date(2:9)+day+"_sheet.prn"
   write_table(fname, "w", [/date(2:9)+day/], "%s ")
   do dd = 0,ndam-1
      write_table(fname, "a", [/dam(dd)+"---------"/], "%s ")
      do vv = 0,nid-1
         write_table(fname, "a", [/var_rg(id(vv),dd,0),"~",var_rg(id(vv),dd,1)/], "%5.2f%s%5.2f")
      end do
   end do

end do
;------------------------------------------------------------------------
end

