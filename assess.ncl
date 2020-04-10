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
;----------------------------------
setfileoption("bin","WriteByteOrder","LittleEndian")
UNDEF = 9.96921e+36
fName = "lwplclwdbz.bin"
fName2 = "cloud_char.bin"
nv = 9 
 var1 = new((/nv,ntimes,nlat,nlon/),float,UNDEF)
 var2 = new((/3,ntimes,nlat,nlon/),float,UNDEF)
do vv = 0,nv-1
   var1(vv,:,:,:) = fbindirread(fName, vv, (/ntimes,nlat,nlon/), "float")
end do
nrh = 0
do t = 0,ntimes-1
  var2(0,t,:,:) = fbindirread(fName2, nrh, (/nlat,nlon/), "float") ;thick
  var2(1,t,:,:) = fbindirread(fName2, nrh+1, (/nlat,nlon/), "float") ;top
  var2(2,t,:,:) = fbindirread(fName2, nrh+2, (/nlat,nlon/), "float") ;base
  nrh = nrh + 3
end do
var1(4-1,:,:,:) = 10.0^var1(4-1,:,:,:) 
lat_c = (/136,80,74/)  ; y
lon_c = (/ 93,71,71/)  ; x
 var1_avg = new((/3,nv,ntimes/),float,UNDEF)
 var2_avg = new((/3,3,ntimes/),float,UNDEF)
do dam = 0,3-1
 do t = 0,ntimes-1
    do vv  = 0,nv-1
       stat = stat_dispersion(var1(vv,t,lat_c(dam)-5:lat_c(dam)+5,lon_c(dam)-5:lon_c(dam)+5),False)
       ;printVarSummary(stat)
       if (vv.eq.4) then
          var1_avg(dam,vv,t) = stat(2) ;min.
       else
          var1_avg(dam,vv,t) = stat(8) ;median
       end if
       delete(stat)
    end do
    do vv = 0,3-1
       stat = stat_dispersion(var2(vv,t,lat_c(dam)-5:lat_c(dam)+5,lon_c(dam)-5:lon_c(dam)+5),False)
       ;printVarSummary(stat)
       var2_avg(dam,vv,t) = stat(8) ;median
       delete(stat)
    end do
 end do
end do
var1_avg(:,3,:) = log10(var1_avg(:,3,:))
print(var1_avg(0,0,:))
;----------------------------------
 score = new((/3,ntimes/),float,UNDEF)
 score = 0.0
do dam = 0,3-1
 do t = 0,ntimes-1
    do vv  = 0,nv-1
       if (vv.eq.0) then ;LWP
          if(var1_avg(dam,vv,t).ge.50.0.and.var1_avg(dam,vv,t).lt.100.0) then
             score(dam,t) = score(dam,t) + 1.0
          else if(var1_avg(dam,vv,t).ge.100.0.and.var1_avg(dam,vv,t).lt.200.0) then
             score(dam,t) = score(dam,t) + 3.0
          else if(var1_avg(dam,vv,t).ge.200.0) then
             score(dam,t) = score(dam,t) + 5.0
          end if
          end if
          end if
       else if (vv.eq.1) then ;LCL
          if(var1_avg(dam,vv,t).le.1500.0.and.var1_avg(dam,vv,t).gt.1000.0) then
             score(dam,t) = score(dam,t) + 1.0
          else if(var1_avg(dam,vv,t).le.1000.0.and.var1_avg(dam,vv,t).gt.500.0) then
             score(dam,t) = score(dam,t) + 3.0
          else if(var1_avg(dam,vv,t).le.500.0) then
             score(dam,t) = score(dam,t) + 5.0
          end if
          end if
          end if
       else if (vv.eq.2) then ;w 1000 m
            var1_avg(dam,vv,t)=var1_avg(dam,vv,t)*100.0
          if(var1_avg(dam,vv,t).ge.0.0.and.var1_avg(dam,vv,t).lt.1.0) then
             score(dam,t) = score(dam,t) + 1.0
          else if(var1_avg(dam,vv,t).ge.1.0.and.var1_avg(dam,vv,t).lt.3.0) then
             score(dam,t) = score(dam,t) + 3.0
          else if(var1_avg(dam,vv,t).ge.3.0) then
             score(dam,t) = score(dam,t) + 5.0
          end if
          end if
          end if
       else if (vv.eq.3) then ;dBZ
          if(var1_avg(dam,vv,t).ge.5.0.and.var1_avg(dam,vv,t).lt.10.0) then
             score(dam,t) = score(dam,t) + 1.0
          else if(var1_avg(dam,vv,t).ge.10.0.and.var1_avg(dam,vv,t).lt.25.0) then
             score(dam,t) = score(dam,t) + 3.0
          else if(var1_avg(dam,vv,t).ge.25.0) then
             score(dam,t) = score(dam,t) + 5.0
          end if
          end if
          end if
       else if (vv.eq.4) then ;ctt
          if(var1_avg(dam,vv,t).le.0.0.and.var1_avg(dam,vv,t).gt.-5.0) then
             score(dam,t) = score(dam,t) + 1.0
          else if(var1_avg(dam,vv,t).le.-5.0.and.var1_avg(dam,vv,t).gt.-15.0) then
             score(dam,t) = score(dam,t) + 3.0
          else if(var1_avg(dam,vv,t).le.-15.0) then
             score(dam,t) = score(dam,t) + 5.0
          end if
          end if
          end if
       else if (vv.eq.5) then ;wind speed 1000 m
          if(var1_avg(dam,vv,t).le.15.0.and.var1_avg(dam,vv,t).gt.10.0) then
             score(dam,t) = score(dam,t) + 1.0
          else if(var1_avg(dam,vv,t).le.10.0.and.var1_avg(dam,vv,t).gt.5.0) then
             score(dam,t) = score(dam,t) + 3.0
          else if(var1_avg(dam,vv,t).le.5.0) then
             score(dam,t) = score(dam,t) + 5.0
          end if
          end if
          end if
       else if (vv.eq.6) then ;RH 1500 m
          if(var1_avg(dam,vv,t).ge.50.0.and.var1_avg(dam,vv,t).lt.60.0) then
             score(dam,t) = score(dam,t) + 1.0
          else if(var1_avg(dam,vv,t).ge.60.0.and.var1_avg(dam,vv,t).lt.70.0) then
             score(dam,t) = score(dam,t) + 3.0
          else if(var1_avg(dam,vv,t).ge.70.0) then
             score(dam,t) = score(dam,t) + 5.0
          end if
          end if
          end if
       else if (vv.eq.7) then ;RH 5000 m
          if(var1_avg(dam,vv,t).ge.50.0.and.var1_avg(dam,vv,t).lt.60.0) then
             score(dam,t) = score(dam,t) + 1.0
          else if(var1_avg(dam,vv,t).ge.60.0.and.var1_avg(dam,vv,t).lt.70.0) then
             score(dam,t) = score(dam,t) + 3.0
          else if(var1_avg(dam,vv,t).ge.70.0) then
             score(dam,t) = score(dam,t) + 5.0
          end if
          end if
          end if
       else if (vv.eq.8) then ;rainfall rate
          if(var1_avg(dam,vv,t).ge.0.1.and.var1_avg(dam,vv,t).lt.0.3) then
             score(dam,t) = score(dam,t) + 1.0
          else if(var1_avg(dam,vv,t).ge.0.3.and.var1_avg(dam,vv,t).lt.0.5) then
             score(dam,t) = score(dam,t) + 3.0
          else if(var1_avg(dam,vv,t).ge.0.5) then
             score(dam,t) = score(dam,t) + 5.0
          end if
          end if
          end if
          

       end if
       end if
       end if
       end if
       end if
       end if
       end if
       end if
       end if
    end do
    if(var2_avg(dam,0,t).gt.0.0.and.var2_avg(dam,0,t).lt.1000.0) then
       score(dam,t) = score(dam,t) + 1.0
    else if(var2_avg(dam,0,t).ge.1000.0.and.var2_avg(dam,0,t).lt.2000.0) then
       score(dam,t) = score(dam,t) + 3.0
    else if(var2_avg(dam,0,t).ge.2000.0) then
       score(dam,t) = score(dam,t) + 5.0
    end if
    end if
    end if

    score(dam,t) = score(dam,t)/10.0
    print("scores = "+score(0,t))

 end do
end do ; dam
 do t = 2,3;ntimes-1
    print(t+","+var1_avg(0,:,t))
    print("---")
    print(t+","+var2_avg(0,:,t))
    print("---------------")
 end do
 File_out = "score.bin"
 system("rm "+File_out)
 setfileoption("bin","WriteByteOrder","LittleEndian")
 fbindirwrite(File_out,score(:,:))
 do vv = 0,nv-1
    fbindirwrite(File_out,var1_avg(:,vv,:))
 end do
    fbindirwrite(File_out,var2_avg(:,0,:))
;------------------------------------------------------------------------
end
