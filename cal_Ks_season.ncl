; Compute zonal stationary waves Ks in the spherical coordinates
; See eq. (1)-(2) Henderson et al., (2017)
; Choose your season of interest! (input monthly)
; by Sandro Lubis,  Sept 2014 
; GEOMAR
;===============================================================================
;===============================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;Some constants
gc=290         ;  Gas constant
ga=9.80665   ;  Gravitational acceleration
re=6378388   ;  Radius of the earth
sclhgt=8000. ;  Scale height



f = addfile("u200.mon.1979_2021.nc", "r")   ;U200


uvar  	= f->var131(:,0,:,:)  

umean	= month_to_season(uvar, (/"JJA"/))  ; choose seasons
umean@_FillValue = -1.e+21
umean=where(umean.lt.0,umean@_FillValue,umean)

printVarSummary(uvar)
printVarSummary(umean)


lat     = dble2flt( f->lat)
lon    = dble2flt(f->lon)
nlat   = dimsizes(lat)
nlon  = dimsizes(lon)

; pi
pi = atan(1.0)*4.

; cosine
coslat = cos(lat(:)*pi/180.)


; beta (first term)
beta =  2.*2.*pi/(60.*60.*24.)*coslat  / re
beta!0 = "lat"
beta&lat = lat
beta@_FillValue = -1.e+21


betatmp    = conform_dims(dimsizes(umean),beta,1)
coslattmp  = conform_dims(dimsizes(umean),coslat,1)




; U*(cos*lat)
ucoslat= umean*coslattmp
;printMinMax(ucoslat,True)


;d(U*(cos*lat)^2)/dy / coslat
dUdy = center_finite_diff_n(ucoslat, lat*pi/180., False,0,1) / coslattmp
;printMinMax(dUdy,True)

;ddU/ddy
ddUdyy = center_finite_diff_n(dUdy, lat*pi/180., False,0,1)/re^2
;printMinMax(ddUdyy,True)

;betam
betam = betatmp - ddUdyy
copy_VarCoords(umean,betam)

;beta/uvar
r=betam/umean
r@_FillValue=beta@_FillValue
r=where(r .le. 0, beta@_FillValue, r)   ; wave evanescence

;Ks
 ks=re*coslattmp*sqrt(r)   

copy_VarCoords(umean,ks)
ks@long_name = "zonal stationary wavenumber"
ks@_FillValue = -1.e+21
printVarSummary(ks)
printMinMax(ks,True)

ks=where(ks.gt.50,ks@_FillValue,ks) ; retain only k=1-50 
printVarSummary(ks)



;===============================================================================
; Save Ks

setfileoption("nc","Format","LargeFile")  ;NetCDF4Classic
system("rm " + "Ks.season.nc")
output = addfile("Ks.season.nc","c")
output->ks    		= ks
output->betam     = betam
output->um        = umean



print("Done!!, see you again :)")
