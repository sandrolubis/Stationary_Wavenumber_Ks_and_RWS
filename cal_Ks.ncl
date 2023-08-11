; Compute zonal stationary waves Ks in the spherical coordinates
; See eq. (1)-(2) Henderson et al., (2017)
; by Sandro Lubis,  Sept 2014 
; GEOMAR
;===============================================================================
;===============================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;Some constants
gc=290       ;  Gas constant
ga=9.80665   ;  Gravitational acceleration
re=6378388   ;  Radius of the earth
sclhgt=8000. ;  Scale height



f = addfile("u200.daily.1979_2021.nc", "r")   ;U200


dum 	      = cd_calendar(f->time,0)
months 	    = toint(dum(:,1))
imonths     = ind(months.ge.6 .and. months.le.8)  ; chose season JJA
              
uvar  	= f->var131(imonths,0,:,:)  

umean	= dim_avg_n_Wrap(uvar,0)  ; Ubar

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


betatmp    = conform_dims(dimsizes(uvar),beta,1)
coslattmp  = conform_dims(dimsizes(uvar),coslat,1)
umtmp=conform_dims(dimsizes(uvar),umean,(/1,2/))
copy_VarCoords(uvar,umtmp)

delete(umean)

; U*(cos*lat)
ucoslat= umtmp*coslattmp
copy_VarCoords(uvar,ucoslat)

;d(U*(cos*lat)^2)/dy / coslat
dUdy = center_finite_diff_n(ucoslat, lat*pi/180., False,0,1) / coslattmp

;ddU/ddy
ddUdyy = center_finite_diff_n(dUdy, lat*pi/180., False,0,1)/re^2
;printMinMax(ddUdyy,True)

;betam
betam = betatmp - ddUdyy
copy_VarCoords(uvar,betam)

;beta/uvar
r=betam/umtmp
r@_FillValue=beta@_FillValue
r=where(r .le. 0, beta@_FillValue, r)   ; wave evanescence

;Ks
 ks=re*coslattmp*sqrt(r)   

copy_VarCoords(uvar,ks)
ks@long_name = "zonal stationary wavenumber"
ks@_FillValue = -1.e+21
printVarSummary(ks)
printMinMax(ks,True)

ks=where(ks.gt.50,ks@_FillValue,ks) ; retain only k=1-50 

;===============================================================================
; Save QG-PV Terms and n2 into one file:

;setfileoption("nc","Format","LargeFile")  ;NetCDF4Classic
;system( "rm " + "Ks."+year+".nc" )
;output = addfile("Ks."+year+".nc","c")
system( "rm " + "Ks.U200_1979_2021.nc" )
output = addfile("Ks.U200_1979_2021.nc","c")
output->ks    		= ks
output->betam     = betam
output->um        = umtmp



print("Done!!, see you again :)")
