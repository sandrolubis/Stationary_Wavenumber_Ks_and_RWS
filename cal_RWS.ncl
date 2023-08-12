; A script to Compute Rossby wave source from the daily mean input data.
; by Sandro Lubis, PhD Student 2015
; GEOMAR, Germany - Kiel

; Formulation:
; Compute components of rossby wave source: 
;
; S = -eta * div - (uchi * etax + vchi * etay)
;
; eta = absolute vorticity, 
; div = divergence,
; uchi & vchi = irrotational (divergent) wind components, 
; etax & etay = gradients of absolute vorticity.


; Method:
; Spherical Harmonics! 

; Notes:
; Input array must be on a global grid  
; Latitude must be in ascending order
; Fixed Grid not Gaussian!

;=====================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


do year=1979,2020

  f1 = addfile ("./data/u200."+year+".nc", "r")          ;U at 200
  f2 = addfile ("./data/v200."+year+".nc", "r")          ;V at 200
  
  
  dum       = cd_calendar(f1->time,0)
  months  = toint(dum(:,1))
  imonths = ind(months.ge.6 .and. months.le.8)  ; chose season JJA 


;=====================================================================


  U     = f1->var131    (:,0,::-1,:)   ; in ascending latitude order.  
  V     = f2->var132    (:,0,::-1,:)             

  printVarSummary(V)
  printVarSummary(U)

;=====================================================================
; Coriolis Parameter

lat    = tofloat(V&lat)

; pi
pi = atan(1.0)*4.

; Coriolis Force


f     =  2.*2.*pi/(60.*60.*24.)*sin(pi/180. * lat(:))  ; Coriolis parameter
f!0   = "lat"
f&lat =  lat
f@_FillValue = -1.e+21
 
ftmp      = conform_dims(dimsizes(V),f,1)   ; f in lat
copy_VarCoords(V,ftmp)
printVarSummary(ftmp)

;============== divergent (irrotational) wind components ==============


  div     = uv2dvF (U,V) 	      ; divergence
  uvd     = dv2uvF (div)        ; divergent  wind components
  uchi    = uvd(0,:,:,:)
  vchi    = uvd(1,:,:,:)

  copy_VarCoords(V,div)

  printVarSummary(div)
  printVarSummary(uchi)
  printVarSummary(vchi)
  

;============== Vorticity and Potential Vorticity =====================  

 vort     = uv2vrF(U,V)   ; relative vorticity

 eta = vort + ftmp        ; absolute vorticity
 printMinMax(eta,True)

 copy_VarCoords(V,eta)
 printVarSummary(eta)


etax=eta
etay=eta
gradsf (eta, etax, etay)   ; decompose eta into x,y comps.


;============== Rossby Wave Source =====================  

VS = -1*eta*div                     ; Vortex Stretching Term

AV = -(uchi*etax + vchi*etay)       ; Advection of Abs. Vorticity by Divergent Flow

RW = VS + AV

 copy_VarCoords(V,VS)
 copy_VarCoords(V,AV)
 copy_VarCoords(V,RW)

  VS@long_name  = "Vortex Stretching"
  AV@long_name  = "Advection of Abs. Vorticity"
  RW@long_name  = "Rossby Wave Source"
  div@long_name = "Absolute Vorticity"
  eta@long_name = "Wind Divergence"

 printVarSummary(VS)
 printVarSummary(AV)
 printVarSummary(RW)

 printMinMax(RW,True)


;==========smoothing (as you wish!)======== 

wrf_smooth_2d(VS,2)
wrf_smooth_2d(AV,2)
wrf_smooth_2d(RW,2)

;============== Save ===================== 

  setfileoption( "nc", "Format", "LargeFile" )             ; LargeFile or NetCDF4Classic
  system   ("/bin/rm -f output/RWS."+year+".nc")            ; remove any pre-existing file
  fout       = addfile("output/RWS."+year+".nc","c")        ; open output netCDF file
  fout->VS     = VS
  fout->AV     = AV
  fout->RWS    = RW
  ;fout->ETA    = eta
  ;fout->DIV    = div
  ;fout->time_bnds = f1->time_bnds



delete(dum)
delete(months)
delete(imonths)

delete(U)
delete(V)
delete(f)
delete(ftmp)
delete(div)
delete(uvd)
delete(uchi)
delete(vchi)
delete(vort)
delete(eta)
delete(etax)
delete(etay)
delete(VS)
delete(RW)
delete(AV)


  end do



