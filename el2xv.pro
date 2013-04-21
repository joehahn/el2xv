pro el2xv, Mcen, a, e, I, O, w, M, x, y, z, vx, vy, vz, mass=mass, error=error

;Convert orbit elements a,e,I,O,w,M to cartesian coordinates and velocities
;x,y,z,vx,vy,vx. See xv2el.pro for a description of the inputs. This routine
;also calls kepler.pro to solve Kepler's equation, and the optional error
;parameter is the tolerance in that solution and defaults to error=1d-8.
;Formulas used here probably come from textbook by Murray & Dermott or Danby.


;Set default error tolerance
if (keyword_set(error) ne 1) then error = 1d-8

;Calculate the rotation matrices
so = sin(O)
co = cos(O)
sp = sin(w)
cp = cos(w)
si = sin(I)
ci = cos(I)
d11 =  cp*co - sp*so*ci
d12 =  cp*so + sp*co*ci
d13 =  sp*si
d21 = -sp*co - cp*so*ci
d22 = -sp*so + cp*co*ci
d23 =  cp*si

;Storage
zero = 0*a
xx = zero
yy = zero
vxx = zero
vyy = zero

;Choose appropriate GM
GM = Mcen
if (keyword_set(mass) eq 1) then GM = Mcen + mass

;Loop over all elements in input array
for j = 0l, n_elements(a) - 1 do begin 

  ;for elliptic elements
  if (a[j] ge 0d) then begin 

    ;solve kepler's eqn.
    EA = kepler( M[j], e[j], error)
    cE = cos(EA)
    sE = sin(EA)

    ;coordinates in orbit frame
    e_root = sqrt( 1d - e[j]^2 )
    xx[j] = a[j]*( cE - e[j] )
    yy[j] = a[j]*e_root*sE 
    r = a[j]*( 1d - e[j]*cE )
    na2 = sqrt( GM*a[j] )
    vxx[j] = -na2*sE/r 
    vyy[j] = na2*e_root*cE/r 

  endif

  ;for hyperbolic elements
  if (a[j] lt 0d) then begin 

    ;solve kepler's eqn.
    F = hyper_kepler( M[j], e[j], error)
    chF = cosh(F)
    shF = sinh(F)

    ;coordinates in orbit frame
    sqe = sqrt( e[j]^2 - 1d ) 
    sqgma = sqrt( -GM*a[j] )
    ri = -1d/( a[j]*( e[j]*chF - 1d ) )
    xx[j] = -a[j]*( e[j] - chF )
    yy[j] = -a[j]*sqe*shF 
    vxx[j] = -ri*sqgma*shF 
    vyy[j] = ri*sqgma*sqe*chF

  endif

endfor

;rotate to reference frame
x  = d11*xx  + d21*yy
y  = d12*xx  + d22*yy
z  = d13*xx  + d23*yy
vx = d11*vxx + d21*vyy
vy = d12*vxx + d22*vyy
vz = d13*vxx + d23*vyy

end

