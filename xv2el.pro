pro xv2el, Mcen, x, y, z, vx, vy, vz, a, e, I, O, w, M, mass=mass

;Convert input Cartesian coordinates x,y,z and velocities vx,vy,vz into 
;Keplerian orbit elements a,e,I,O,w,M = semimajor axis, eccentricity,
;inclination, longitude of ascending node, argument of periapse, and mean anomaly.
;The central body's mass is Mcen. Here the gravitational constant G=1,
;so when Mcen=1 then units are such that masses are in solar units,
;distances are in units of astronomical units (AU), and velocities are in
;AU/(2*!pi*year), and the output elements are in units AU and radians.
;The GM of the system is Msun unless the mass keyword (which is the mass of the
;orbiting bodies) is used and then GM = Mcen + mass. Formulas used here probably
;came from Murray & Dermott's or maybe Danby's textbook on planetary dynamics.
;These codes do sucessfully execute in the Gnu equivalent to IDL, GDL.
;
;Joe Hahn
;jhahn@spacescience.org
;April 20, 2013


;Define TINY = small floating point number whose precision is determined by 
;whether the inputs are floats or doubles, so when certain quantities are < TINY
;then the orbit is assumed to have e=0 or I=0 as appropriate.
sz = size(x)
Ndims = sz[0]
arr_code = sz[Ndims+1]
if (arr_code eq 4) then TINY = 2e-7
if (arr_code eq 5) then TINY = 2d-15

;inclination I
hx = y*vz - z*vy
hy = z*vx - x*vz
hz = x*vy - y*vx
hxy = sqrt(hx^2 + hy^2)
h2 = hx^2 + hy^2 + hz^2
h = sqrt(h2)
I = x*0
r = sqrt(x^2 + y^2 + z^2)
j = where(r gt TINY)
if (j[0] ne -1) then I[j] = acos(hz[j]/h[j])

;set mu = GM = Mcen unless the mass keyword is set
mu = Mcen
if (keyword_set(mass)) then mu = Mcen + mass

;eccentricity e
zero = x*0
ex = zero
if (j[0] ne -1) then ex[j] = ( vy[j]*hz[j] - vz[j]*hy[j] )/mu[j] - x[j]/r[j]
ey = zero
if (j[0] ne -1) then ey[j] = ( vz[j]*hx[j] - vx[j]*hz[j] )/mu[j] - y[j]/r[j]
ez = zero
if (j[0] ne -1) then ez[j]=( vx[j]*hy[j] - vy[j]*hx[j] )/mu[j] - z[j]/r[j]
e = sqrt(ex^2 + ey^2 + ez^2)

;longitude of ascending node O
O = x
if (j[0] ne -1) then O[j] = atan(hx[j], -hy[j])

;argument of periapse w
ec = ex*cos(O) + ey*sin(O)
es = zero
j = where(I gt TINY)
if (j[0] ne -1) then es[j] = ez[j]/sin(I[j])
w = atan(es, ec)

;semimajor axis a
v2 = vx^2 + vy^2 + vz^2
a = zero
j = where(r gt TINY)
if (j[0] ne -1) then a[j]=1d/( 2d/r[j]-v2[j]/mu[j] )

;mean anomaly M
M = zero

;M for elliptic orbits
rv = x*vx + y*vy + z*vz
k = where(rv lt 0d)
ec = zero
ec[j] = 1d - r[j]/a[j]
j = where(a gt 0d)
if (j[0] ne -1) then begin 
  es[j] = rv[j]/sqrt( abs(a[j])*mu[j] ) 
  Ea = zero 
  Ea[j] = atan(es[j], ec[j]) 
  M[j] = Ea[j] - es[j] 
endif

;M for hyperbolic orbits
j = where(a lt 0d)
if (j(0) ne -1) then begin
  ZZ = zero
  ZZ[j] = ec[j]/( e[j]>1d )
  F = zero
  F[j] = alog( abs(ZZ[j] + sqrt(abs(ZZ[j]^2 - 1d))) )
  if (k[0] ne -1) then F[j] = -F[j]
  M[j] = e[j]*sinh(F[j]) - F[j]
endif

;M for parabolic orbits
j = where((a eq 0d) and (r gt TINY))
if (j[0] ne -1) then begin
  q = zero
  q(j) = (0.5d)*h2[j]/mu[j]
  tau = zero
  tau[j] = sqrt( r[j]/q[j] - 1d )
  if (j[0] ne -1) then tau[j] = -tau[j]
  M[j] = sqrt(2d)*( tau[j]+(tau[j]^3)/3d )
endif

;arrange angles to be between 0 and 2*Pi
twopi = 2*!dpi
j = where(O lt 0d)
if (j[0] ne -1) then O[j] = O[j] + tp
j = where(w lt 0d)
if (j[0] ne -1) then w[j] = w[j] + tp
j = where(M lt 0d)
if (j[0] ne -1) then M[j] = M[j] + tp

return
end
