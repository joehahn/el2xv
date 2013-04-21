;test_el2xv.pro

;The following IDL script tests the el2xv.pro and xv2el.pro procedures, and it also
;illustrates their usage. First convert Keplerian orbit elements into cartesian
;coordinates and velocities using the el2xv procedure, then convert those
;coordinates back into orbit elements, and verify that the original orbit
;is recovered to within numerical errors. Units are AU and AU/(2*Pi*year).
;To execute, start IDL and enter @test_xv2el . These codes do sucessfully
;execute in the Gnu equivalent to IDL which is GDL.
;
;Joe Hahn
;jhahn@spacescience.org
;April 20, 2013


;Choose an orbit.
a = 1.5d
e = 0.3d
I = 0.1d
O = 1d
w = 1.5d
M = 2d

;Choose masses of central body and the orbiting body.
Mcen = 1d
mass = 1d-6

;Convert orbit elements into cartesian coordinates and velocities, then
;back to orbit elements, assuming GM = 1.
el2xv, Mcen, a, e, i, O, w, M, x, y, z, vx, vy, vz
xv2el, Mcen, x, y, z, vx, vy, vz, a2, e2, I2, O2, w2, M2

;Confirm that a=a2, e=e2 etc within numerical precision
print, a, a2, a - a2
print, e, e2, e - e2
print, I, I2, I - I2
print, O, O2, O - O2
print, w, w2, w - w2
print, M, M2, M - M2
print, ''

;Repeat the above assuming GM = 1 + mass
el2xv, Mcen, a, e, i, O, w, M, x, y, z, vx, vy, vz, mass=mass
xv2el, Mcen, x, y, z, vx, vy, vz, a2, e2, I2, O2, w2, M2, mass=mass
print, a, a2, a - a2
print, e, e2, e - e2
print, I, I2, I - I2
print, O, O2, O - O2
print, w, w2, w - w2
print, M, M2, M - M2
print, ''

;Now choose a circular orbit with a=1 and convert orbit elements to coordinates with GM=1.
;When GM=1, time increments by 2*Pi every orbit with a=1, so the 
;orbit period T=2*Pi*a/v should be 2*Pi and thus speed v should be =1.
a = 1d
e = 0d
el2xv, Mcen, a, e, i, O, w, M, x, y, z, vx, vy, vz
v = sqrt(vx^2 + vy^2 + vz^2)
print, vx, vy, vz, v

