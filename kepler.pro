function kepler, M, e, err

;This IDL function returns the solution x=eccentric anomaly to Kepler's elliptic equation
;f(x) = x - e*sin(x) - M = 0 via Halley's method to an accuracy of |f(x)|<err.
;Inputs are the mean anomaly M, eccentricity e, and the allowed error err.
;I suspect this code came from Danby's textbook, will look that up someday...
;This code does work in GNU :)
;
;Joe Hahn, jhahn@spacescience.org, 21 April 2013.

;Initial guess
twopi = 2*!dpi
M = M - twopi*fix(M/twopi)
s = M*0d + 1d
sinM = sin(M)
j = where(sinM lt 0d)
if (j[0] ne -1) then s[j] = -1d
x = M + s*(0.85d)*e
f = M*0d + 1d

max_count = 15
for p = 0l, n_elements(M) - 1 do begin
  count = 0
  while (( abs(f[p]) gt err ) and (count lt max_count)) do begin
    es = e[p]*sin(x[p])
    ec = e[p]*cos(x[p])
    f[p] = x[p] - es - M[p]
    df = 1d - ec
    ddf = es
    dddf = ec
    d1 = -f[p]/df
    d2 = -f[p]/(df + d1*ddf/2d)
    d3 = -f[p]/(df + d2*ddf/2d + d2*d2*dddf/6d)
    x[p]=x[p] + d3
    count = count + 1
  endwhile
endfor

if (count ge max_count) then print,'** kepler.pro failed to converge! **'

return,x
end

