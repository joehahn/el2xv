;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function kepler,M,e,err

;returns the solution x=eccentric anomaly to kepler's elliptic equation
; f(x)=x-e*sin(x)-M=0 via Halley's method
;to an accuracy of |f(x)|<err.


;initial guess
tp=6.283185307179586d0
M=M-tp*fix(M/tp)
s=M*0.d0+1.d0
sinM=sin(M)
j=where(sinM lt 0.0d)
if (j(0) ne -1) then s(j)=-1.d0
x=M+s*(0.85d0)*e
f=M*0.0d0+1.d0

max_count=15
for p=0l,n_elements(M)-1 do begin
  count=0
  while ((abs(f(p)) gt err) and (count lt max_count)) do begin
    es=e(p)*sin(x(p))
    ec=e(p)*cos(x(p))
    f(p)=x(p)-es-M(p)
    df=1.0d0-ec
    ddf=es
    dddf=ec
    d1=-f(p)/df
    d2=-f(p)/(df+d1*ddf/2.0d0)
    d3=-f(p)/(df+d2*ddf/2.0d0+d2*d2*dddf/6.0d0)
    x(p)=x(p)+d3
    count=count+1
  endwhile
  ;print, count
endfor

if (count ge max_count) then $
  print,'** kepler.pro failed to converge! **'

return,x
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


