;+
;
; Name:       fitsse
;
; Purpose:    function to fit HI elongation tracks using Davies et al. 2012
;             Apj Equations 8 and 9 with fixed transient half width lambda
;
; Parameters: X
;
; Keywords:  	-
;
;
; Called by fitall_sse.pro
;
;
; History:    March 2012
;
; Author:     Christian Moestl
;             SSL Berkeley
;-



FUNCTION fitsse_bauer, X

common myfit3,xueber3,yueber3, dst3, scnameueber3, lambda

degtorad=!dpi/180;

phi=X[0];  % constant angle phi measured FROM OBSERVER
v=X[1];  % constant velocity
t=xueber3-X[2]; (also deltat zur anfangszeit t0, die ja X[2] ist


;because of the acos the angle phi has to be positive
if scnameueber3 eq 'A' then phi=-phi;


;Davies et al. 2012 formula:
c=sin(lambda)
a=((dst3*(1+c))/(v*t))-cos(phi)
b=sin(phi)
fit=acos( (-b*c+a*sqrt((a^2+b^2-c^2)))/((a^2+b^2))      ) /degtorad


sizey=size(yueber3)

residue=double(0);

for i=0,sizey(1)-1 do begin
   if finite(fit(i)) eq 0 then fit(i)=0  ;if NaNs make the residue very high  (do this by  setting values to zero)
     residue=residue+(abs(yueber3(i))-abs(fit(i)) )^2;
;    residue=residue+( abs(yueber2(i)-fit(i)) )^2;
endfor

return, residue


END
