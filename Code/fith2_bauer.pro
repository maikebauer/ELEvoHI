; $Id:
;+
;
; Name:       fith
;
; Purpose:    function to fit HI elongation tracks using Lugaz formula
;
; Parameters: X
;
; Keywords:  	-
;
;
; Called by fithm.pro
;
;
; History:    May 2010/ update Feb 2011
;
; Author:     Christian Moestl
;             Space Research Institute, Austrian Academy of Sciences
;-



FUNCTION fith2_bauer, X

common myfit2,xueber2,yueber2, dst, scnameueber

degtorad=!dpi/180;

b=X[0];  % constant angle beta measured FROM OBSERVER
v=X[1];  % constant velocity
t=xueber2-X[2]; (also deltat zur anfangszeit t0, die ja X[2] ist


;because of the acos the angle beta has to be positive
if scnameueber eq 'A' then b=-b;


;my version
a=((2*dst)/(v*t))-cos(b)
c=sin(b)
fit=-acos( (-c^2+a*sqrt((c^2)*(-1+a^2+c^2)))/(c*(a^2+c^2))      ) /degtorad


sizey=size(yueber2)

residue=double(0);
for i=0,sizey(1)-1 do begin
   if finite(fit(i)) eq 0 then fit(i)=0  ;if NaNs make the residue very high  (do this by  setting values to zero)

     residue=residue+( abs(yueber2(i))-abs(fit(i)) )^2;
;    residue=residue+( abs(yueber2(i)-fit(i)) )^2;

endfor

return, residue


END
