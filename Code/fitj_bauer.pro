; $Id:
;+
;
; Name:       fitj
;
; Purpose:    function to fit HI elongation tracks using Sheeleys formula
;
; Parameters: X
;
; Keywords:  	-
;
;
; Called by fitelongation.pro
;
;
; History:    June 20009
;
; Author:     Christian Möstl
;             Space Research Institute, Austrian Academy of Sciences
;-




FUNCTION fitj_bauer, X

common myfit,xueber,yueber, H0ueber


;after Sheeley 2008 ApJ
;X=[beta, vr, t0]; radians km/s sec
degtorad=!dpi/180;

delta=90*degtorad-X[0]; delta wird v. Plane of Sky weg gemessen
rho=X[1]*(xueber-X[2])/H0ueber;

fit=atan(rho*cos(delta)/(1-rho*sin(delta) ) )/degtorad

;window,4
;plot, xueber,yueber
;window, 5
;plot, xueber,fit


;stop

sizey=size(yueber)
residue=0;
for i=0,sizey(1)-1 do begin
  residue=residue+( abs(yueber(i))-abs(fit(i)) )^2;
endfor
; print, residue
return, residue

;minimierung zwischen datenpunkten und theoretischer Kurve...


END
