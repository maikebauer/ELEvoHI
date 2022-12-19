;+
;
; Name:       fitall_sse_out
;
; Purpose:    fit elongation vs. time to infer constant direction, speed, launch and arrival times (for ICMEs, CIRs)
;           using fitting methods which assume different geometries of the transient front
;             (e.g. Sheeley et al. 1999 ApJ, Rouillard et al. 2008 GRL, Moestl et al. 2009,2010,2011; Davies et al. 2012 ApJ)
;         to be used within the SolarSoft SATPLOT package
;             The arrival times at a spacecraft position are corrected for apex/flank hits (Moestl et al. 2011; Moestl and Davies, 2013)
;
; Parameters:  file_in... the name of the satplot .ht file which contains the elongation-time data of the transient
;
; Keywords:   all optional:
;               cut...  two integers indicating how many data points of the CME
;                   are not used at the beginning and end of the tracking data;
;                 default=[0,0]
;               dir...  direction measured from the observer in degree,
;                       default=+60 degree for STEREO-B and -60 degree for STEREO-A
;               lambda..half width of the transient within the SSE model in degree (Davies et al. 2012 ApJ)
;                       default=45 degree
;
; Calling sequence: results=fitall_sse_out(filename, cut=[cutfront, cutback], dir=direction, lambda=setlambda)
;
; Example: results=fitall_sse_out('satplot__20110213_120000__20110218_120000__pa270_d10_B.ht', dir=70, lambda=50)
;
; Side effects: uses spice to get Earth and STEREO positions
;               calls fitj.pro, fith.pro, fitsse.pro for minimization with different models of the transient front
;               calls fit_geometry.pro for visualizing spacecraft positions and transient directions
;               produces an output file for results (.txt) and figures of the fits and resulting CME geometry (.eps, .jpg)
;
; History:    12 Feb 2012: started rewriting to be used with satplot
;             7 Mar 2012: continued, now works with .ht files from satplot
;             12 Mar 2012: streamlined the input process
;             29 Mar 2012: included SSE fitting
;             31 Mar 2012: included geometry visualization
;             13 May 2013: included output of results into txt file
;             29 January 2014: fixed rounding errors in terminal and graphical outputs
;
; Authors:    Christian Moestl, SSL Berkeley, USA and University of Graz, Austria
;             Tanja Rollett - University of Graz, Austria
;
;-

function fitall_sse_tro, cut=cut, start=start, bflag=bflag, dir=dir, lambda=setlambda, single=single, tr_num=tr_num,  _extra = extra_keywords, sc=sc

;---------------------------------------
;read in .ht files produced by satplot
;str = STRARR(200)
;
;OPENR, 10, file_in
;   dummy = ''
;   i = 0L
;   WHILE NOT(EOF(10)) DO BEGIN
;      READF, 10, dummy
;      str[i] = dummy
;      i = i + 1
;   ENDWHILE
;CLOSE, 10
;str = str[0:i-1]
;searchsc=strmid(str,13,1);
;spacecraft= searchsc(3);
;a = where( strmid(str,0,1) eq '#' )
;a = max(a)+1
;str = str[a:*]
;
;track_date=strmid(str,9,19)
;track_y=float(strmid(str,30,5))
;inst=strmid(str,35,4)
;
;final structure containing the times, elongations, instrument, spacecraft
;tr = {track_date: track_date, track_y: track_y, inst: inst, sc:spacecraft}
;---------------------------------------

if KEYWORD_SET(single) then single=1 else single=0

if sc eq 'Solar Orbiter' then begin
  dir = '/nas/helio/data/SolarOrbiter/HItracks/mabauer/'
  file_in = dir + start + '.sav'
  bflag=''
endif else begin

  if not single then begin
    dir = '/nas/helio/data/STEREO/HItracks/mabauer/'
    file_in = dir + start + '_' + sc + '_' + bflag + '.sav'
  endif
  if single then begin
    dir = '/nas/helio/data/STEREO/HItracks/mabauer/Single/'
    file_in = dir + start + '_A_' + bflag + tr_num + '.sav'
  endif
endelse

restore, file_in

radtodeg=180./!dpi;

tr=track

;cu=15

time=tr.track_date
epsilon=tr.elon;track_y
scname=tr.sc;
;tr.date=0

if sc eq 'Solar Orbiter' then scname = 'Solar Orbiter'

; set parameters for the imaging observatory
if scname eq 'A' then begin
    othersc='B'
    angle_signum=1;
     ;IF KEYWORD_SET(dir) THEN BEGIN
     ; initialphi=dir;
     ;ENDIF ELSE BEGIN
    initialphi=60;
endif

if scname eq 'B' then begin
    othersc='A'
    angle_signum=1
    ;IF KEYWORD_SET(dir) THEN BEGIN
    ;  initialphi=dir;
    ; ENDIF ELSE BEGIN
    initialphi=+60;
endif

if scname eq 'Solar Orbiter' then begin
    othersc='A'
    angle_signum=1
    ;IF KEYWORD_SET(dir) THEN BEGIN
    ;  initialphi=dir;
    ; ENDIF ELSE BEGIN
    initialphi=+60;
endif
;;;___________________________________________-

;use correct sign for A (-1) and B (+1)
epsilon=epsilon*angle_signum



; IF KEYWORD_SET(cut) THEN BEGIN
;     time2=strarr(s[1]-cut[0]-cut[1])
;     epsilon2=fltarr(s[1]-cut[0]-cut[1])
;     time2=time(0+cut[0]:s[1]-cut[1]-1)
;     epsilon2=epsilon(0+cut[0]:s[1]-cut[1]-1)
;     time=0
;     epsilon=0
;     ;    print, time2, epsilon2, size(time2)
;     time=time2
;     epsilon=epsilon2
;     print, 'Data points cutted from track for fitting:', cut
;     print, '   '
; ENDIF

if keyword_set(cut) then begin
  time = time[cut[0]:cut[1]]
  epsilon = epsilon[cut[0]:cut[1]]
endif

s=size(time)

;get STEREO positions
pos1=get_stereo_lonlat(time(s[1]/2), scname, system='HEE')
pos2=get_stereo_lonlat(time(s[1]/2), scname, system='HEEQ')

posother=get_stereo_lonlat(time(s[1]/2), othersc, system='HEE')
posother2=get_stereo_lonlat(time(s[1]/2), othersc, system='HEEQ')

;get position of the Earth
posearth=get_stereo_lonlat(time(s[1]/2), 'earth', system='HEE')


print,'Spacecraft positions on ', time(s[1]/2)
print,'   '
print,'IMAGING OBSERVATORY:'
print,'Position of ', scname,' in HEEQ long/lat:',pos1[1]*radtodeg, pos1[2]*radtodeg
print,'Position of ', scname,' in HEE long/lat:',pos2[1]*radtodeg, pos2[2]*radtodeg
print,'  '
print,'IN SITU OBSERVATORY:'
print,'Position of STEREO-',othersc,' in HEEQ long/lat:',posother[1]*radtodeg, posother[2]*radtodeg
print,'Position of STEREO-',othersc,' in HEE long/lat:',posother2[1]*radtodeg, posother2[2]*radtodeg
print,'   '
print,'elongation of first track data point:  ',min(abs(epsilon)), ' on ', time[0], ' UT'
print,'elongation of last track data point:   ',max(abs(epsilon)), ' on ', time[s[1]-1], ' UT'

print, '         '


tes=anytim(time)-anytim(time[0]); time to seconds for fitting

;convert time back from seconds to UTC
tzero=anytim(time[0])



;________________________________________________
;---------------FITTING

degtorad=!dpi/180;
H0=pos1[0];Imaging observatory to Sun in km
;------------------
;starting points of the minimization
phi=initialphi*degtorad;
vr=500;
tinit=tes[0]-86400; launch time guess: first HI1 observation minus one day
;-------------------

;Fixed Phi fitting
X=[phi, vr, tinit];
common myfit,xueber,yueber, dst
xueber=tes;
yueber=epsilon;
dst=H0;
RES=0;
RES=AMOEBA(1.0e-5,FUNCTION_NAME='fitj',FUNCTION_VALUE=values,P0=X,scale=[1,100,1e4]);


;Harmonic Mean fitting
X=[phi, vr, tinit];
common myfit2, xueber2, yueber2, dst2, scnameueber2
xueber2=tes;
yueber2=epsilon;
scnameueber2=scname;
dst2=H0;
RESH=0;
RESH=AMOEBA(1.0e-5,FUNCTION_NAME='fith2',FUNCTION_VALUE=values,P0=X,scale=[1,100,1e4]);




X=[phi, vr, tinit];
common myfit3, xueber3, yueber3, dst3, scnameueber3, lambda
xueber3=tes;
yueber3=epsilon;
scnameueber3=scname;
dst3=H0;
;fix width in sse fits from beginning
IF KEYWORD_SET(setlambda) THEN BEGIN
       lambda=setlambda*degtorad;
     ENDIF ELSE BEGIN
       lambda=45*degtorad;
     ENDELSE

;SSE fitting
RESS=0;
RESS=AMOEBA(1.0e-5,FUNCTION_NAME='fitsse',FUNCTION_VALUE=values,P0=X,scale=[1,100,1e4]);

;FITTING RESULTS

print, 'STARTING POINT ', X(0)/degtorad, X(1), X(2)
print, 'Fixed- Phi Minimization:    ', RES[0]/degtorad,'  ', RES[1],  '   ',RES[2]
print, 'Harmonic Mean Minimization: ', RESH[0]/degtorad,'  ', RESH[1],  '   ',RESH[2]
print, 'SSE Minimization:           ', RESS[0]/degtorad,'  ', RESS[1],  '   ',RESS[2]
print, 'SSE half width lambda:      ', lambda*radtodeg



;**********************************************************
;------------- FP DIRECTIONS
earthangle=RES[0]*radtodeg+pos1(1)*radtodeg
otherangle=(pos1(1)*radtodeg-posother(1)*radtodeg )+RES[0]*radtodeg
print, '          '
print, 'Direction: positive longitude angles -> solar west (right)'
print, 'FP:'
print, 'Angle to Earth   ', num2str(earthangle, FORMAT='(F10.2)')
print, 'Angle to STEREO-',othersc,'   ', num2str(otherangle,FORMAT='(F10.2)')
print, '  '

;------------- HM DIRECTIONS
earthangleh=RESH[0]*radtodeg+pos1(1)*radtodeg
otherangleh=(pos1(1)*radtodeg-posother(1)*radtodeg )+RESH[0]*radtodeg
print, 'HM:'
print, 'Angle to Earth   ', num2str(earthangleh, FORMAT='(F10.2)')
print, 'Angle to STEREO-',othersc,'  ', num2str(otherangleh,FORMAT='(F10.2)')
print, ' '


;------------- SSE DIRECTIONS
earthangles=RESS[0]*radtodeg+pos1(1)*radtodeg
otherangles=(pos1(1)*radtodeg-posother(1)*radtodeg )+RESS[0]*radtodeg
print, 'SSE: width lambda=', num2str(lambda*radtodeg,FORMAT='(F10.1)')
print, 'ICME Angle to Earth   ', num2str(earthangles, FORMAT='(F10.2)')
print, 'ICME Angle to STEREO-',othersc,'  ', num2str(otherangles,FORMAT='(F10.2)')
;************************************************
print, ' '



;------FP Arrival time calculation
;launch time
tl0=RES(2)

;Arrival time: constant speed Sun to 1 AU, launch on Sun center
arrivaltime_st_sec=tzero+tl0+posother(0)/RES[1]; in seconds
arrst=strmid(anytim(arrivaltime_st_sec, /vms),0,17); as date

arrivaltime_earth_sec=tzero+tl0+posearth(0)/RES[1]
arrearth=strmid(anytim(arrivaltime_earth_sec, /vms),0,17)

arrivaltime_ace_sec=arrivaltime_earth_sec-1.5*1e6/RES[1];
arrace=strmid(anytim(arrivaltime_ace_sec, /vms),0,17)


;------HM Arrival time calculation (Mšstl et al. 2011)
;launch time
tl0h=RESH(2)

;check if the spacecraft is hit at all: -90, +90 around direction
if  abs(otherangleh) lt 90 then begin
 ViHM_other=RESH[1]*cos(otherangleh/radtodeg)
 arrivaltime_st_sech=tzero+tl0h+posother(0)/(ViHM_other); in seconds
 arrsth=strmid(anytim(arrivaltime_st_sech, /vms),0,17); as date
endif else begin
 print, 'STEREO-',othersc,' is not hit by the HM circle'
endelse

;check if L1/Earth are hit: -90, +90 around direction
if abs(earthangleh) lt 90 then begin
 ;speed correction a la Moestl et al. 2011: flank speed is apex speed x cos(delta)
 ViHM_earth=RESH[1]*cos(earthangleh/radtodeg)
 ;with the speed correction there is no need for extra
 ;arrival time correction (Moestl et al. 2011, ApJ)
 arrivaltime_earth_sech=tzero+tl0h+posearth(0)/(ViHM_earth)
 arrearthh=strmid(anytim(arrivaltime_earth_sech, /vms),0,17)
 arrivaltime_ace_sech=arrivaltime_earth_sech-1.5*1e6/(ViHM_earth);
 arraceh=strmid(anytim(arrivaltime_ace_sech, /vms),0,17)
endif else begin
 print, 'Earth is not hit by the HM circle'
endelse


;------SSE Arrival time calculation  (Moestl and Davies, 2012)
;launch time
tl0s=RESS(2)
;other STEREO spacecraft
;check if the spacecraft is hit at all: spacecraft in between -lambda, +lambda of apex direction
 if abs(otherangles) lt (lambda*radtodeg)  then begin
   ; Arrival time correction is implemented through speed correction
   ; so the speed of the apex is is replaced with the arrival speed; see
   ; Equation 18 in Moestl and Davies, 2012, Solar Physics
   ViSSE_other=(RESS[1]*(cos(otherangles*degtorad)+sqrt(sin((lambda))^2-sin((otherangles*degtorad))^2))/(1+sin(lambda)));
   arrivaltime_st_secs=tzero+tl0s+posother(0)/(ViSSE_other); in seconds
   arrsts=strmid(anytim(arrivaltime_st_secs, /vms),0,17); arrival time as date
   endif else begin
  print, 'STEREO-',othersc,' is not hit by the SSE circle'
 endelse

;arrival time for Earth/L1
 if abs(earthangles) lt (lambda*radtodeg) then begin
   ViSSE_earth=(RESS[1]*(cos(earthangles*degtorad)+sqrt(sin((lambda))^2-sin((earthangles*degtorad))^2))/(1+sin(lambda)));
   arrivaltime_earth_secs=tzero+tl0h+posearth(0)/(ViSSE_earth)
   arrearths=strmid(anytim(arrivaltime_earth_secs, /vms),0,17);
   arrivaltime_ace_secs=arrivaltime_earth_secs-1.5*1e6/(RESS[1]*cos(earthangles/radtodeg));
   arraces=strmid(anytim(arrivaltime_ace_secs, /vms),0,17)
 endif else begin
  print, 'Earth is not hit by the SSE circle'
 endelse
;---------------------------------------



print, '         '




print, 'RESULTS SUMMARY'
print, '         '
;_____________________________________________

print, 'FPF: '
print, 'phi =     ', num2str(round(RES[0]*radtodeg),FORMAT='(I10)')
print, 'speed =   ', num2str(round(RES[1]),FORMAT='(I10)')
launchtime=strmid(anytim(tzero+tl0, /vms),0,17)
print, 'launch time solar center:  ', launchtime
print, 'arrivaltime L1:           ', arrace
print, 'arrivaltime Earth:         ', arrearth
print, 'arrivaltime STEREO-',othersc,':      ', arrst
elements=size(yueber)
normalized_residue=fitj(RES)/elements(1)
print, 'fitting residue FP: ', num2str(normalized_residue,FORMAT='(F10.5)')
print, '  '

;_____________________________________________

print, 'HMF: '
print, 'phi =     ', num2str(round(RESH[0]*radtodeg),FORMAT='(I10)')
print, 'speed =   ', num2str(round(RESH[1]),FORMAT='(I10)')
launchtimeh=strmid(anytim(tzero+tl0h, /vms),0,17)
print, 'launch time solar center:  ', launchtimeh
;;;;check if the Earth is hit by the HM circle
if abs(earthangleh) lt 90 then begin
  print, 'arrivaltime L1:           ', arraceh
  print, 'arrival speed L1:         ', num2str(round(ViHM_earth),FORMAT='(I10)')
  print, 'arrivaltime Earth:         ', arrearthh
  endif else begin
   arrearthh='no hit';
   arraceh='no hit';
   ViHM_earth='no hit';
   print, 'L1 and Earth are not hit by the transient.  '
  endelse
;;;;;check if the other STEREO is hit by the HM circle
if abs(otherangleh) lt 90 then begin
   print, 'arrivaltime STEREO-',othersc,':      ', arrsth
   print, 'arrival speed STEREO-',othersc,':    ', num2str(round(ViHM_other),FORMAT='(I10)')
   endif else begin
     arrsth='no hit';
     ViHM_other='no hit';
   endelse

normalized_residueh=fith2(RESH)/elements(1)
print, 'fitting residue HM: ',num2str(normalized_residueh,FORMAT='(F10.5)')
print, '  '



;_____________________________________________


print, 'SSEF with half width = ', num2str(lambda*radtodeg, FORMAT='(F10.1)'),' degree'
print, 'phi =     ', num2str(round(RESS[0]*radtodeg),FORMAT='(I10)')
print, 'speed =   ', num2str(round(RESS[1]),FORMAT='(I10)')
launchtimes=strmid(anytim(tzero+tl0s, /vms),0,17)
print, 'launch time solar center:  ', launchtimes
;;;;;check if Earth is hit by the SSE circle
if abs(earthangles) lt (lambda*radtodeg) then begin
  print, 'arrivaltime L1:           ', arraces
  print, 'arrival speed L1:         ', num2str(round(ViSSE_earth),FORMAT='(I10)')
  print, 'arrivaltime Earth:         ', arrearths
  endif else begin
   arrearths='no hit';
   arraces='no hit';
   ViSSE_earth='no hit';
   print, 'L1 and Earth are not hit by the transient.  '
  endelse
;;;;;check if the other STEREO is hit by the SSE circle
if abs(otherangles) lt (lambda*radtodeg) then begin
   print, 'arrivaltime STEREO-',othersc,':      ', arrsts
   print, 'arrival speed STEREO-',othersc,':      ', num2str(round(ViSSE_other),FORMAT='(I10)')
   endif else begin
     arrsts='no hit';
     ViSSE_other='no hit';
     print, 'STEREO-',othersc,' is not hit by the transient.'
   endelse
 normalized_residues=fitsse(RESS)/elements(1)
print, 'fitting residue SSE: ', num2str(normalized_residues,FORMAT='(F10.5)')
print, '  '






;_________________________________ PLOT OUTPUT


;make interpolated time array so the fitting curve is smooth
timeinterp=findgen(101)
tesinterp=findgen(101)
timestep=(anytim(time(elements(1)-1))-anytim(time(0)))/100;
 for i=0,100 do begin
  timeinterp(i)=anytim(time(0))+i*timestep;
  tesinterp(i)=i*timestep;
 endfor

;FP: synthetic elongation time function from fits

delta=90*degtorad-RES[0]
;rho=RES[1]*(tes-RES[2])/H0;
rho=RES[1]*(tesinterp-RES[2])/H0;
fitmin=atan(rho*cos(delta)/(1-rho*sin(delta) ) )/degtorad


;HM: synthetic elongation time function from fits (Mšstl et al. 2011)
phi=RESH[0]
;because of the problems with the acosine one needs to switch the beta back to negative
;so that its positive in fith.pro
if scname eq 'A' then phi=-RESH[0]
v=RESH[1]
t=(tesinterp-RESH[2])
a=(2*dst)/(v*t)-cos(phi)
b=sin(phi)
fitminh=-acos((-b+a*sqrt(a^2+b^2-1))/(a^2+b^2) ) /degtorad


;SSE synthetic function with Davies et al. 2012 formula:
phi=RESS[0]
if scname eq 'A' then phi=-RESS[0]
v=RESS[1]
t=(tesinterp-RESS[2])
c=sin(lambda)
a=((dst3*(1+c))/(v*t))-cos(phi)
b=sin(phi)
fitmins=acos( (-b*c+a*sqrt((a^2+b^2-c^2)))/((a^2+b^2))      ) /degtorad









;--------------------------------------------------------------------------
;draw figure in window for jpeg output

set_plot,'X'

window,1,  xsize=1500,ysize=400, retain=2,xpos=0,ypos=500
!p.multi=[0,3,1]
!p.background=255
!p.color=0


;********left panel fixed phi
 ymax=max(tr.elon)+10

 ;plotting time - elongation
 utplot, time, abs(epsilon), psym=4, ytitle='Elongation in degree', charsize=1.5,  $
 background=255,color=0,title=[scname+' HI1/2 FP track fitting'],yrange=[0,round(ymax)], $
 _extra = extra_keywords;
 ;plotting the fit
 outplot, anytim(timeinterp,/vms), abs(fitmin), color=0, thick=2, psym=0;, xrange=[xmin, xmax],yrange=[0,ymax]

 xyouts,tes(2),ymax-4,['Angle to Earth = '$
 +num2str(earthangle,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=1.2


xyouts,tes(2),ymax-8,['to '+scname+' = '$
 +num2str(RES[0]/degtorad,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=1.2

xyouts,tes(2),ymax-12,['to STEREO-'+othersc+' = '$
 +num2str(otherangle,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=1.2


xyouts,tes(2),ymax-17,['V = '+num2str(round(RES[1]),FORMAT='(I14)')+ $
' km/s'], /DATA, Color=0, charsize=1.2

xyouts,tes(2),ymax-23,['launch time: '+launchtime], /DATA, Color=0, charsize=1.2

xyouts,tes(2),ymax-26,['Arrival L1: '+arrace], /DATA, Color=0, charsize=1.2

xyouts,tes(2),ymax-29,['Arrival ST'+othersc+': '+arrst], /DATA, Color=0, charsize=1.2

xyouts,tesinterp(n_elements(tesinterp)-10),6,$
  ['Fitting Residue: '+num2str(normalized_residue, format='(F5.3)')], $
  /DATA, Color=0, charsize=1.2, alignment=1

;plot HI1 borders
hi12border=[18.7, 18.7, 18.7]
hi1border=[4, 4, 4]
hi1borderup=[24, 24, 24]
outplot, time(0:3), hi12border, linestyle=2, color=0
outplot, time(0:3), hi1border, linestyle=1, color=0
outplot, time(0:3), hi1borderup, linestyle=1, color=0


;********middle panel harmonic mean
 ;plotting time - elongation
 utplot, time, abs(epsilon), psym=4, ytitle='Elongation in degree', charsize=1.5,  $
 background=255,color=0,title=[scname+' HI1/2 HM track fitting'],yrange=[0,round(ymax)], $
 _extra = extra_keywords;
 ;plotting the fit
 outplot, anytim(timeinterp,/vms), abs(fitminh), color=0, thick=2, psym=0;, xrange=[xmin, xmax],yrange=[0,ymax]

 ;Angles
xyouts,tes(2),ymax-4,['Angle to Earth = '$
 +num2str(earthangleh,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=1.2
xyouts,tes(2),ymax-8,['to '+scname+' = '$
 +num2str(RESH[0]/degtorad,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=1.2
xyouts,tes(2),ymax-12,['to STEREO-'+othersc+' = '$
 +num2str(otherangleh,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=1.2
;apex speed and times
xyouts,tes(2),ymax-17,['V = '+num2str(round(RESH[1]),FORMAT='(I14)')+ $
' km/s'], /DATA, Color=0, charsize=1.2
xyouts,tes(2),ymax-23,['launch time: '+launchtimeh], /DATA, Color=0, charsize=1.2
xyouts,tes(2),ymax-26,['Arrival L1: '+arraceh], /DATA, Color=0, charsize=1.2
xyouts,tes(2),ymax-29,['Arrival ST'+othersc+': '+arrsth], /DATA, Color=0, charsize=1.2
xyouts,tesinterp(n_elements(tesinterp)-10),6,$
  ['Fitting Residue: '+num2str(normalized_residueh, format='(F5.3)')], $
  /DATA, Color=0, charsize=1.2, alignment=1

;plot HI1 borders
hi12border=[18.7, 18.7, 18.7]
hi1border=[4, 4, 4]
hi1borderup=[24, 24, 24]
outplot, time(0:3), hi12border, linestyle=2, color=0
outplot, time(0:3), hi1border, linestyle=1, color=0
outplot, time(0:3), hi1borderup, linestyle=1, color=0

;--------------------------------------------------------------

;********right panel SSE
 ;plotting time - elongation
 utplot, time, abs(epsilon), psym=4, ytitle='Elongation in degree', charsize=1.5,  $
 background=255,color=0,title=[scname+' HI1/2 SSE track fitting'],yrange=[0,round(ymax)], $
 _extra = extra_keywords;
 ;plotting the fit
 outplot, anytim(timeinterp,/vms), abs(fitmins), color=0, thick=2, psym=0;, xrange=[xmin, xmax],yrange=[0,ymax]

 ;Angles
xyouts,tes(2),ymax-4,['Angle to Earth = '$
 +num2str(earthangles,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=1.2
xyouts,tes(2),ymax-8,['to '+scname+' = '$
 +num2str(RESS[0]/degtorad,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=1.2
xyouts,tes(2),ymax-12,['to STEREO-'+othersc+' = '$
 +num2str(otherangles,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=1.2
;apex speed and times
xyouts,tes(2),ymax-17,['V = '+num2str(round(RESS[1]),FORMAT='(I14)')+ $
' km/s'], /DATA, Color=0, charsize=1.2
xyouts,tes(2),ymax-23,['launch time: '+launchtimes], /DATA, Color=0, charsize=1.2
xyouts,tes(2),ymax-26,['Arrival L1: '+arraces], /DATA, Color=0, charsize=1.2
xyouts,tes(2),ymax-29,['Arrival ST'+othersc+': '+arrsts], /DATA, Color=0, charsize=1.2
xyouts,tesinterp(n_elements(tesinterp)-10),6,$
  ['Fitting Residue: '+num2str(normalized_residues, format='(F5.3)')], $
  /DATA, Color=0, charsize=1.2, alignment=1

;plot HI1 borders
hi12border=[18.7, 18.7, 18.7]
hi1border=[4, 4, 4]
hi1borderup=[24, 24, 24]
outplot, time(0:3), hi12border, linestyle=2, color=0
outplot, time(0:3), hi1border, linestyle=1, color=0
outplot, time(0:3), hi1borderup, linestyle=1, color=0

;make jpeg output
;filejpg=dir+'/fitall.jpg'
;x2jpeg, filejpg



;---------GEOMETRY PLOT, output as jpeg and eps with filename +'_geometry'
;the average observation date is used
;fit_geometry, earthangle,earthangleh,earthangles,lambda*radtodeg,time(s(1)/2), file_in=file_in
;-------------------------------



;___________________________________________________

;make similar eps output

set_plot,'PS'

thickness=3
chsize=0.6;
device, /encapsulated, $
filename=dir+bflag+'_HIfits'+'.eps',$
xsize=30, ysize=30*400/1500, /color, bits_per_pixel=8

!p.multi=[0,3,1]
!p.background =255
!p.color = 0

;********left panel fixed phi
 ; plot elongation vs. time
 utplot, time, abs(epsilon), psym=4, ytitle='Elongation in degree', charsize=1.5,  $
    background=255,color=0,title=[scname+' HI1/2 FP track fitting'],yrange=[0,round(ymax)], $
    charthick=thickness, thick=thickness;, xrange=[launchtime, time(n_elements(time)-1)]
 ;plot for the fit
 outplot, anytim(timeinterp,/vms), abs(fitmin), color=0, thick=thickness, psym=0;, xrange=[xmin, xmax],yrange=[0,ymax]

 xyouts,tes(2),ymax-4,['Angle to Earth = '$
  +num2str(earthangle,FORMAT='(F10.1)')+' deg'], $
  /DATA, Color=0, charsize=chsize, charthick=thickness

 xyouts,tes(2),ymax-8,['to '+scname+' = '$
  +num2str(RES[0]/degtorad,FORMAT='(F10.1)')+' deg'], $
  /DATA, Color=0, charsize=chsize, charthick=thickness

 xyouts,tes(2),ymax-12,['to STEREO-'+othersc+' = '$
  +num2str(otherangle,FORMAT='(F10.1)')+' deg'], $
  /DATA, Color=0, charsize=chsize, charthick=thickness

 xyouts,tes(2),ymax-17,['V = '+num2str(round(RES[1]),FORMAT='(I14)')+ $
 ' km/s'], /DATA, Color=0, charsize=chsize, charthick=thickness

 xyouts,tes(2),ymax-23,['launch time: '+launchtime], /DATA, Color=0, charsize=chsize, charthick=thickness
 xyouts,tes(2),ymax-26,['Arrival L1: '+arrace], /DATA, Color=0, charsize=chsize, charthick=thickness
 xyouts,tes(2),ymax-29,['Arrival ST'+othersc+': '+arrst], /DATA, Color=0, charsize=chsize, charthick=thickness
 xyouts,tesinterp(n_elements(tesinterp)-10),6,['Fitting Residue: '+num2str(normalized_residue, $
  format='(F5.3)')], /DATA, Color=0, charsize=chsize, charthick=thickness,alignment=1

;indicate borders of HI1/HI2 field of view
hi12border=[ 18.7, 18.7, 18.7]
hi1border=[4, 4, 4]
hi1borderup=[24, 24, 24]
outplot, time(0:3), hi12border, linestyle=2, color=0, thick=thickness
outplot, time(0:3), hi1border, linestyle=1, color=0, thick=thickness
outplot, time(0:3), hi1borderup, linestyle=1, color=0, thick=thickness




;********middle panel harmonic mean

 ;plotting time - elongation
 utplot, time, abs(epsilon), psym=4, ytitle='Elongation in degree', charsize=1.5,  $
 background=255,color=0,title=[scname+' HI1/2 HM track fitting'],yrange=[0,round(ymax)], $
 _extra = extra_keywords, charthick=thickness, thick=thickness
 ;plotting the fit
 outplot, anytim(timeinterp,/vms), abs(fitminh), color=0, thick=thickness, psym=0;, xrange=[xmin, xmax],yrange=[0,ymax]


 xyouts,tes(2),ymax-4,['Angle to Earth = '$
 +num2str(earthangleh,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=chsize, charthick=thickness


xyouts,tes(2),ymax-8,['to '+scname+' = '$
 +num2str(RESH[0]/degtorad,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=chsize, charthick=thickness

xyouts,tes(2),ymax-12,['to STEREO-'+othersc+' = '$
 +num2str(otherangleh,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=chsize,charthick=thickness


xyouts,tes(2),ymax-17,['V = '+num2str(round(RESH[1]),FORMAT='(I14)')+ $
' km/s'], /DATA, Color=0, charsize=chsize, charthick=thickness

xyouts,tes(2),ymax-23,['launch time: '+launchtimeh], /DATA, Color=0, charsize=chsize, charthick=thickness
xyouts,tes(2),ymax-26,['Arrival L1: '+arraceh], /DATA, Color=0, charsize=chsize,charthick=thickness
xyouts,tes(2),ymax-29,['Arrival ST'+othersc+': '+arrsth], /DATA, Color=0, charsize=chsize, charthick=thickness
xyouts,tesinterp(n_elements(tesinterp)-10),6,['Fitting Residue: '+num2str(normalized_residues, $
  format='(F5.3)')], /DATA, Color=0, charsize=chsize, charthick=thickness,alignment=1

;plot HI1 borders
hi12border=[18.7, 18.7, 18.7]
hi1border=[4, 4, 4]
hi1borderup=[24, 24, 24]
outplot, time(0:3), hi12border, linestyle=2, color=0, thick=thickness
outplot, time(0:3), hi1border, linestyle=1, color=0, thick=thickness
outplot, time(0:3), hi1borderup, linestyle=1, color=0, thick=thickness



;********right panel self-similar expansion (SSE)

 ;plotting time - elongation
 utplot, time, abs(epsilon), psym=4, ytitle='Elongation in degree', charsize=1.5,  $
 background=255,color=0,title=[scname+' HI1/2 SSE track fitting'],yrange=[0,round(ymax)], $
 _extra = extra_keywords, charthick=thickness, thick=thickness
 ;plotting the fit
 outplot, anytim(timeinterp,/vms), abs(fitmins), color=0, thick=thickness, psym=0;, xrange=[xmin, xmax],yrange=[0,ymax]

xyouts,tes(2),ymax-4,['Angle to Earth = '$
 +num2str(earthangles,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=chsize, charthick=thickness

xyouts,tes(2),ymax-8,['to '+scname+' = '$
 +num2str(RESS[0]/degtorad,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=chsize, charthick=thickness

xyouts,tes(2),ymax-12,['to STEREO-'+othersc+' = '$
 +num2str(otherangles,FORMAT='(F10.1)')+' deg'], $
 /DATA, Color=0, charsize=chsize,charthick=thickness

xyouts,tes(2),ymax-17,['V = '+num2str(round(RESS[1]),FORMAT='(I14)')+ $
' km/s'], /DATA, Color=0, charsize=chsize, charthick=thickness

xyouts,tes(2),ymax-23,['launch time: '+launchtimes], /DATA, Color=0, charsize=chsize, charthick=thickness
xyouts,tes(2),ymax-26,['Arrival L1: '+arraces], /DATA, Color=0, charsize=chsize,charthick=thickness
xyouts,tes(2),ymax-29,['Arrival ST'+othersc+': '+arrsts], /DATA, Color=0, charsize=chsize, charthick=thickness
xyouts,tesinterp(n_elements(tesinterp)-10),6,['Fitting Residue: '+num2str(normalized_residues, $
  format='(F5.3)')], /DATA, Color=0, charsize=chsize, charthick=thickness,alignment=1

;plot HI1 borders
hi12border=[18.7, 18.7, 18.7]
hi1border=[4, 4, 4]
hi1borderup=[24, 24, 24]
outplot, time(0:3), hi12border, linestyle=2, color=0, thick=thickness
outplot, time(0:3), hi1border, linestyle=1, color=0, thick=thickness
outplot, time(0:3), hi1borderup, linestyle=1, color=0, thick=thickness

device, /close




;make a structure with the results
results={ $
         ;Fixed-Phi --------------------------
         fp_launchtime:launchtime,           $
         fp_arrival_earth:arrearth,          $
   fp_arrival_L1:arrace,               $
         fp_arrival_stereo:arrst,            $
   fp_angle_observer:RES(0)*radtodeg,  $
   fp_angle_earth: earthangle,         $
         fp_angle_stereo:otherangle,         $
   fp_speed:RES(1),                    $
   fp_residue:normalized_residue,      $
         ;Harmonic Mean  ---------------------
         hm_launchtime:launchtimeh,          $
         hm_arrival_earth:arrearthh,         $
         hm_arrival_L1:arraceh,              $
         hm_arrival_stereo:arrsth,           $
   hm_angle_observer:RESH(0)*radtodeg, $
   hm_angle_earth:earthangleh,         $
         hm_angle_stereo:otherangleh,        $
         hm_speed_apex:RESH(1),              $
         hm_speed_earth:ViHM_earth,          $
         hm_speed_stereo:ViHM_other,         $
         hm_residue:normalized_residueh,     $
         ;SSE --------------------------------
   sse_width:lambda*radtodeg,          $
         sse_launchtime:launchtimes,         $
         sse_arrival_earth:arrearths,        $
   sse_arrival_L1:arraces,             $
         sse_arrival_stereo:arrsts,          $
   sse_angle_observer:RESS(0)*radtodeg,$
   sse_angle_earth:earthangles,        $
         sse_angle_stereo:otherangles,       $
         sse_speed_apex:RESS(1),             $
         sse_speed_earth:ViSSE_earth,        $
         sse_speed_stereo:ViSSE_other,       $
         sse_residue:normalized_residues }


set_plot,'X'
!p.multi=[0,1,1]
!p.background = 0
!p.color = 255



; ;make txt output file
; ;-------------------------------------------------------------------------------

; journal_file='results_'+strmid(file_in,0,53)+'.txt';

; print, 'Results written in current directory to txt file:', journal_file

; openw, 1, journal_file

; printf,1,'Spacecraft positions '
; printf,1,'Mean position on ', tr.track_date(s(1)/2) , ' (average date of the elongation values)'
; printf,1,'IMAGING OBSERVATORY:'
; printf,1,'STEREO-', scname, ' in HEEQ long/lat:', pos1(1)*radtodeg, pos1(2)*radtodeg
; printf,1,'STEREO-', scname, ' in HEE long/lat:',  pos2(1)*radtodeg, pos2(2)*radtodeg
; printf,1,'       '
; printf,1,'IN SITU OBSERVATORY:'
; printf,1,'STEREO-', othersc, ' in HEEQ long/lat:',  posother(1)*radtodeg, posother(2)*radtodeg
; printf,1,'STEREO-', othersc, ' in HEE long/lat:',  posother2(1)*radtodeg, posother2(2)*radtodeg
; printf,1,'   '
; printf,1,'STEREO separation HEE   = ', (abs(posother2(1))+abs(pos2(1)))*radtodeg,' degree '
; printf,1,'       '
; printf,1,'elongation of first track data point:  ',min(abs(epsilon)), ' degree on ', tr.track_date(0), ' UT'
; printf,1,'elongation of last track data point:   ',max(abs(epsilon)), ' degree on ', tr.track_date(s(1)-1), ' UT'
; printf,1,'       '
; printf,1,'-------------------------------------       '

; printf, 1,'FPF results'
; printf, 1,'launch Sun center:   ',results.fp_launchtime, ' UT'
; printf, 1,'arrival Earth:       ',results.fp_arrival_earth, ' UT'
; printf, 1,'arrival L1:          ',results.fp_arrival_L1, ' UT'
; printf, 1,'arrival STEREO:      ',results.fp_arrival_stereo, ' UT'
; printf, 1,'angle observer:      ',num2str(results.fp_angle_observer, FORMAT='(F10.1)'), ' degree'
; printf, 1,'angle Earth/L1:      ',num2str(results.fp_angle_earth, FORMAT='(F10.1)'),' degree'
; printf, 1,'angle STEREO:        ',num2str(results.fp_angle_stereo, FORMAT='(F10.1)'), ' degree'
; printf, 1,'speed:               ',results.fp_speed, ' km/s'
; printf, 1,'fit residue:         ',results.fp_residue
; printf, 1,'   '


; printf, 1,'HMF results'
; printf, 1,'launch Sun center:   ',results.hm_launchtime, ' UT'
; printf, 1,'arrival Earth:       ',results.hm_arrival_earth, ' UT'
; printf, 1,'arrival L1:          ',results.hm_arrival_L1, ' UT'
; printf, 1,'arrival STEREO:      ',results.hm_arrival_stereo, ' UT'
; printf, 1,'angle observer:      ',num2str(results.hm_angle_observer, FORMAT='(F10.1)'), ' degree'
; printf, 1,'angle Earth:         ',num2str(results.hm_angle_earth, FORMAT='(F10.1)'),' degree'
; printf, 1,'angle STEREO:        ',num2str(results.hm_angle_stereo, FORMAT='(F10.1)'), ' degree'
; printf, 1,'speed apex:          ',results.hm_speed_apex, ' km/s'
; printf, 1,'speed Earth/L1:      ',results.hm_speed_earth, ' km/s'
; printf, 1,'speed STEREO  :      ',results.hm_speed_stereo, ' km/s'
; printf, 1,'fit residue:         ',results.hm_residue
; printf, 1,'  '

; printf, 1,'SSEF results'
; printf, 1,'half width:          ',num2str(results.sse_width, FORMAT='(F10.2)'),' degree'
; printf, 1,'launch Sun center:   ',results.sse_launchtime, ' UT'
; printf, 1,'arrival Earth:       ',results.sse_arrival_earth, ' UT'
; printf, 1,'arrival L1:          ',results.sse_arrival_L1, ' UT'
; printf, 1,'arrival STEREO:      ',results.sse_arrival_stereo, ' UT'
; printf, 1,'angle observer:      ',num2str(results.sse_angle_observer, FORMAT='(F10.1)'), ' degree'
; printf, 1,'angle Earth:         ',num2str(results.sse_angle_earth, FORMAT='(F10.1)'),' degree'
; printf, 1,'angle STEREO:        ',num2str(results.sse_angle_stereo, FORMAT='(F10.1)'), ' degree'
; printf, 1,'speed apex:          ',results.sse_speed_apex, ' km/s'
; printf, 1,'speed Earth/L1:      ',results.sse_speed_earth, ' km/s'
; printf, 1,'speed STEREO  :      ',results.sse_speed_stereo, ' km/s'
; printf, 1,'fit residue:         ',results.sse_residue
; close, 1

;---------------------------------------------------------------------------------

;return the structure with all results
return, results


end


