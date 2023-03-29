;+
; general ellipse evolution model
;
; Authors: Christian Moestl, Tanja Amerstorfer
; last update 2015 Feb 23
;
; takes input from parameter file
;
; example call: elevo, 'filein.txt'  filein.txt is in folder elevo_events
;
; output:
;       pred: Prediction for the different targets
;       elevo_kin: kinematics of the differents runs
;
; calls elevo_elong
;       elliptical_geometry_for_movies_mars_new
;       elevo_analytic
;       SPICE
;
; History:    2019/08: added elevo_kin (Juergen Hinterreiter)
;-



pro elevo, dir, pred, elevo_kin, runnumber

common vinit, gammaparam, background_wind

;constants
AU=149597871 ;km
r_sun=6.957d5; km

print, '                      '
print, '============================'
print, 'ELEvo CME arrival prediction'
print, '============================'

;************general controls************

;for movie
startframe=0;
endframe=400;

;read parameters from control file
fnam=dir+'elevo_input.txt'

str = STRARR(200)
OPENR, 10, fnam
   dummy = ''
   i = 0L
   WHILE NOT(EOF(10)) DO BEGIN
      READF, 10, dummy
      str[i] = dummy
      i = i + 1
   ENDWHILE
CLOSE, 10

elevo_plot_title=str[22]

;this is the initial time for the CME, used for spacecraft positions + initial time for drag model

suntime=str[12]

tinit=suntime

timegrid=1440  ;number of points for 10 days with 10 min resolution
;timegrid=2500

;equidistant grid for DBM times, with 10 min resolution
;time
tdrag=dblarr(timegrid)
;speed
vdrag=dblarr(n_elements(tdrag))
;distance
rdrag=dblarr(n_elements(tdrag))

for i=0, n_elements(tdrag)-1 do begin
   tdrag[i]=anytim(tinit[0])+i*10*60.
endfor

tdrag_str =  anytim(tdrag, out_style='vms')


soloAvailable = 0
pspAvailable = 0
bepiAvailable = 0

if anytim(suntime) gt anytim('01-Mar-2020 00:00:00') then soloAvailable = 1
if anytim(suntime) gt anytim('01-Jan-2019 00:00:00') then pspAvailable = 1
if (anytim(suntime) gt anytim('21-Oct-2018 00:00:00')) and (anytim(suntime) lt anytim('01-Nov-2025 00:00:00')) then bepiAvailable = 1


;spacecraft and planet positions (cartesian, lon/lat)
pos_sta=get_sunspice_lonlat(tdrag_str, 'A', system='HEE')
pos_sta_xy=get_sunspice_coord(tdrag_str, 'A', system='HEE', /AU)

pos_stb=get_sunspice_lonlat(tdrag_str, 'B', system='HEE')
pos_stb_xy=get_sunspice_coord(tdrag_str, 'B', system='HEE', /AU)

pos_earth=get_sunspice_lonlat(tdrag_str, 'Earth', system='HEE')
pos_earth_xy=get_sunspice_coord(tdrag_str, 'Earth', system='HEE', /AU)

pos_vex=get_sunspice_lonlat(tdrag_str, 'Venus', system='HEE')
pos_vex_xy=get_sunspice_coord(tdrag_str, 'Venus', system='HEE', /AU)

pos_mars=get_sunspice_lonlat(tdrag_str, 'Mars', system='HEE')
pos_mars_xy=get_sunspice_coord(tdrag_str, 'Mars', system='HEE', /AU)

if soloAvailable eq 1 then pos_solo=get_sunspice_lonlat(tdrag_str, 'solo', system='HEE')
if soloAvailable eq 1 then pos_solo_xy=get_sunspice_coord(tdrag_str, 'solo', system='HEE', /AU)

if pspAvailable eq 1 then pos_psp=get_sunspice_lonlat(tdrag_str, 'psp', system='HEE')
if pspAvailable eq 1 then pos_psp_xy=get_sunspice_coord(tdrag_str, 'psp', system='HEE', /AU)

if bepiAvailable eq 1 then pos_bepi=get_sunspice_lonlat(tdrag_str, 'BepiColombo-MPO', system='HEE')
if bepiAvailable eq 1 then pos_bepi_xy=get_sunspice_coord(tdrag_str, 'BepiColombo-MPO', system='HEE', /AU)

;Position of MESSENGER is Mercury after 01/04/2012, load CSPICE Kernel to determine position before
if anytim(suntime) gt anytim('2004-08-03T00:00:00') and anytim(suntime) lt anytim('2012-04-01T00:00:00') then begin

	utc = anytim2utc(tdrag_str, /ccsds) 

  file = concat_dir( getenv('STEREO_SPICE_OTHER'), $
                     'msgr_20040803_20120401_od051.bsp')
  get_stereo_spice_range, file, tai0, tai1, scid, /tai
  cspice_furnsh, file
  tai = utc2tai(utc)

  edate = tai2utc(tai0 > tai < tai1, /ccsds)

  if (tai[-1] ge tai0) and (tai[0] le tai1) then begin
    mess = get_stereo_coord(utc, /au, system='HEE', 'Messenger')

  end else begin
    mu = 2.2032*10^4
    state = get_stereo_coord(edate, 'Mercury', system='HAE')

    et = dblarr(n_elements(edate))

    for i=0, n_elements(edate)-1 do begin

      cspice_utc2et, edate[i], et_val
      et[i] = et_val

    endfor

    cspice_oscelt, state, et, mu, elts

    for i=0, n_elements(edate)-1 do begin

      cspice_utc2et, utc[i], et_val
      et[i] = et_val

    endfor

    cspice_conics, elts, et, mess
    convert_stereo_coord, utc, mess, 'HAE', 'HEE'

    mess = mess/au
  endelse

  cspice_unload, file

  cspice_reclat, mess[0:2, *], radius, longitude, latitude


  mesp_sz = size(mess)
  pos_mes_xy=fltarr(3, mesp_sz[2])
  pos_mes=fltarr(3, mesp_sz[2])

  pos_mes_xy[0, *]=mess[0, *]*au
  pos_mes_xy[1, *]=mess[1, *]
  pos_mes_xy[2, *]=mess[2, *]

  pos_mes[0, *]=radius*au
  pos_mes[1, *]=longitude
  pos_mes[2, *]=latitude
endif else begin
    pos_mes=get_sunspice_lonlat(tdrag_str, 'Mercury', system='HEE')
    pos_mes_xy = get_sunspice_coord(tdrag_str, 'Mercury', system='HEE', /AU)
endelse

;ellipse parameters
;direction of apex relative to Earth at tinit
halfwidth=float(str[2])
aspectratio=float(str[4])
direction=float(str[6])
f=1/aspectratio

;apex directions relative to SC
;can be positive or negative, depending on SC position
;delta is negative if CME is to the left (solar East) of the SC
;this takes into account that SC can switch from positive to negative longitude

delta_A = fltarr(n_elements(pos_sta[1, *]))

for i=0,n_elements(pos_sta[1, *])-1 do begin
  if SIGNUM(direction) EQ SIGNUM(pos_sta[1,i]) then delta_A[i] = direction-pos_sta[1, i]/!dtor
  if SIGNUM(direction) NE SIGNUM(pos_sta[1,i]) then delta_A[i] = direction-(pos_sta[1, i]/!dtor+360*SIGNUM(direction))
endfor

delta_B = fltarr(n_elements(pos_stb[1, *]))

for i=0,n_elements(pos_stb[1, *])-1 do begin
  if SIGNUM(direction) EQ SIGNUM(pos_stb[1,i]) then delta_B[i] = direction-pos_stb[1, i]/!dtor
  if SIGNUM(direction) NE SIGNUM(pos_stb[1,i]) then delta_B[i] = direction-(pos_stb[1, i]/!dtor+360*SIGNUM(direction))
endfor

delta_V = fltarr(n_elements(pos_vex[1, *]))

for i=0,n_elements(pos_vex[1, *])-1 do begin
  if SIGNUM(direction) EQ SIGNUM(pos_vex[1,i]) then delta_V[i] = direction-pos_vex[1, i]/!dtor
  if SIGNUM(direction) NE SIGNUM(pos_vex[1,i]) then delta_V[i] = direction-(pos_vex[1, i]/!dtor+360*SIGNUM(direction))
endfor

delta_MES = fltarr(n_elements(pos_mes[1, *]))

for i=0,n_elements(pos_mes[1, *])-1 do begin
  if SIGNUM(direction) EQ SIGNUM(pos_mes[1,i]) then delta_MES[i] = direction-pos_mes[1, i]/!dtor
  if SIGNUM(direction) NE SIGNUM(pos_mes[1,i]) then delta_MES[i] = direction-(pos_mes[1, i]/!dtor+360*SIGNUM(direction))
endfor

if soloAvailable eq 1 then begin
  delta_SOLO = fltarr(n_elements(pos_solo[1, *]))

  for i=0,n_elements(pos_solo[1, *])-1 do begin
  if SIGNUM(direction) EQ SIGNUM(pos_solo[1,i]) then delta_SOLO[i] = direction-pos_solo[1, i]/!dtor
  if SIGNUM(direction) NE SIGNUM(pos_solo[1,i]) then delta_SOLO[i] = direction-(pos_solo[1, i]/!dtor+360*SIGNUM(direction))
  endfor
endif

if pspAvailable eq 1 then begin
  delta_PSP = fltarr(n_elements(pos_psp[1, *]))

  for i=0,n_elements(pos_psp[1, *])-1 do begin
  if SIGNUM(direction) EQ SIGNUM(pos_psp[1,i]) then delta_PSP[i] = direction-pos_psp[1, i]/!dtor
  if SIGNUM(direction) NE SIGNUM(pos_psp[1,i]) then delta_PSP[i] = direction-(pos_psp[1, i]/!dtor+360*SIGNUM(direction))
  endfor
endif

if bepiAvailable eq 1 then begin
  delta_BEPI = fltarr(n_elements(pos_bepi[1, *]))

  for i=0,n_elements(pos_bepi[1, *])-1 do begin
  if SIGNUM(direction) EQ SIGNUM(pos_bepi[1,i]) then delta_BEPI[i] = direction-pos_bepi[1, i]/!dtor
  if SIGNUM(direction) NE SIGNUM(pos_bepi[1,i]) then delta_BEPI[i] = direction-(pos_bepi[1, i]/!dtor+360*SIGNUM(direction))
  endfor
endif

;delta_E constant -> position CME apex relative to Earth is constant
delta_E = fltarr(n_elements(pos_earth[1, *]))
  
for i=0,n_elements(pos_earth[1, *])-1 do begin
  delta_E[i] = direction
endfor

print, '------------------------------------'
print, 'Spacecraft separation from CME apex at tinit:'
print, 'STEREO-A:', delta_A[0], format='(A, 2x, F6.1)'
print, 'STEREO-B:', delta_B[0], format='(A, 2x, F6.1)'
print, 'Wind:', delta_E[0], format='(A, 5x, F6.1)'
print, 'MESSENGER: ', delta_MES[0], format='(A, 2x, F6.1)'
print, 'Venus', delta_V[0], format='(A, 5x, F6.1)'
print, '------------------------------------'

print, '------------------------------------'
print, 'Ellpise parameters:'
print, 'Direction from Earth: ', direction, format='(A, 1x, F5.1)'
print, 'Half width: ', halfwidth, format='(A, 12x, I2)'
print, 'Aspect ratio (a/b): ', 1./f, format='(A, 4x, F4.2)'
print, '------------------------------------'

;parameters needed for drag based model
rinit=float(str[10])            ;initial distance [Rsun]
vinit=float(str[14])            ;initial speed [km/s]
tinit=suntime                   ;initial time [UT]
background_wind=float(str[16])  ;background wind speed [km/s]
gammaparam=float(str[18])       ;E-07/km

print, '------------------------------------'
print, 'DBM parameters:'
print, 'Initial distance: ', rinit, format='(A, 12x, F6.2)'
print, 'Initial time: ', anytim(tinit, /vms), format='(A, 3x, A17)'
print, 'Initial speed: ', vinit, format='(A, 15x, F7.2)'
print, 'Drag parameter [E-07/km]: ', gammaparam, format='(A, 5x, F5.2)'
print, 'Background solar wind speed:', background_wind, format='(A, 2x, I4)'
print, '------------------------------------'

;create 1-D DBM kinematic for ellipse apex with
;constant drag parameter and constant background solar wind speed

;then use Vrsnak et al. 2013 equation 5 for v(t), 6 for r(t)

;acceleration or deceleration
;Note that the sign in the dbm equation is derived from fitting - not from comparing vinit to solar wind speed.
;It is possible that the sign in the equation is negative although vinit > sw_speed because all data points are taken into account - not only the initial speed.
accsign=1

;if vinit lt background_wind then accsign=-1 else accsign=1

tinitnum=anytim(tinit)

for i=0, n_elements(tdrag)-1 do begin

  ;distance in km
  rdrag[i]=(accsign/(gammaparam*1e-7))*alog(1+(accsign*(gammaparam*1e-7)*((vinit-background_wind)*(tdrag[i]-tinitnum))))+background_wind*(tdrag[i]-tinitnum)+Rinit*r_sun
  ;convert from km to AU
  rdrag[i]=rdrag[i]/AU;
  ;speed in km/s
  vdrag[i]=(vinit-background_wind)/(1+(accsign*(gammaparam*1e-7)*((vinit-background_wind)*(tdrag[i]-tinitnum))))+background_wind

  if finite(rdrag[i]) eq 0 then begin
    print, 'Sign of gamma does not fit to vinit and w! Check dbmfit.pro!'
    stop
  endif
endfor

steps=float(str[29]); steps are equidistant in time

f1=10
f2=f1+1*steps
f3=f1+2*steps
f4=f1+3*steps
f5=f1+4*steps

tdrag=anytim(tdrag, /vms)
t_plot=[tdrag[f1],tdrag[f2],tdrag[f3],tdrag[f4],tdrag[f5]]

color1 = 0
color2 = 50  ;blue
color3 = 170 ;green
color4 = 185 ;yellow
color5 = 120 ;red

colors=[color1,color2,color3,color4,color5]

;Turn -eps plot on/off
plot=1

if plot then begin 

  if runnumber eq 1 then begin

    !P.MULTI=0

    set_plot,'ps'

    fileplot=dir+str[27]

    device, filename=fileplot, xsize=15, ysize=15, /inches, /color, bits_per_pixel=8, /encapsulated
    print, 'plotting .eps figure...'
    loadct, 5, /silent
  endif
    ;decide which timesteps are plotted

    steps=float(str[29]); steps are equidistant in time

    f1=10
    f2=f1+1*steps
    f3=f1+2*steps
    f4=f1+3*steps
    f5=f1+4*steps

    framenumber=4 ;for which of the 0-4 frames the general parameters are given

    figure_frames=[f1,f2,f3,f4,f5]
    R_plot=[rdrag[f1],rdrag[f2],rdrag[f3],rdrag[f4],rdrag[f5]]
    V_plot=[vdrag[f1],vdrag[f2],vdrag[f3],vdrag[f4],vdrag[f5]]
    tdrag=anytim(tdrag, /vms)
    t_plot=[tdrag[f1],tdrag[f2],tdrag[f3],tdrag[f4],tdrag[f5]]

    s=size(figure_frames)

    ;plot ellipse propagating as eps:

  ;Earth distance in AU
  earth_dist=pos_earth/AU

  ;---------------------------------

  ;sun position with ssc_plot_where xrange and yrange from -1.1 to 1.1, and window 1000 by 1000
  sun=[0.5205,0.508]
  sun=[0.5205,0.509]

  ;earth position in normal coordinates
  earth=[sun[0], 0.245]

  ;scaling 1 AU to normal coordinates
  AUscale=sun[1]-earth[1]

  ;plot s/c positions

  if runnumber eq 1 then begin

    framenumber = 0
    ssc_plot_where_elevo, t_plot[framenumber], /WHITE_BG, xrange=[1.5, -1.5], yrange=[-1.5,1.5], /yst, /xst, thick=5, font=1, charsize=3, pos=[0.12,0.1,0.9,0.9], /norm, /mess, /psp, /solo

    ;spacecraft positions
    ;plots, [sun[0], earth[0]], [sun[1], earth[1]],  $
    ;	color=0, linestyle=0, /NORMAL, thick=2

    plots, [0, earth_dist[1]], [0, earth_dist[0]],  $
    	color=0, linestyle=0, /data, thick=2

    ;variable AUscale is distance in normal coordinates

    loadct, 5, /silent
    fpcolor = 120 ;red
    hmcolor = 50  ;blue
    ssecolor= 180 ;green

    ;set bounding and central line for ellipse similar to black for all
    loadct, 0, /silent
    elcolor=0
  endif

  earthangle=direction
  lambda=halfwidth

  sun=[0,0]

  ;ellipse directions
  ;central
  ;auscale=earthdist
  ssdir=[sin(earthangle*!dtor),cos(earthangle*!dtor)]*earth_dist[0]*(R_plot[s[1]-1])

  if runnumber eq 1 then begin

    plots, [sun[0], ssdir[0]], [sun[1], ssdir[1]],  $
    		color=elcolor, linestyle=0, /data, thick=5

    ;width positive
    ssdirplus=[sin((earthangle+lambda)*!dtor),cos((earthangle+lambda)*!dtor)]*earth_dist[0]*(R_plot[s[1]-1])
    plots, [sun[0], ssdirplus[0]], [sun[1], ssdirplus[1]],  $
    		color=elcolor, linestyle=0, /data, thick=5

    ;width negative
    ssdirminus=[sin((earthangle-lambda)*!dtor),cos((earthangle-lambda)*!dtor)]*earth_dist[0]*(R_plot[s[1]-1])
    plots, [sun[0], ssdirminus[0]], [sun[1], ssdirminus[1]],  $
    		color=elcolor, linestyle=0, /data, thick=5

    loadct, 5, /silent

   endif

  color1 = 0
  color2 = 50  ;blue
  color3 = 170 ;green
  color4 = 185 ;yellow
  color5 = 120 ;red

  colors=[color1,color2,color3,color4,color5]

  for i=0,s[1]-1  do begin

    ;set ellipse color different for each timestep
    elcolor=colors[i]

    ;draw ellipse, ;R_plot[i] is apex
    ;f=1/aspectratio  ;f=b/a
    theta=atan(f^2*tan(lambda*!dtor))
    omega=sqrt(cos(theta)^2*(f^2-1)+1)
    ;if this factor is set to other than 1 one can make very wide ellipses around the Sun
    ;if necessary
    factor=1  ; another possible free parameter
    b=R_plot[i]*omega*sin(lambda*!dtor)/(cos(lambda*!dtor-theta)+omega*sin(lambda*!dtor))*factor


    a=b/f
    c=R_plot[i]-b


    ;****speed and distance from Sun of given point along ellipse front****

    frame_ind = UINT(figure_frames[i])

    ;get distance and speed of point along delta of STEREO-A:
    dvalue=elevo_analytic(R_plot[i], aspectratio, halfwidth, delta_A[frame_ind])

    if finite(dvalue) then begin
  		deltaspeed_A=dvalue/R_plot[i]*V_plot[i]
    endif

    ;get distance and speed of point along delta of STEREO-B:
    dvalue=elevo_analytic(R_plot[i], aspectratio, halfwidth, delta_B[frame_ind])

    if finite(dvalue) then begin
  		deltaspeed_B=dvalue/R_plot[i]*V_plot[i]
    endif

    ;get distance and speed of point along delta of Venus:
    dvalue=elevo_analytic(R_plot[i], aspectratio, halfwidth, delta_V[frame_ind])

    if finite(dvalue) then begin
  		deltaspeed_V=dvalue/R_plot[i]*V_plot[i]
    endif

    ;get distance and speed of point along delta of MESSENGER:
    ;if anytim(suntime) gt anytim('2004-08-03T00:00:00') and anytim(suntime) lt anytim('2012-04-01T00:00:00') then begin
    dvalue=elevo_analytic(R_plot[i], aspectratio, halfwidth, delta_MES[frame_ind])

    if finite(dvalue) then begin
  	   deltaspeed_MES=dvalue/R_plot[i]*V_plot[i]
    endif
    ;endif
    ;stop
    ;get distance and speed of point along delta of Earth:
    dvalue=elevo_analytic(R_plot[i], aspectratio, halfwidth, delta_E[frame_ind])

    if finite(dvalue) then begin
  		deltaspeed_E=dvalue/R_plot[i]*V_plot[i]
    endif


    if soloAvailable eq 1 then begin
      ;get distance and speed of point along delta of Solar Orbiter:
      dvalue=elevo_analytic(R_plot[i], aspectratio, halfwidth, delta_SOLO[frame_ind])

      if finite(dvalue) then begin
        deltaspeed_SOLO=dvalue/R_plot[i]*V_plot[i]
      endif
    endif

    if pspAvailable eq 1 then begin
      ;get distance and speed of point along delta of Parker Solar Probe:
      ;delta_PSP = angle between CME apex and SC
      dvalue=elevo_analytic(R_plot[i], aspectratio, halfwidth, delta_PSP[frame_ind])

      if finite(dvalue) then begin
        deltaspeed_PSP=dvalue/R_plot[i]*V_plot[i]
      endif
    endif

    if bepiAvailable eq 1 then begin
      ;get distance and speed of point along delta of BEPI:
      dvalue=elevo_analytic(R_plot[i], aspectratio, halfwidth, delta_BEPI[frame_ind])

      if finite(dvalue) then begin
        deltaspeed_BEPI=dvalue/R_plot[i]*V_plot[i]
      endif
    endif


    ;get apex position minus b for ellipse center
    xangle=sin((earthangle)*!dtor)
    yangle=cos((earthangle)*!dtor)
    ellipse_center=[xangle,yangle]*c ;this is in data coordinates, so in AU
    ;draw this particular ellipse

    ;the sun is at data = 0,0 coordinates..
    if runnumber eq 1 then begin
      tvellipse_elevo, b, a, ellipse_center[0], ellipse_center[1], 90-earthangle, thick=4, /data, color=elcolor
    endif
  endfor
endif


lambda=halfwidth
theta=atan(f^2*tan(lambda*!dtor))
omega=sqrt(cos(theta)^2*(f^2-1)+1)

;if this factor is set to other than 1 one can make very wide ellipses around the Sun
;if necessary
factor=1  ; another possible free parameter

b = fltarr(N_ELEMENTS(rdrag))
a = fltarr(N_ELEMENTS(rdrag))
c = fltarr(N_ELEMENTS(rdrag))

for i=0,N_ELEMENTS(rdrag)-1 do begin
  b[i]=rdrag[i]*omega*sin(lambda*!dtor)/(cos(lambda*!dtor-theta)+omega*sin(lambda*!dtor))*factor

  a[i]=b[i]/f
  c[i]=rdrag[i]-b[i]
endfor

;*****add arrival times for each spacecraft that is hit by the ellipse*****

print, ' '
print, 'ARRIVAL TIMES at planets, spacecraft:'
print, '-------------------------------------'


tars=anytim(tdrag)-anytim(tdrag[0])

;calculate all distances of the ellipse in direction of MESSENGER

d_MES=elevo_analytic(rdrag, aspectratio, halfwidth, delta_MES)

v_MES=deriv(tars, d_MES*au)

;check where the heliocentric distance of STEREO-A is less than the ellipse distance
index_hit_MES=sc_ellipse_coords(direction, a, b, c, pos_mes_xy)
;take first value of these indices = arrival time at STEREO-A within drag time resolution (10 minutes)


if finite(index_hit_MES) then begin
  arrival_MES=anytim(tdrag[index_hit_MES[0]], /ccsds)
  arrival_speed_MES=v_MES[index_hit_MES[0]]
  print, '---------------------------------------------'
  print, 'Arrival at MESSENGER [UT]: ', anytim(arrival_MES, /vms), format='(A, 4x, A17)'
  print, 'Arrival speed at MESSENGER [km/s]: ', arrival_speed_MES, format='(A, 1x, I4)'
  print, '---------------------------------------------'
endif else begin
  print, '---------------------------------------------'
  print, 'MESSENGER: No hit!'
  print, '---------------------------------------------'
  arrival_MES=!Values.F_nan
  arrival_speed_MES=!Values.F_nan
endelse

;stop

;same for Venus

d_VEX=elevo_analytic(rdrag, aspectratio, halfwidth, delta_V)

v_VEX=deriv(tars, d_VEX*au)

;check where the heliocentric distance of VEX is less than the ellipse distance
index_hit_VEX=sc_ellipse_coords(direction, a, b, c, pos_vex_xy)
;take first value of these indices = arrival time at VEX within drag time resolution (10 minutes)

if finite(index_hit_VEX) then begin
  arrival_VEX=anytim(tdrag[index_hit_VEX[0]], /ccsds)
  arrival_speed_VEX=v_VEX[index_hit_VEX[0]]
  print, '---------------------------------------------'
  print, 'Arrival at Venus Express [UT]: ', anytim(arrival_VEX, /vms), format='(A, 4x, A17)'
  print, 'Arrival speed at Venus Express [km/s]: ', arrival_speed_VEX, format='(A, 1x, I4)'
endif else begin
  print, '---------------------------------------------'
  print, '---------------------------------------------'
  print, 'Venus Express: No hit!'
  print, '---------------------------------------------'
  arrival_VEX=!Values.F_nan
  arrival_speed_VEX=!Values.F_nan
endelse

;same for STEREO-A
d_A=elevo_analytic(rdrag, aspectratio, halfwidth, delta_A)

v_A=deriv(tars, d_A*au)

;Calculate SC coordinates relative to ellipse center
;index_hit_SC = first value where sc lies within ellipse defined by ELEvoHI
index_hit_A=sc_ellipse_coords(direction, a, b, c, pos_sta_xy)

if finite(index_hit_A) then begin
  arrival_A=anytim(tdrag[index_hit_A[0]], /ccsds)
  arrival_speed_A=v_A[index_hit_A[0]]
  print, '---------------------------------------------'
  print, 'Arrival at STEREO-A [UT]: ', anytim(arrival_A, /vms), format='(A, 4x, A17)'
  print, 'Arrival speed at STEREO-A [km/s]: ', arrival_speed_A, format='(A, 1x, I4)'
  print, '---------------------------------------------'
endif else begin
  print, '---------------------------------------------'
  print, 'STEREO-A: No hit!'
  print, '---------------------------------------------'
  arrival_A=!Values.F_nan
  arrival_speed_A=!Values.F_nan
endelse

;same for B
d_B=elevo_analytic(rdrag, aspectratio, halfwidth, delta_B)

v_B=deriv(tars, d_B*au)

index_hit_B=sc_ellipse_coords(direction, a, b, c, pos_stb_xy)

if finite(index_hit_B) then begin
  arrival_B=anytim(tdrag[index_hit_B[0]], /ccsds)
  arrival_speed_B=v_B[index_hit_B[0]]

  print, '---------------------------------------------'
  print, 'Arrival at STEREO-B [UT]: ', anytim(arrival_B, /vms), format='(A, 4x, A17)'
  print, 'Arrival speed at STEREO-B [km/s]: ', arrival_speed_B, format='(A, 1x, I4)'
  print, '---------------------------------------------'
endif else begin
  print, '---------------------------------------------'
  print, 'STEREO-B: No hit!'
  print, '---------------------------------------------'
  arrival_B=!Values.F_nan
  arrival_speed_B=!Values.F_nan
endelse

;same for Earth (now Wind!)
d_W=elevo_analytic(rdrag, aspectratio, halfwidth, delta_E)


v_W=deriv(tars, d_W*au)

index_hit_W=sc_ellipse_coords(direction, a, b, c, pos_earth_xy)

if finite(index_hit_W) then begin
  arrival_W=anytim(tdrag[index_hit_W[0]], /ccsds)
  arrival_speed_W=v_W[index_hit_W[0]]
  print, '---------------------------------------------'
  print, 'Arrival time at Wind [UT]: ', anytim(arrival_W, /vms), format='(A, 4x, A17)'
  print, 'Arrival speed at Wind [km/s]:', arrival_speed_W, format='(A, 1x, I4)'
  print, '---------------------------------------------'
endif else begin
  print, '---------------------------------------------'
  print, 'Wind: No hit!'
  print, '---------------------------------------------'
  arrival_W=!Values.F_nan
  arrival_speed_W=!Values.F_nan
endelse

;same for Solar Orbiter

if soloAvailable eq 1 then begin

	d_SOLO=elevo_analytic(rdrag, aspectratio, halfwidth, delta_SOLO)

	v_SOLO=deriv(tars, d_SOLO*au)

	;check where the heliocentric distance of VEX is less than the ellipse distance
	index_hit_SOLO=sc_ellipse_coords(direction, a, b, c, pos_solo_xy)
	;take first value of these indices = arrival time at VEX within drag time resolution (10 minutes)


	if finite(index_hit_SOLO) then begin
    arrival_SOLO=anytim(tdrag[index_hit_SOLO[0]], /ccsds)
    arrival_speed_SOLO=v_SOLO[index_hit_SOLO[0]]
	  print, '---------------------------------------------'
	  print, 'Arrival at Solar Orbiter [UT]: ', anytim(arrival_SOLO, /vms), format='(A, 4x, A17)'
	  print, 'Arrival speed at Solar Orbiter [km/s]: ', arrival_speed_SOLO, format='(A, 1x, I4)'
	endif else begin
	  print, '---------------------------------------------'
	  print, '---------------------------------------------'
	  print, 'Solar Orbiter: No hit!'
	  print, '---------------------------------------------'
	  arrival_SOLO=!Values.F_nan
	  arrival_speed_SOLO=!Values.F_nan
	endelse
endif


;same for Parker Solar Probe
if pspAvailable eq 1 then begin
  ;calculates all ellipse distances over longer time frame
	d_PSP=elevo_analytic(rdrag, aspectratio, halfwidth, delta_PSP)

	v_PSP=deriv(tars, d_PSP*au)
	;check where the heliocentric distance of PSP is less than the ellipse distance
	index_hit_PSP=sc_ellipse_coords(direction, a, b, c, pos_psp_xy)
	;take first value of these indices = arrival time at PSP within drag time resolution (10 minutes)


	if finite(index_hit_PSP) then begin
    arrival_PSP=anytim(tdrag[index_hit_PSP[0]], /ccsds)
    arrival_speed_PSP=v_PSP[index_hit_PSP[0]]
	  print, '---------------------------------------------'
	  print, 'Arrival at Parker Solar Probe [UT]: ', anytim(arrival_PSP, /vms), format='(A, 4x, A17)'
	  print, 'Arrival speed at Parker Solar Probe [km/s]: ', arrival_speed_PSP, format='(A, 1x, I4)'
	endif else begin
	  print, '---------------------------------------------'
	  print, '---------------------------------------------'
	  print, 'Parker Solar Probe: No hit!'
	  print, '---------------------------------------------'
	  arrival_PSP=!Values.F_nan
	  arrival_speed_PSP=!Values.F_nan
	endelse
endif

;same for BEPI
if bepiAvailable eq 1 then begin
	d_BEPI=elevo_analytic(rdrag, aspectratio, halfwidth, delta_BEPI)

	v_BEPI=deriv(tars, d_BEPI*au)

	;check where the heliocentric distance of VEX is less than the ellipse distance
	index_hit_BEPI=sc_ellipse_coords(direction, a, b, c, pos_bepi_xy)
	;take first value of these indices = arrival time at VEX within drag time resolution (10 minutes)


	if finite(index_hit_BEPI) then begin
    arrival_BEPI=anytim(tdrag[index_hit_BEPI[0]], /ccsds)
    arrival_speed_BEPI=v_BEPI[index_hit_BEPI[0]]
	  print, '---------------------------------------------'
	  print, 'Arrival at BEPI [UT]: ', anytim(arrival_BEPI, /vms), format='(A, 4x, A17)'
	  print, 'Arrival speed at BEPI [km/s]: ', arrival_speed_BEPI, format='(A, 1x, I4)'
	endif else begin
	  print, '---------------------------------------------'
	  print, '---------------------------------------------'
	  print, 'BEPI: No hit!'
	  print, '---------------------------------------------'
	  arrival_BEPI=!Values.F_nan
	  arrival_speed_BEPI=!Values.F_nan
	endelse
endif

;Add arrival time predictions and CME parameters to plot
if plot then begin
  if runnumber eq 1 then begin ;eq 1
    timesx=0.57
    dragx=0.7

    XYOUTS, !y.crange[1]+2.4, 1, [strmid(t_plot[0],0,11)+'  '+strmid(t_plot[0],12,5)], /data, $
    	charsize=3, charthick=2, alignment=0, font=1, color=colors[0]

    XYOUTS, !y.crange[1]+2.4, 1.15, [strmid(t_plot[1],0,11)+'  '+strmid(t_plot[1],12,5)], /data, $
      charsize=3, charthick=2, alignment=0, font=1, color=colors[1]

    XYOUTS, !y.crange[1]+2.4, 1.30, [strmid(t_plot[2],0,11)+'  '+strmid(t_plot[2],12,5)], /data, $
      charsize=3, charthick=2, alignment=0, font=1, color=colors[2]

    XYOUTS, !y.crange[1]+2.4, 1.45, [strmid(t_plot[3],0,11)+'  '+strmid(t_plot[3],12,5)], /data, $
      charsize=3, charthick=2, alignment=0, font=1, color=colors[3]

    ;last time is at the place with the general ones
    XYOUTS, !y.crange[1]+2.4,1.60, [strmid(t_plot[4],0,11)+'  '+strmid(t_plot[4],12,5)], /data, $
    	charsize=3, charthick=2, alignment=0, font=1, color=colors[4]

    ;DBM parameters upper left
    startx=-1.6
    starty=-1.55

    text4=['Launch time at '+strtrim(string(fix(rinit)),2)+' Rs:']
    XYOUTS, startx,starty, text4, /data, charsize=3, charthick=2,alignment=0, font=1
    text5=[strtrim(strmid(anytim(tinit,/vms),0,17),2)]
    XYOUTS, startx,starty+0.15, text5, /data, charsize=3, charthick=2,alignment=0, font=1

    text1=['Initial speed: '+num2str(vinit,FORMAT='(I4)')+ ' km/s']
    XYOUTS, startx,starty+0.3, text1, /data, charsize=3, charthick=2,alignment=0, font=1
    text2=['Gamma: '+num2str(gammaparam,FORMAT='(F5.2)')]
    XYOUTS, startx,starty+0.45, text2, /data, charsize=3, charthick=2,alignment=0, font=1
    text3=['Background wind: '+num2str(background_wind,FORMAT='(I4)')+ ' km/s']
    XYOUTS, startx,starty+0.6, text3, /data, charsize=3, charthick=2,alignment=0, font=1

    ;general parameters

    dirtext=['Direction '+num2str(earthangle, FORMAT='(F5.1)') + ' deg']
    XYOUTS, startx,starty+0.75, dirtext, /data, charsize=3, charthick=2,alignment=0, font=1

    widthtext=['Half width '+num2str(lambda, FORMAT='(F4.1)')+ ' deg']
    XYOUTS, startx,starty+0.9, widthtext, /data, charsize=3, charthick=2,alignment=0, font=1

    aspecttext=['Aspect ratio '+num2str(aspectratio,FORMAT='(F4.2)')]
    XYOUTS, startx,starty+1.05, aspecttext, /data, charsize=3, charthick=2,alignment=0, font=1

    ;title
    XYOUTS, 0,-1.85, elevo_plot_title, /data, $
    	charsize=4, charthick=3, alignment=0.5, font=1

    device, /close

    loadct, 0, /silent

    set_plot,'X'

  endif
endif
pred = {Wind_time:string(0),  $
	    Wind_speed:float(0.), $
	    STA_time:string(0),	  $
        STA_speed:float(0.),  $
        STB_time:string(0),   $
        STB_speed:float(0.),  $
		MES_time:string(0),   $
		MES_speed:float(0.),   $
		VEX_time:string(0),  $
		VEX_speed:float(0.), $
		SOLO_time:string(0),  $
		SOLO_speed:float(0.),$
		PSP_time:string(0),  $
		PSP_speed:float(0.), $
		BEPI_time:string(0),  $
		BEPI_speed:float(0.)}


pred.Wind_time = arrival_W
pred.Wind_speed= arrival_speed_W
pred.STA_time = arrival_A
pred.STA_speed= arrival_speed_A
pred.STB_time = arrival_B
pred.STB_speed= arrival_speed_B
pred.MES_time = arrival_MES
pred.MES_speed= arrival_speed_MES
pred.VEX_time = arrival_VEX
pred.VEX_speed= arrival_speed_VEX

pred.SOLO_time=!Values.F_nan
pred.SOLO_speed=!Values.F_nan
pred.PSP_time=!Values.F_nan
pred.PSP_speed=!Values.F_nan
pred.BEPI_time=!Values.F_nan
pred.BEPI_speed=!Values.F_nan

if soloAvailable eq 1 then begin
  pred.SOLO_time = arrival_SOLO
  pred.SOLO_speed= arrival_speed_SOLO
endif
if pspAvailable eq 1 then begin
  pred.PSP_time = arrival_PSP
  pred.PSP_speed= arrival_speed_PSP
endif
if bepiAvailable eq 1 then begin
  pred.BEPI_time = arrival_BEPI
  pred.BEPI_speed= arrival_speed_BEPI
endif


elevo_kin = {all_apex_r:dblarr(n_elements(rdrag)),    $
      all_apex_t:strarr(n_elements(rdrag)),   $
      all_apex_lat:fltarr(n_elements(rdrag)), $
        all_apex_lon:fltarr(n_elements(rdrag)), $
        all_apex_s:strarr(n_elements(rdrag)), $
        all_apex_f:fltarr(n_elements(rdrag)),$
        all_apex_w:fltarr(n_elements(rdrag))}

elevo_kin.all_apex_r = rdrag
elevo_kin.all_apex_t = anytim(tdrag, /ccsds)
elevo_kin.all_apex_lat = fltarr(n_elements(rdrag))
elevo_kin.all_apex_lon = fltarr(n_elements(rdrag))+direction
elevo_kin.all_apex_s = make_array(n_elements(rdrag), /string, value=sc)
elevo_kin.all_apex_f = fltarr(n_elements(rdrag))+f
elevo_kin.all_apex_w  = fltarr(n_elements(rdrag))+lambda
;if isa(frame_ind) then stop
END
