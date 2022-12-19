;+
;
; Name:       constrain_elevohi_results
;
; Purpose:    Constraining predictions amde using ELEvoHI with arrival times at one or more spacecraft
;
; Calling sequence: constrain_elevohi_results, date='date', sc_track = 'SC', time_constraints = [t1, t2, t3,...], bflag='bflag'
;
; Parameters (input):
;             date.................'YYYYMMDD' first measurement in HI
;             sc_track.............Spacecraft data used to make time-elongaiton track
;             time_constraints.....array of time constraints in hours for each spacecraft, order must be same as in elevohi_input.txt
;
; Keywords:
;             bflag......set to 'science' or 'beacon' when using STEREO data, leave blank when using other spacecraft

PRO constrain_elevohi_results, date, sc_track, time_constraints, bflag=bflag

path=getenv('ELEvoHI_DIR')


eventdate=date
eventdateSC = eventdate+'_'+sc_track

if KEYWORD_SET(bflag) then eventdateSC = eventdate+'_'+sc_track+'_'+bflag

if not KEYWORD_SET(bflag) then bflag = ''

dir=path+'PredictedEvents/'+eventdateSC+'/'

fnam=dir+'elevohi_input.txt'

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

insitu_targets = []
insitu_times = []

for i=27, n_elements(str[27:-1])-1 do begin
  if not STRMATCH(str[i], '') then begin
    insitu_targets = [insitu_targets, str[i]]
    insitu_times = [insitu_times, str[i+1]]
    i = i+1
  endif
endfor

restore, dir + 'results/' + insitu_targets[0] + '/arrivaltimes.sav'

earth_index = WHERE(insitu_targets EQ 'Earth')
REMOVE, earth_index, insitu_times
REMOVE, earth_index, insitu_targets

sz_1 = SIZE(insitu_targets)
sz_2 = SIZE(arrivaltimes)

arr_times = MAKE_ARRAY(sz_1[1], sz_2[1], /STRING)

for i=0,n_elements(insitu_targets)-1 do begin
    restore, dir + 'results/' + insitu_targets[i] + '/arrivaltimes.sav', /verb
    arr_times[i, *] = arrivaltimes
endfor

sz_3 = SIZE(arr_times)

dt_times = MAKE_ARRAY(sz_1[1], sz_2[1], /FLOAT)

for i=0, sz_3[1]-1 do begin
  for j=0, sz_2[1]-1 do begin
    if finite(arr_times[i, j]) then begin
      dt_times[i, j] = (anytim(arr_times[i, j]) - anytim(insitu_times[i]))/3600.
    endif else begin
      dt_times[i, j] = arr_times[i, j]
    endelse
  endfor
endfor

sz_4 = SIZE(time_constraints)

best_indices = []

for i=0, sz_4[1]-1 do begin
  indices = WHERE(ABS(dt_times[i, *]) LT time_constraints[i])
  if i > 0 then begin
    match, indices, best_indices, suba, subb
    best_indices = [best_indices, indices[suba]]
  endif else begin
    best_indices = [best_indices, indices]
  endelse
endfor

best_indices = best_indices[UNIQ(best_indices, sort(best_indices))]

csvname='/home/mabauer/ELEvoHI/PredictedEvents/' + eventdateSC + '/best_indices.csv'
WRITE_CSV, csvname, best_indices

restore, path+'Code/ASCII_template.sav'
eELEvoHI=read_ascii(dir+'eELEvoHI_results.txt', template=temp)

eELEvoHI_old = eELEvoHI

undefine, eelevohi

new_struc = {PHI:fltarr(n_elements(best_indices)),$
    F:fltarr(n_elements(best_indices)),$
    LAMBDA:fltarr(n_elements(best_indices)),$
    ELONGATION_MIN:fltarr(n_elements(best_indices)),$
    ELONGATION_MAX:fltarr(n_elements(best_indices)),$
    STARTCUT:fltarr(n_elements(best_indices)),$
    ENDCUT:fltarr(n_elements(best_indices)),$
    BG_SW_SPEED:fltarr(n_elements(best_indices)),$
    GAMMA:fltarr(n_elements(best_indices)),$
    TINIT:strarr(n_elements(best_indices)),$
    RINIT:fltarr(n_elements(best_indices)),$
    VINIT:fltarr(n_elements(best_indices)),$
    MEAN_RESIDUAL:fltarr(n_elements(best_indices)),$
    ARRTIME_MES:strarr(n_elements(best_indices)),$
    ARRSPEED_MES:fltarr(n_elements(best_indices)),$
    ARRTIME_VEX:strarr(n_elements(best_indices)),$
    ARRSPEED_VEX:fltarr(n_elements(best_indices)),$
    ARRTIME_EARTH:strarr(n_elements(best_indices)),$
    ARRSPEED_EARTH:fltarr(n_elements(best_indices)),$
    ARRTIME_STA:strarr(n_elements(best_indices)),$
    ARRSPEED_STA:fltarr(n_elements(best_indices)),$
    ARRTIME_STB:strarr(n_elements(best_indices)),$
    ARRSPEED_STB:fltarr(n_elements(best_indices)),$
    ARRTIME_SOLO:strarr(n_elements(best_indices)),$
    ARRSPEED_SOLO:fltarr(n_elements(best_indices)),$
    ARRTIME_PSP:strarr(n_elements(best_indices)),$
    ARRSPEED_PSP:fltarr(n_elements(best_indices)),$
    DT_MES:fltarr(n_elements(best_indices)),$
    DT_VEX:fltarr(n_elements(best_indices)),$
    DT_EARTH:fltarr(n_elements(best_indices)),$
    DT_STA:fltarr(n_elements(best_indices)),$
    DT_STB:fltarr(n_elements(best_indices)),$
    DT_SOLO:fltarr(n_elements(best_indices)),$
    DT_PSP:fltarr(n_elements(best_indices))}

new_struc.PHI = eELEvoHI_old.PHI[best_indices]
new_struc.F = eELEvoHI_old.F[best_indices]
new_struc.LAMBDA = eELEvoHI_old.LAMBDA[best_indices]
new_struc.ELONGATION_MIN = eELEvoHI_old.ELONGATION_MIN[best_indices]
new_struc.ELONGATION_MAX = eELEvoHI_old.ELONGATION_MAX[best_indices]
new_struc.STARTCUT = eELEvoHI_old.STARTCUT[best_indices]
new_struc.ENDCUT = eELEvoHI_old.ENDCUT[best_indices]
new_struc.BG_SW_SPEED = eELEvoHI_old.BG_SW_SPEED[best_indices]
new_struc.GAMMA = eELEvoHI_old.GAMMA[best_indices]
new_struc.TINIT = eELEvoHI_old.TINIT[best_indices]
new_struc.RINIT = eELEvoHI_old.RINIT[best_indices]
new_struc.VINIT = eELEvoHI_old.VINIT[best_indices]
new_struc.MEAN_RESIDUAL = eELEvoHI_old.MEAN_RESIDUAL[best_indices]
new_struc.ARRTIME_MES = eELEvoHI_old.ARRTIME_MES[best_indices]
new_struc.ARRSPEED_MES = eELEvoHI_old.ARRSPEED_MES[best_indices]
new_struc.ARRTIME_VEX = eELEvoHI_old.ARRTIME_VEX[best_indices]
new_struc.ARRSPEED_VEX = eELEvoHI_old.ARRSPEED_VEX[best_indices]
new_struc.ARRTIME_EARTH = eELEvoHI_old.ARRTIME_EARTH[best_indices]
new_struc.ARRSPEED_EARTH = eELEvoHI_old.ARRSPEED_EARTH[best_indices]
new_struc.ARRTIME_STA = eELEvoHI_old.ARRTIME_STA[best_indices]
new_struc.ARRSPEED_STA = eELEvoHI_old.ARRSPEED_STA[best_indices]
new_struc.ARRTIME_STB = eELEvoHI_old.ARRTIME_STB[best_indices]
new_struc.ARRSPEED_STB = eELEvoHI_old.ARRSPEED_STB[best_indices]
new_struc.ARRTIME_SOLO = eELEvoHI_old.ARRTIME_SOLO[best_indices]
new_struc.ARRSPEED_SOLO = eELEvoHI_old.ARRSPEED_SOLO[best_indices]
new_struc.ARRTIME_PSP = eELEvoHI_old.ARRTIME_PSP[best_indices]
new_struc.ARRSPEED_PSP = eELEvoHI_old.ARRSPEED_PSP[best_indices]
new_struc.DT_MES = eELEvoHI_old.DT_MES[best_indices]
new_struc.DT_VEX = eELEvoHI_old.DT_VEX[best_indices]
new_struc.DT_EARTH = eELEvoHI_old.DT_EARTH[best_indices]
new_struc.DT_STA = eELEvoHI_old.DT_STA[best_indices]
new_struc.DT_STB = eELEvoHI_old.DT_STB[best_indices]
new_struc.DT_SOLO = eELEvoHI_old.DT_SOLO[best_indices]
new_struc.DT_PSP = eELEvoHI_old.DT_PSP[best_indices]

eelevohi = new_struc

save, eelevohi, filename = dir + 'eELEvoHI_New.sav'

END
