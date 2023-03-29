;restore, '/home/tamerstorfer/ElEvoHI_events/20101103_newer/all_results_mass.sav', /verb
;restore, 'ELEvoHI/PredictedEvents/20200415_A_final/results/SOLO/arrivaltimes.sav', /verb
restore, '/home/mabauer/ELEvoHI/PredictedEvents/20220311_A_science/results/SOLO/arrivaltimes.sav', /verb
solo_arrtimes = arrivaltimes
restore, '/home/mabauer/ELEvoHI/PredictedEvents/20220311_A_science/results/Earth/arrivaltimes.sav', /verb
;restore, 'ELEvoHI/PredictedEvents/20200415_A_final/results/Earth/arrivaltimes.sav', /verb
; restore, 'ELEvoHI/PredictedEvents/20200415_A_track1/results/Earth/sw.sav', /verb
; restore, 'ELEvoHI/PredictedEvents/20200415_A_track1/results/Earth/gamma.sav', /verb
; restore, 'ELEvoHI/PredictedEvents/20200415_A_track1/results/Earth/transittimes.sav', /verb
; restore, 'ELEvoHI/PredictedEvents/20200415_A_track1/results/Earth/vinit.sav', /verb
; restore, 'ELEvoHI/PredictedEvents/20200415_A_track1/results/Earth/rinit.sav', /verb
; restore, 'ELEvoHI/PredictedEvents/20200415_A_track1/results/Earth/prediction.sav', /verb
; restore, 'ELEvoHI/PredictedEvents/20200415_A_track1/results/Earth/labels.sav', /verb
; restore, 'ELEvoHI/PredictedEvents/20200415_A_track1/results/Earth/plottimes.sav', /verb
; restore, 'ELEvoHI/PredictedEvents/20200415_A_track1/results/Earth/residuals.sav', /verb



;SOLO = anytim('2020-04-19T05:07')
;Earth = anytim('2020-04-20T01:34')
SOLO = anytim('2022-03-11T19:52')
Earth = anytim('2022-03-13T10:11')

dt_earth = dblarr(n_elements(arrivaltimes))
dt_solo = dblarr(n_elements(arrivaltimes))

for i=0, n_elements(arrivaltimes)-1 do begin
    if finite(arrivaltimes[i]) then begin
        dt_earth[i] = (anytim(arrivaltimes[i]) - Earth)/3600.
        print, anytim(arrivaltimes[i])
    endif else begin
        dt_earth[i] = arrivaltimes[i]
    endelse
endfor

for i=0, n_elements(solo_arrtimes)-1 do begin
    if finite(solo_arrtimes[i]) then begin
        dt_solo[i] = (anytim(solo_arrtimes[i]) - SOLO)/3600.
    endif else begin
        dt_solo[i] = solo_arrtimes[i]
    endelse
endfor


best = where(abs(dt_solo) lt 4)


eventdateSC = '20220311_A_science'

path=getenv('ELEvoHI_DIR')
dir=path+'PredictedEvents/'+eventdateSC+'/'

restore, path+'Code/ASCII_template.sav'
eELEvoHI=read_ascii(dir+'eELEvoHI_results.txt', template=temp)

eELEvoHI_old = eELEvoHI

undefine, eelevohi

new_struc = {PHI:fltarr(n_elements(best)),$
    F:fltarr(n_elements(best)),$
    LAMBDA:fltarr(n_elements(best)),$
    ELONGATION_MIN:fltarr(n_elements(best)),$
    ELONGATION_MAX:fltarr(n_elements(best)),$
    STARTCUT:fltarr(n_elements(best)),$
    ENDCUT:fltarr(n_elements(best)),$
    BG_SW_SPEED:fltarr(n_elements(best)),$
    GAMMA:fltarr(n_elements(best)),$
    TINIT:strarr(n_elements(best)),$
    RINIT:fltarr(n_elements(best)),$
    VINIT:fltarr(n_elements(best)),$
    MEAN_RESIDUAL:fltarr(n_elements(best)),$
    ARRTIME_MES:strarr(n_elements(best)),$
    ARRSPEED_MES:fltarr(n_elements(best)),$
    ARRTIME_VEX:strarr(n_elements(best)),$
    ARRSPEED_VEX:fltarr(n_elements(best)),$
    ARRTIME_EARTH:strarr(n_elements(best)),$
    ARRSPEED_EARTH:fltarr(n_elements(best)),$
    ARRTIME_STA:strarr(n_elements(best)),$
    ARRSPEED_STA:fltarr(n_elements(best)),$
    ARRTIME_STB:strarr(n_elements(best)),$
    ARRSPEED_STB:fltarr(n_elements(best)),$
    ARRTIME_SOLO:strarr(n_elements(best)),$
    ARRSPEED_SOLO:fltarr(n_elements(best)),$
    ARRTIME_PSP:strarr(n_elements(best)),$
    ARRSPEED_PSP:fltarr(n_elements(best)),$
    DT_MES:fltarr(n_elements(best)),$
    DT_VEX:fltarr(n_elements(best)),$
    DT_EARTH:fltarr(n_elements(best)),$
    DT_STA:fltarr(n_elements(best)),$
    DT_STB:fltarr(n_elements(best)),$
    DT_SOLO:fltarr(n_elements(best)),$
    DT_PSP:fltarr(n_elements(best))}

new_struc.PHI = eELEvoHI_old.(0)[best]
new_struc.F = eELEvoHI_old.(1)[best]
new_struc.LAMBDA = eELEvoHI_old.(2)[best]
new_struc.ELONGATION_MIN = eELEvoHI_old.(3)[best]
new_struc.ELONGATION_MAX = eELEvoHI_old.(4)[best]
new_struc.STARTCUT = eELEvoHI_old.(5)[best]
new_struc.ENDCUT = eELEvoHI_old.(6)[best]
new_struc.BG_SW_SPEED = eELEvoHI_old.(7)[best]
new_struc.GAMMA = eELEvoHI_old.(8)[best]
new_struc.TINIT = eELEvoHI_old.(9)[best]
new_struc.RINIT = eELEvoHI_old.(10)[best]
new_struc.VINIT = eELEvoHI_old.(11)[best]
new_struc.MEAN_RESIDUAL = eELEvoHI_old.(12)[best]
new_struc.ARRTIME_MES = eELEvoHI_old.(13)[best]
new_struc.ARRSPEED_MES = eELEvoHI_old.(14)[best]
new_struc.ARRTIME_VEX = eELEvoHI_old.(15)[best]
new_struc.ARRSPEED_VEX = eELEvoHI_old.(16)[best]
new_struc.ARRTIME_EARTH = eELEvoHI_old.(17)[best]
new_struc.ARRSPEED_EARTH = eELEvoHI_old.(18)[best]
new_struc.ARRTIME_STA = eELEvoHI_old.(19)[best]
new_struc.ARRSPEED_STA = eELEvoHI_old.(20)[best]
new_struc.ARRTIME_STB = eELEvoHI_old.(21)[best]
new_struc.ARRSPEED_STB = eELEvoHI_old.(22)[best]
new_struc.ARRTIME_SOLO = eELEvoHI_old.(23)[best]
new_struc.ARRSPEED_SOLO = eELEvoHI_old.(24)[best]
new_struc.ARRTIME_PSP = eELEvoHI_old.(25)[best]
new_struc.ARRSPEED_PSP = eELEvoHI_old.(26)[best]
new_struc.DT_MES = eELEvoHI_old.(27)[best]
new_struc.DT_VEX = eELEvoHI_old.(28)[best]
new_struc.DT_EARTH = eELEvoHI_old.(29)[best]
new_struc.DT_STA = eELEvoHI_old.(30)[best]
new_struc.DT_STB = eELEvoHI_old.(31)[best]
new_struc.DT_SOLO = eELEvoHI_old.(32)[best]
new_struc.DT_PSP = eELEvoHI_old.(33)[best]


eelevohi = new_struc

save, eelevohi, filename = '/home/mabauer/ELEvoHI/PredictedEvents/' + eventdateSC + '/eELEvoHI_New.sav'

fname='/home/mabauer/ELEvoHI/PredictedEvents/' + eventdateSC + '/best_indices.csv'
WRITE_CSV, fname, best

print, 'done'

end
