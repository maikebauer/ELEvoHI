

;read in ht files of SATPLOT measurements

PRO track_statistics_maike, event, sc, make_dir=make_dir, save_file=save_file, no_cor=no_cor, bflag=bflag
;add files to pro at beginning

if KEYWORD_SET(bflag) then bflag=bflag else bflag=''

 IF KEYWORD_SET(make_dir) THEN BEGIN
   cd, '/home/mabauer/hievents/'
   mk_dir, event
   cd
   mk_dir, '/home/mabauer/hievents/' + event
   cd, ''
   print, '...make new event directory...'
 ENDIF

!P.MULTI=0


if (sc eq 'A') or (sc eq 'B') then begin

  path = '/nas/helio/data/STEREO/Events/jplot/' + sc + '/' + bflag + '/hi1hi2/' + STRMID(event, 0, 4) + '/Tracks/' + event + '/'

  infile = FILE_SEARCH(path + '*' + sc + '*.sav')

endif

if (sc eq 'SolarOrbiter') then begin

  path = '/nas/helio/data/SolarOrbiter/Events/' + STRMID(event, 0, 4) + '/Tracks/' + event + '/'

  infile = FILE_SEARCH(path + '*.sav')

endif

if (sc eq 'PSP') then begin

  path = '/nas/helio/data/PSP/Events/' + STRMID(event, 0, 4) + '/Tracks/' + event + '/'

  infile = FILE_SEARCH(path + '*.sav')

endif

;path = '/nas/helio/data/STEREO/Events/Old/' + event + '/jplot/' + bflag + '/'

;infile = files
RESTORE, infile[0]

;read_satplot_ht, infile[0], track

time_axis = dblarr(30)
time_axis[0] = anytim(track.track_date[0])
date=anytim(track.track_date)
sz = SIZE(date)
dt = (date[N_ELEMENTS(date)-1] - date[0])/29.
print, dt
FOR i=1, 29 DO BEGIN

  time_axis[i] = time_axis[i-1]+dt

ENDFOR

ntracks = n_elements(infile)

tracky = FLTARR(ntracks,N_ELEMENTS(time_axis))

tracky[0,*] = interpol(track.elon, date, time_axis)

psym_style = [1,2,3,4,5]
window, 1

utplot, anytim(date, /vms), track.elon;, psym=psym_style[0]



FOR i=1, ntracks-1 DO BEGIN
k=i+1
a=strmid(string(k),7, 1)

;read_satplot_ht, infile[i], track
RESTORE, infile[i]
date=anytim(track.track_date)

tracky[i,*] = interpol(track.elon, date, time_axis)
outplot, anytim(time_axis, /vms), tracky[i,*];, psym=psym_style[i]


ENDFOR

;mean + stddev
y_mean = FLTARR(N_ELEMENTS(tracky[0,*]))
y_stdd1 = FLTARR(N_ELEMENTS(tracky[0,*]))
y_stdd = FLTARR(N_ELEMENTS(tracky[0,*]))

FOR i=0, N_ELEMENTS(tracky[0,*])-1 DO BEGIN

  y_mean[i] = MEAN(tracky[*,i], /NaN)
  y_stdd1[i] = STDDEV(tracky[*,i], /NaN)


ENDFOR
print, time_axis
print, y_mean
loadct, 4
outplot, anytim(time_axis, /vms), y_mean, psym=1, symsize=4, thick=3, color=180

loadct, 0

if (sc eq 'A') or (sc eq 'B') then begin
if keyword_set(no_cor) then begin

hi1 = 0

lasthi = where(y_mean ge 20)
hi1 = lasthi[0]

window, 2
utplot, anytim(time_axis, /vms), y_mean, psym=1

print, 'Standarddeviation:', y_stdd1
print, '*******************'
print, 'Last element of HI1:'
print, hi1
print, '*******************'

err_hi1  =  mean(y_stdd1[0:hi1])

size_ystdd1 = SIZE(y_stdd1)

IF hi1 EQ size_ystdd1[1]-1 THEN BEGIN

  hi2 = hi1
  err_hi2 = mean(y_stdd1[0:hi2])


ENDIF ELSE BEGIN

  hi2 = hi1+1
  err_hi2  = mean(y_stdd1[hi2:*])


ENDELSE

print, 'error in degree:'
print, '****************'
print, 'error hi1:', err_hi1
print, 'error hi2:', err_hi2

y_stdd[0:hi1] = err_hi1
y_stdd[hi2:*] = err_hi2

endif else begin

cor2 = 0
hi1 = 0

lasthi = where(y_mean ge 20)
hi1 = lasthi[0]

window, 2
utplot, anytim(time_axis, /vms), y_mean, psym=1
print, 'Standarddeviation:', y_stdd1
print, 'Last element of COR2:'
read, cor2
print, '*******************'
print, 'Last element of HI1:'
print, hi1
print, '*******************'

err_cor2 = mean(y_stdd1[0:cor2])
err_hi1  =  mean(y_stdd1[cor2+1:hi1])
err_hi2  = mean(y_stdd1[hi2:*])

print, 'error in degree:'
print, '****************'
print, 'error cor2:', err_cor2
print, 'error hi1:', err_hi1
print, 'error hi2:', err_hi2

y_stdd[0:cor2] = err_cor2
y_stdd[cor2+1:hi1] = err_hi1
y_stdd[hi2:*] = err_hi2

endelse
endif
print, 'Mean Standarddeviation:', mean(y_stdd)

;define structure for saving


track = {track_date:STRARR(N_ELEMENTS(time_axis)),    $
         elon:FLTARR(N_ELEMENTS(y_mean)),    $
         elon_stdd:FLTARR(N_ELEMENTS(y_stdd1)),    $
         sc:STRING('')}

track.track_date  =  anytim(time_axis, /vms)
track.elon  =  y_mean
track.elon_stdd  =  y_stdd1
track.sc  =  sc


IF KEYWORD_SET(save_file) THEN BEGIN

  if (sc eq 'A') or (sc eq 'B') then begin

    test = FILE_EXIST('/nas/helio/data/STEREO/HItracks/mabauer/')

    IF NOT test THEN FILE_MKDIR, '/nas/helio/data/STEREO/HItracks/mabauer/'

    file= '/nas/helio/data/STEREO/HItracks/mabauer/' + event + '_' + sc + '_' + bflag + '.sav'

  endif

  if (sc eq 'SolarOrbiter') then begin

    test = FILE_EXIST('/nas/helio/data/SolarOrbiter/HItracks/mabauer/')

    IF NOT test THEN FILE_MKDIR, '/nas/helio/data/SolarOrbiter/HItracks/mabauer/'

    file= '/nas/helio/data/SolarOrbiter/HItracks/mabauer/' + event + '.sav'

  endif

  if (sc eq 'PSP') then begin

    test = FILE_EXIST('/nas/helio/data/PSP/HItracks/mabauer/')

    IF NOT test THEN FILE_MKDIR, '/nas/helio/data/PSP/HItracks/mabauer/'

    file= '/nas/helio/data/PSP/HItracks/mabauer/' + event + '.sav'

  endif


  print, file
  save, track, filename=file
  help, track
  print, '...save file...'
ENDIF

END
