PRO read_input, start, nevents, sc, bflag = bflag

if KEYWORD_SET(bflag) then bflag=bflag else bflag=''

start = STRING(start)

if (sc eq 'SolarOrbiter') then begin

  path = '/nas/helio/data/SolarOrbiter/Events/' + STRMID(start, 0, 4) + '/Tracks/' + start + '/'
  trc_file = FILE_SEARCH(path + '*.csv')

endif

if (sc eq 'A') or (sc eq 'B') then begin

  path = '/nas/helio/data/STEREO/Events/jplot/' + sc + '/' + bflag + '/hi1hi2/' + STRMID(start, 0, 4) + '/Tracks/' + start + '/'
  trc_file = FILE_SEARCH(path + '*' + sc + '*.csv')

endif

FOR i=0, nevents-1 DO BEGIN

  trc_data = READ_CSV(trc_file[i], HEADER=trcHeader, $
  N_TABLE_HEADER=0, TABLE_HEADER=trcTableHeader)

  track_date = trc_data.FIELD1
  elon = trc_data.FIELD2
  sc_fld = trc_data.FIELD3
  elon_stdd = trc_data.FIELD4

  filflag_s = FILE_TEST(path + STRTRIM(STRING(i+1), 2) + '*.sav')
  IF filflag_s THEN FILE_DELETE, path + STRTRIM(STRING(i+1), 2) + '*.sav'

  track = {track_date:track_date, elon:elon, sc:sc_fld, elon_stdd:elon_stdd}

  SAVE, FILENAME = path + 'track_' + STRTRIM(STRING(i+1), 2) + '_' + STRTRIM(start, 2) + '_' + sc + '.sav', track

ENDFOR


END
