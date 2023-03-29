;Calculate direciton of CME apex relative to Eartha t one point in time
;Takes into account SC position and pointing
FUNCTION get_orientation, phi, sc_pos, pos_E, crval

;calculate direction from Earth at tinit

if sc_pos[1] ge 0 and crval ge 0 then begin
  sep=abs(pos_E[1]-sc_pos[1])/!dtor
  delta_Earth=sep+phi
endif

if sc_pos[1] ge 0 and crval lt 0 then begin
  sep=abs(pos_E[1]-sc_pos[1])/!dtor
  delta_Earth=sep-phi
endif

if sc_pos[1] lt 0 and crval ge 0 then begin
  sep=abs(pos_E[1]-sc_pos[1])/!dtor
  delta_Earth=-(sep-phi)
endif

if sc_pos[1] lt 0 and crval lt 0 then begin
  sep=abs(pos_E[1]-sc_pos[1])/!dtor
  delta_Earth=-(sep+phi)
endif

RETURN, delta_Earth

END