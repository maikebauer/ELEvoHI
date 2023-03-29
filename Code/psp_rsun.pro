;get distance of PSP from Sun
FUNCTION psp_rsun, datetime, sc, DISTANCE=distance

heed=get_sunspice_lonlat(datetime,sc,system='hee')
distance=heed[0,*]
r= (6.95508e5 * 648d3 / !dpi ) / distance
IF n_elements(r) GT 1 THEN return,r ELSE return,r[0]

END
