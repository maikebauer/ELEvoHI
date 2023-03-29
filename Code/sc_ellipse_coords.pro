;converts spacecraft coordinates into coordinates relative to ellipse center
;calculates if spacecraft is withi n ellipse defined by CME parameters
FUNCTION sc_ellipse_coords, direction, a, b, c, sc_pos


angle = 180-direction

;get x and y angle of ellipse
;direction relative to Earth
xangle=sin((direction)*!dtor)
yangle=cos((direction)*!dtor)

;calculate ellipse center
ellipse_center = fltarr(2, n_elements(c))
ellipse_center[1, *] = xangle*c
ellipse_center[0, *] = yangle*c

;calculate position of SC relative to ellipse center
sin_angle = sin(!pi-angle*!dtor)
cos_angle = cos(!pi-angle*!dtor)

xc = sc_pos[0, *] - ellipse_center[0, *]
yc = sc_pos[1, *] - ellipse_center[1, *]

;calculate points along semi-major and -minor axis
;Formula is just calculation for new points after rotation
xct = xc * cos_angle - yc * sin_angle
yct = xc * sin_angle + yc * cos_angle

;ellipse equaiton
;if LE 1 -> within ellipse defined by ELEvoHI
rad_cc = (xct^2/a^2) + (yct^2/b^2)

index_inside_ellipse = MIN(WHERE(rad_cc LE 1))

;return smallest index where hit criterium is met
if index_inside_ellipse ne -1 then begin
	index_hit = index_inside_ellipse
endif else begin
	index_hit = !Values.F_NAN
endelse

RETURN, index_hit
END