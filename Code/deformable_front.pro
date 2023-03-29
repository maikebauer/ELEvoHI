;+
;
; Name:       deformable_front
;
; Purpose:    makes the front of the CME deformable by adjusting the kinematics to the ambient solar wind conditions
;              for each ensemble member a .sav-file is stored with all the needed information
;              at the end of elevohi, all these .sav-files are combined to one final .sav file ('frontDataAll.sav')
;
; Calling sequence: deformable_front, lambda, f, phi, kappa, tinit, fitend, swspeed, drag_parameter, defFrontStartTime, speedEndCut, sc, bgswdata, bgswTimeNum, runnumber, resdir
;
; Parameters (input):
;              bgsw: which ambient solar wind is used ('HUX' or 'HUXt')
;              lambda: angular half width within ecliptic in degree
;              f: inverse ellipse aspect ratio (b/a)
;              phi: apex direction from the HI observer in degree
;              kappa: angular half width with respect to the latitude
;              tinit: initial time of the CME
;              fitend: distance of the last fit from DBM-fitting
;              swspeed: solar wind speed from DBM-fitting
;              drag_parameter: drag parameter from DBM-fitting
;              defFrontStartTime: start time for the deformable front
;              speedEndCut: speed of the solar wind at fitend
;              sc: HI observer ['A' or 'B']
;              bgswdata: ambient solar wind speed from WSA/HUX model combination
;              bgswTimeNum: time for which the ambient solar wind speed was created
;              runnumber: number of run in the ensemble mode
;              resdir: directory for the .sav-files
;
; History:    2021/03: created (Juergen Hinterreiter)
;
; Authors:    Tanja Amerstorfer & Christian Moestl & Juergen Hinterreiter
;             Space Research Institute, Austrian Academy of Sciences
;              Graz, Austria
; -

pro deformable_front, bgsw, lambda, f, phi, kappa, tinit, fitend, swspeed, drag_parameter, defFrontStartTime, speedEndCut, sc, bgswdata, bgswTimeNum, runnumber, resdir
    runtimeStart = systime(/seconds)
    au=double(149597870.)
    r_sun=double(695700.)


    numberOfEllPoints = 101
    timeResolution = 15 ; in minutes

    maxRRun = 300
    if strupcase(bgsw) eq 'HUXT' then begin
        timeResolution = 0
        while timeResolution lt 15 do begin
            timeResolution = timeResolution + bgswdata.dt/60
        endwhile
        maxRRun = fix(max(bgswData.r))
    endif

    if strupcase(bgsw) eq 'EUHFORIA' then maxRRun = fix(max(bgswData.r))

    print, timeResolution
    runDays = 5
    numberTimeArray = runDays*24*60/timeResolution


    if finite(tinit) ne 0 and tinit ne 0 then begin
        
        ;get SC positions
        if (sc eq 'A') or (sc eq 'B') then begin
            pos_E=get_sunspice_lonlat(tinit, 'Earth', system='HEE')
            pos_A=get_sunspice_lonlat(tinit, 'Ahead', system='HEE')
            pos_B=get_sunspice_lonlat(tinit, 'Behind', system='HEE')
        endif

        if (sc eq 'Solar_Orbiter') then begin
            pos_E=get_sunspice_lonlat(tinit, 'Earth', system='HEE')
            pos_SolO=get_sunspice_lonlat(tinit, 'Solar_Orbiter', system='HEE')
        endif

        if (sc eq 'PSP') then begin
            pos_E=get_sunspice_lonlat(tinit, 'Earth', system='HEE')
            pos_PSP=get_sunspice_lonlat(tinit, 'PSP', system='HEE')
        endif

        print, phi

        ;calculate direction from Earth
        if sc eq 'A' then begin
          dir_E = get_orientation(phi, pos_A, pos_E, crval)
        endif

        if sc eq 'B' then begin
          dir_E = get_orientation(phi, pos_B, pos_E, crval)
        endif

        if sc eq 'Solar_Orbiter' then begin
          dir_E = get_orientation(phi, pos_SolO, pos_E, crval)
        endif

        if sc eq 'PSP' then begin
          dir_E = get_orientation(phi, pos_PSP, pos_E, crval)
        endif

        print, 'runnumber: ', runnumber

        distAU = double(fitend) ; [AU]
        distRSun = double(fitend*au/r_sun) ; [R_sun]
        distKM = double(fitend*au) ; [km]
        areaFact = double(1);(0.25*au)*(0.25*au)
        earthangle = -1 * dir_e

        theta=atan(f^2*tan(lambda*!dtor))
        omega=sqrt(cos(theta)^2*(f^2-1)+1)
        ;if this factor is set to other than 1 one can make very wide ellipses around the Sun
        ;if necessary
        factor=1  ; another possible free parameter
        b=distRSun*omega*sin(lambda*!dtor)/(cos(lambda*!dtor-theta)+omega*sin(lambda*!dtor))*factor
        a=b/f
        c=distRSun-b

        xangle=sin((earthangle)*!dtor)
        yangle=cos((earthangle)*!dtor)
        ellipse_center=[xangle,yangle]*c

        pos_ang = 90 - earthangle; angle with repect to Earth
        xc = ellipse_center[0]; x center of ellipse
        yc = ellipse_center[1]; y center of ellipse
        rmin = a;
        rmax = b;


        kap = kappa*!dtor
        cLat = distRSun*tan(kap)

        csArea = a*cLat*!pi*r_sun*r_sun ; [km^2]


        rhoSW = double(get_bgsw_density(fitend, swspeed, mass_dens=1)) ; [g/km^3]
        dp = abs(double(drag_parameter)) ; [km^-1]
        dp = double(drag_parameter) ; [km^-1]
        mass = csArea*rhoSW/abs(dp) ; [g]

;        print, 'EC: ', eC
        print, 'Model: ', bgsw
        print, 'SWspeed [km/s]: ', swspeed
        print, 'SPEndCut [km/s]: ', speedEndCut
        print, 'rhoSW [g/km^3]: ', rhoSW
        print, 'dp [km^-1]:', dp
        print, 'csArea [km^2]: ', csArea
        print, 'mass [g]: ', mass

        dragNew = csArea*rhoSW/mass
;        print, 'drag new: '
;        print, dragNew

        ;to draw parts of ellipse phi should go from 0 to 135 and then from 225 to 360
        ; min(phiell) has to be 260 max(phiell) has to be 460 (because 360 is direction)
        phiEll=findgen(numberOfEllPoints)*200/(numberofEllPoints-1)+260

        phiEll=phiEll*!dpi/180.

        ;original version
        ;phi = 2*!pi*(findgen(npoints)/(npoints-1))       ;Divide circle into Npoints
        ang = pos_ang/!RADEG                            ;Position angle in radians
        cosang = cos(ang)
        sinang = sin(ang)

        x =  rmax*cos(phiEll)              ;Parameterized equation of ellipse
        y =  rmin*sin(phiEll)

        xprime = xc + x*cosang - y*sinang      ;Rotate to desired position angle
        yprime = yc + x*sinang + y*cosang


        R_ellipse = sqrt(xprime*xprime + yprime*yprime)
        lon_ell = atan(yprime, xprime)

        lonellnew = (lon_ell - lon_ell[n_elements(lon_ell)/2])*180/!dpi
        lonNew = lonellnew + earthangle
        ;print, lonNew

        vbefore = fltarr(n_elements(r_ellipse))
        dragRuns = dblarr(n_elements(r_ellipse))
        densRuns = dblarr(n_elements(r_ellipse))
        swRuns = dblarr(n_elements(r_ellipse))

        vbefore[*] = speedEndcut
		for i = 0, n_elements(r_ellipse)-1 do begin
            sval=(elevo_analytic(r_ellipse[i]/au*r_sun, 1/f, lambda, lonellnew[i], speed=speedEndcut))[1]
            vbefore[i] = sval
        endfor

        rbefore = r_ellipse
        tini = defFrontStartTime
        tinitnum = tini

        tdrag = dblarr(numberTimeArray)
        tdrag[0] = tinitnum

        timeLen = n_elements(tdrag)
        ellLen = n_elements(r_ellipse)

        frontArr = dblarr(timelen, ellLen)
        vArr = dblarr(timelen, ellLen)
        dragArr = dblarr(timelen, ellLen)
        densArr = dblarr(timelen, ellLen)
        swArr = dblarr(timelen, ellLen)

        print, 'tini: ', tini
        for j = 0, n_elements(tdrag)-1 do begin
            tdrag[j] = (j+1)*timeResolution*60. + tini

            ; cLat is the latitduinal extent of the CME front
            ; taken is the radius of the apex to calculate the latitudinal extent
            cLat = rbefore[n_elements(rbefore)/2]*tan(kap)
            areaRun = cLat*a*!pi*r_sun*r_sun
            if max(rbefore) lt maxRRun then begin
                for i=0, n_elements(r_ellipse)-1 do begin
                    if strupcase(bgsw) eq 'HUX' then begin
                        sw_speed = get_bgsw_speed_hux(bgswdata=bgswdata, bgswTimeNum, tinitnum, lonNew[i], rbefore[i], i)
                        sw_density = double(get_bgsw_density(rbefore[i]/au*r_sun, sw_speed, mass_dens=1))
                    endif
                    if strupcase(bgsw) eq 'HUXT' then begin
                        sw_speed = get_bgsw_speed_huxt(bgswdata=bgswdata, tdrag[j], lonNew[i], rbefore[i])
                        sw_density = double(get_bgsw_density(rbefore[i]/au*r_sun, sw_speed, mass_dens=1))
                    endif
                    if strupcase(bgsw) eq 'EUHFORIA' then begin
                        swParam = get_bgsw_speed_euhforia(bgswdata=bgswdata, tdrag[j], lonNew[i], rbefore[i])
                        sw_speed = swParam[0]
                        sw_density = swParam[1]
                    endif

                    accsign = 1
                    gammaparam = areaRun*sw_density/mass
                    background_wind = sw_speed

                    if vbefore[i] lt background_wind then accsign = -1

                    dragRuns[i] = gammaparam*accsign
                    densRuns[i] = sw_density
                    swRuns[i] = background_wind

                    ; calculate the new distance and speed at each point
                    rnew=(accsign/(gammaparam))*alog(1+(accsign*(gammaparam)*((vbefore[i]-background_wind)*(tdrag[j]-tinitnum))))+background_wind*(tdrag[j]-tinitnum)+rbefore[i]*r_sun
                    vnew=(vbefore[i]-background_wind)/(1+(accsign*(gammaparam)*((vbefore[i]-background_wind)*(tdrag[j]-tinitnum))))+background_wind
                    if finite(rnew) eq 0 then begin
                        print, '!!!!!!!!!!!!!!!'
                        print, '!!!!!!!!!!!!!!!'
                        print, '!! nan value !!'
                        print, '!!!!!!!!!!!!!!!'
                        print, '!!!!!!!!!!!!!!!'

                        print, 'i: ', i
                        print, 'gamma param: ', gammaparam
                        print, 'sw dens: ', sw_density
                        print, 'vbefore: ', vbefore[i]
                        print, 'bgsw: ', background_wind
                        print, 'rbefore, ', rbefore[i]

                        print, '!!!!!!!!'
                        print, 'gamma original: ', drag_parameter
                        print, 'dens original: ', rhoSW
                        stop
                     endif
                     rnew = rnew/r_sun

                     rbefore[i] = rnew
                     vbefore[i] = vnew
                endfor

                frontArr[j, *] = rbefore
                vArr[j, *] = vbefore
                dragArr[j, *] = dragRuns
                densArr[j, *] = densRuns
                swArr[j, *] = swRuns
            endif else begin
                j = n_elements(tdrag)-1
            endelse
            tinitnum = tdrag[j]

        endfor

        timearr = tdrag

        noEarthHit = 0
        indMinLon = where(abs(lonnew) eq min(abs(lonnew)))
        if min(abs(lonnew)) gt 1 then noEarthHit = 1
        radsEarth = frontarr[*, indMinLon]
        distEarth = pos_E[0]/r_sun - 1.5e6/r_sun ; correct for L1
        print, distEarth

        indMinRad = where(abs(radsEarth - distEarth) eq min(abs(radsEarth - distEarth)))

        distEarth = radsEarth[indMinRad]
        vEarth = vArr[indMinRad, indMinLon]
        arrtimeEarth = anytim(timearr[indMinRad], /ccsds)
        if 0 eq 1 then begin
            print, 'distEarth: ', radsEarth[indMinRad]
            print, 'Earth Arr: ', arrtimeEarth
            print, 'min Lon: ', lonnew[indMinLon]
            print, 'DistRSun: ', distRSun
            print, 'area: ', csArea
            print, 'mass: ', mass
        endif

        indEarthDirection = indMinLon
        longitude = lonnew
        timearr = anytim(timearr, /ccsds)

        if noEarthHit eq 1 then begin
            indEarthDirection = -1
            arrtimeEarth = -1
            vEarth = -1
        endif

        runArea = csArea
        runMass = mass
        runRho = rhoSW
        distMass = distRSun
        runDP = dp
        filenameFront = 'frontData_'+string(runnumber, format='(I003)')+'.sav'
        save, timearr, frontarr, vArr, dragArr, densArr, swArr, longitude, indEarthDirection, distEarth, arrtimeEarth, vEarth, phi, lambda, f, runarea, runmass, distmass, runRho, runDP, filename = resdir + filenameFront

        if runnumber eq 1 then begin
            print, 'Save front Data'
            filenameFront = 'firstRun_frontData.sav'
            save, timearr, frontarr, vArr, dragArr, densArr, swArr, longitude, indEarthDirection, distEarth, arrtimeEarth, vEarth, phi, lambda, f, runarea, runmass, distmass, runRho, runDP, filename = resdir + filenameFront
        endif

        print, 'single run neededs ', (systime(/seconds)-runtimeStart), ' seconds', format='(A20, F5.1, A9)'
     endif
end
