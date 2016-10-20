; This library contains procedures for reading Bryan Baums's ice cloud optical
; properties. These properties are described by:
;
; Bryan A. Baum, Ping Yang, Andrew J. Heymsfield, Aaron Bansemer, Benjamin H.
; Cole, Aronne Merrelli, Carl Schmitt, and Chenxi Wang. Ice cloud single-
; scattering property models with the full phase matrix at wavelengths from 0.2
; to 100 mm. Journal of Quantitative Spectroscopy and Radiative Transfer,
; 146:123–139, 2014.
;
; The ice cloud optical property data is available at:
;    https://www.ssec.wisc.edu/ice_models/
;
; History:
; ??/??/15, G. McGarragh: Initial implementation.
; 12/10/16, G. McGarragh: Extrapolate scattering properties for effective radii
;    larger than the maximum Baum effective radius.


; procedure read_baum_lambda
;
; Read and interpolate Bryan Baums's spectral ice cloud optical properties as a
; function of wavelength and effective radius.
;
; INPUT ARGUMENTS:
; filename (string) Full path to the satellite imager model NetCDF file.
; wl       (array)  Wavelengths (um) at which to interpolate optical properties.
; Re       (array)  Effective radii  at which to interpolate optical properties.
; theta    (array)  Angles at which to interpolate the phase function.
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; Bext (fltarr) Extinction coefficient   (n_elements(wl), n_elements(Re))
; w    (fltarr) Single scattering albedo (n_elements(wl), n_elements(Re))
; g    (fltarr) Asymmetry parameter      (n_elements(wl), n_elements(Re))
; P    (fltarr) Phase function           (n_elements(theta),
;                                         n_elements(wl), n_elements(Re))

pro read_baum_lambda, filename, wl, Re, Bext, w, g, theta, P

   n_wl    = n_elements(wl)
   n_Re    = n_elements(Re)
   n_theta = n_elements(theta)

   Bext = fltarr(n_wl, n_Re)
   w    = fltarr(n_wl, n_Re)
   g    = fltarr(n_wl, n_Re)
   P    = fltarr(n_theta, n_wl, n_Re)

   nc_id = ncdf_open(filename, /nowrite)

   var_id = ncdf_varid(nc_id, 'wavelengths')
   ncdf_varget, nc_id, var_id, wl_baum

   var_id = ncdf_varid(nc_id, 'effective_diameter')
   ncdf_varget, nc_id, var_id, De

   Re_baum = De / 2.

   extrap = 0
   n_Re_baum = n_elements(Re_baum)
   if Re_baum[n_Re_baum - 1] lt Re[n_Re - 1] then begin
      extrap = 1
      Re_baum = [Re_baum, Re[n_Re - 1]]
      f_Re_extrap = (Re[n_Re - 1] - Re_baum[n_Re_baum - 2]) / $
                    (Re_baum[n_Re_baum - 1] - Re_baum[n_Re_baum - 2])
   endif

   var_id = ncdf_varid(nc_id, 'phase_angles')
   ncdf_varget, nc_id, var_id, theta_baum

   i_wl = interpol(indgen(n_elements(wl_baum)), wl_baum, wl)
   i_Re = interpol(indgen(n_elements(Re_baum)), Re_baum, Re)

   ; Bext
   baum_read_field_lambda, nc_id, 'extinction_efficiency', $
                           wl_baum, wl, i_wl, Re, i_Re, Bext, extrap, f_Re_extrap
   ; w
   baum_read_field_lambda, nc_id, 'single_scattering_albedo', $
                           wl_baum, wl, i_wl, Re, i_Re, w, extrap, f_Re_extrap
   ; g
   baum_read_field_lambda, nc_id, 'asymmetry_parameter', $
                           wl_baum, wl, i_wl, Re, i_Re, g, extrap, f_Re_extrap

   ; Phs
   var_id = ncdf_varid(nc_id, 'p11_phase_function')
   ncdf_varget, nc_id, var_id, temp

   if extrap then begin
      n_wl_baum    = n_elements(temp[*,0,0])
      n_Re_baum    = n_elements(temp[0,*,0])
      n_theta_baum = n_elements(theta_baum)

      temp2 = fltarr(n_wl_baum, n_Re_baum + 1, n_theta_baum)

      temp2[0:n_wl_baum-1, 0:n_Re_baum-1, 0:n_theta_baum-1] = temp

      for i = 0, n_theta_baum - 1 do begin
         for j = 0, n_wl_baum - 1 do begin
            temp2[j, n_Re_baum, i] = (1. - f_Re_extrap) * temp2[j, n_Re_baum - 2, i] + $
                                           f_Re_extrap  * temp2[j, n_Re_baum - 1, i]
         endfor
      endfor
   endif else begin
      temp2 = temp
   endelse

   i_theta = interpol(indgen(n_elements(theta_baum)), theta_baum, theta)

   temp = interpolate(temp2, i_wl, i_Re, i_theta, /grid)

   for i = 0, n_theta - 1 do begin
      P[i,*,*] = temp[*,*,i]
   endfor

   ncdf_close, nc_id
end


pro baum_read_field_lambda, nc_id, name, wl_baum, wl, i_wl, Re, i_Re, x, $
                            extrap, f_Re_extrap

   var_id = ncdf_varid(nc_id, name)
   ncdf_varget, nc_id, var_id, temp

   n_Re = n_elements(Re)

   if extrap then begin
      n_wl_baum = n_elements(temp[*,0])
      n_Re_baum = n_elements(temp[0,*])

      temp2 = fltarr(n_wl_baum, n_Re_baum + 1)

      temp2[0:n_wl_baum-1, 0:n_Re_baum-1] = temp

      for i = 0, n_wl_baum - 1 do begin
         temp2[i, n_Re_baum] = (1. - f_Re_extrap) * temp2[i, n_Re_baum - 2] + $
                                     f_Re_extrap  * temp2[i, n_Re_baum - 1]
      endfor
   endif else begin
      temp2 = temp
   endelse

   x = interpolate(temp2, i_wl, i_Re, /grid)
end


; procedure read_baum_channel
;
; Read and interpolate Bryan Baums's instrument and channel specific ice cloud
; optical properties known as "satellite imager models" for a set of channels as
; a function of effective radius.
;
; INPUT ARGUMENTS:
; filename (string) Full path to the satellite imager model NetCDF file.
; chans    (array)  Instrument channels for which to obtain optical properties.
; Re       (array)  Effective radii at which to interpolate optical properties.
; theta    (array)  Angles at which to interpolate the phase function.
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; Bext (fltarr) Extinction coefficient   (n_elements(chans), n_elements(Re))
; w    (fltarr) Single scattering albedo (n_elements(chans), n_elements(Re))
; g    (fltarr) Asymmetry parameter      (n_elements(chans), n_elements(Re))
; P    (fltarr) Phase function           (n_elements(theta),
;                                         n_elements(chans), n_elements(Re))
;
; HISTORY:
; ??/??/15 G. McGarragh: Initial implementation.

pro read_baum_channel, filename, chans, Re, Bext, w, g, theta, P

   n_chans = n_elements(chans)
   n_Re    = n_elements(Re)
   n_theta = n_elements(theta)

   Bext = fltarr(n_chans, n_Re)
   w    = fltarr(n_chans, n_Re)
   g    = fltarr(n_chans, n_Re)
   P    = fltarr(n_theta, n_chans, n_Re)

   nc_id = ncdf_open(filename, /nowrite)

   var_id = ncdf_varid(nc_id, 'effective_diameter')
   ncdf_varget, nc_id, var_id, De

   Re_baum = De / 2.

   extrap = 0
   n_Re_baum = n_elements(Re_baum)
   if Re_baum[n_Re_baum - 1] lt Re[n_Re - 1] then begin
      extrap = 1
      Re_baum = [Re_baum, Re[n_Re - 1]]
      f_Re_extrap = (Re[n_Re - 1] - Re_baum[n_Re_baum - 2]) / $
                    (Re_baum[n_Re_baum - 1] - Re_baum[n_Re_baum - 2])
   endif

   var_id = ncdf_varid(nc_id, 'phase_angles')
   ncdf_varget, nc_id, var_id, theta_baum

   i_Re = interpol(indgen(n_elements(Re_baum)), Re_baum, Re)

   ; Bext
   baum_read_field_channel, nc_id, 'extinction_efficiency', $
                            chans, Re, i_Re, Bext, extrap, f_Re_extrap
   ; w
   baum_read_field_channel, nc_id, 'single_scattering_albedo', $
                            chans, Re, i_Re, w, extrap, f_Re_extrap
   ; g
   baum_read_field_channel, nc_id, 'asymmetry_parameter', $
                            chans, Re, i_Re, g, extrap, f_Re_extrap

   ; Phs
   var_id = ncdf_varid(nc_id, 'p11_phase_function')
   ncdf_varget, nc_id, var_id, temp

   if extrap then begin
      n_chans_baum = n_elements(temp[*,0,0])
      n_Re_baum    = n_elements(temp[0,*,0])
      n_theta_baum = n_elements(theta_baum)

      temp2 = fltarr(n_chans_baum, n_Re_baum + 1, n_theta_baum)

      temp2[0:n_chans_baum-1, 0:n_Re_baum-1, 0:n_theta_baum-1] = temp

      for i = 0, n_theta_baum - 1 do begin
         for j = 0, n_chans_baum - 1 do begin
            temp2[j, n_Re_baum, i] = (1. - f_Re_extrap) * temp2[j, n_Re_baum - 2, i] + $
                                           f_Re_extrap  * temp2[j, n_Re_baum - 1, i]
         endfor
      endfor
   endif else begin
      temp2 = temp
   endelse

   i_theta = interpol(indgen(n_elements(theta_baum)), theta_baum, theta)

   for i = 0, n_chans - 1 do begin
      P[*,i,*] = transpose(interpolate(temp2[chans[i] - 1,*,*], i_Re, i_theta, /grid))
   endfor

   ncdf_close, nc_id
end


pro baum_read_field_channel, nc_id, name, chans, Re, i_Re, x, extrap, f_Re_extrap

   n_chans = n_elements(chans)

   var_id = ncdf_varid(nc_id, name)
   ncdf_varget, nc_id, var_id, temp

   n_Re = n_elements(Re)

   if extrap then begin
      n_chans_baum = n_elements(temp[*,0])
      n_Re_baum    = n_elements(temp[0,*])

      temp2 = fltarr(n_chans_baum, n_Re_baum + 1)

      temp2[0:n_chans_baum-1, 0:n_Re_baum-1] = temp

      for i = 0, n_chans_baum - 1 do begin
         temp2[i, n_Re_baum] = (1. - f_Re_extrap) * temp2[i, n_Re_baum - 2] + $
                                     f_Re_extrap  * temp2[i, n_Re_baum - 1]
      endfor
   endif else begin
      temp2 = temp
   endelse

   for i = 0, n_chans - 1 do begin
      x[i,*] = interpolate(temp2[chans[i] - 1, *], i_Re, /grid)
   endfor
end
