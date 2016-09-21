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
; ??/??/15 G. McGarragh: Initial implementation.


; procedure read_baum_lambda
;
; Read and interpolate Bryan Baums's spectral ice cloud optical properties as a
; function of wavelength and effective radius.
;
; INPUT ARGUMENTS:
; filename (string) Full path to the satellite imager model NetCDF file.
; wls      (array)  Wavelengths (um) at which to interpolate optical properties.
; r_e      (array)  Effective radii  at which to interpolate optical properties.
; theta    (array)  Angles at which to interpolate the phase function.
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; Bext (fltarr) Extinction coefficient   (n_elements(wls), n_elements(r_e))
; w    (fltarr) Single scattering albedo (n_elements(wls), n_elements(r_e))
; g    (fltarr) Asymmetry parameter      (n_elements(wls), n_elements(r_e))
; P    (fltarr) Phase function           (n_elements(theta),
;                                         n_elements(wls), n_elements(r_e))
;
; HISTORY:
; ??/??/15 G. McGarragh: Initial implementation.

pro read_baum_lambda, filename, wls, r_e, Bext, w, g, theta, P

   n_wls   = n_elements(wls)
   n_r_e   = n_elements(r_e)
   n_theta = n_elements(theta)

   Bext = fltarr(n_wls, n_r_e)
   w    = fltarr(n_wls, n_r_e)
   g    = fltarr(n_wls, n_r_e)
   P    = fltarr(n_theta, n_wls, n_r_e)

   nc_id = ncdf_open(filename, /nowrite)

   var_id = ncdf_varid(nc_id, 'wavelengths')
   ncdf_varget, nc_id, var_id, wls_baum

   var_id = ncdf_varid(nc_id, 'effective_diameter')
   ncdf_varget, nc_id, var_id, D_e

   r_e_baum = D_e / 2.

   var_id = ncdf_varid(nc_id, 'phase_angles')
   ncdf_varget, nc_id, var_id, phase_angles

   i_wls = interpol(indgen(n_elements(wls_baum)), wls_baum, wls)
   i_r_e = interpol(indgen(n_elements(r_e_baum)), r_e_baum, r_e)

   ; Bext
   baum_read_field_lambda, nc_id, 'extinction_coefficient_over_iwc', $
                           wls_baum, wls, i_wls, r_e_baum, r_e, i_r_e, Bext
   ; w
   baum_read_field_lambda, nc_id, 'single_scattering_albedo', $
                           wls_baum, wls, i_wls, r_e_baum, r_e, i_r_e, w
   ; g
   baum_read_field_lambda, nc_id, 'asymmetry_parameter', $
                           wls_baum, wls, i_wls, r_e_baum, r_e, i_r_e, g

   ; Phs
   var_id = ncdf_varid(nc_id, 'p11_phase_function')
   ncdf_varget, nc_id, var_id, temp

   i_theta = interpol(indgen(n_elements(phase_angles)), phase_angles, theta)

   temp = interpolate(temp, i_wls, i_r_e, i_theta, /grid)

   for i = 0, n_theta - 1 do begin
      P[i,*,*] = temp[*,*,i]
   endfor

   ncdf_close, nc_id
end


pro baum_read_field_lambda, nc_id, name, wls_baum, wls, i_wls, $
                            r_e_baum, r_e, i_r_e, x

   var_id = ncdf_varid(nc_id, name)
   ncdf_varget, nc_id, var_id, temp

   x = interpolate(temp, i_wls, i_r_e, /grid)
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
; r_e      (array)  Effective radii at which to interpolate optical properties.
; theta    (array)  Angles at which to interpolate the phase function.
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; Bext (fltarr) Extinction coefficient   (n_elements(chans), n_elements(r_e))
; w    (fltarr) Single scattering albedo (n_elements(chans), n_elements(r_e))
; g    (fltarr) Asymmetry parameter      (n_elements(chans), n_elements(r_e))
; P    (fltarr) Phase function           (n_elements(theta),
;                                         n_elements(chans), n_elements(r_e))
;
; HISTORY:
; ??/??/15 G. McGarragh: Initial implementation.

pro read_baum_channel, filename, chans, r_e, Bext, w, g, theta, P

   n_chans = n_elements(chans)
   n_r_e   = n_elements(r_e)
   n_theta = n_elements(theta)

   Bext = fltarr(n_chans, n_r_e)
   w    = fltarr(n_chans, n_r_e)
   g    = fltarr(n_chans, n_r_e)
   P    = fltarr(n_theta, n_chans, n_r_e)

   nc_id = ncdf_open(filename, /nowrite)

   var_id = ncdf_varid(nc_id, 'effective_diameter')
   ncdf_varget, nc_id, var_id, D_e

   r_e_baum = D_e / 2.

   var_id = ncdf_varid(nc_id, 'phase_angles')
   ncdf_varget, nc_id, var_id, phase_angles

   i_r_e = interpol(indgen(n_elements(r_e_baum)), r_e_baum, r_e)

   ; Bext
   baum_read_field_channel, nc_id, 'extinction_coefficient_over_iwc', $
                            chans, r_e_baum, r_e, i_r_e, Bext
   ; w
   baum_read_field_channel, nc_id, 'single_scattering_albedo', $
                            chans, r_e_baum, r_e, i_r_e, w
   ; g
   baum_read_field_channel, nc_id, 'asymmetry_parameter', $
                            chans, r_e_baum, r_e, i_r_e, g

   ; Phs
   var_id = ncdf_varid(nc_id, 'p11_phase_function')
   ncdf_varget, nc_id, var_id, temp

   i_theta = interpol(indgen(n_elements(phase_angles)), phase_angles, theta)

   for i = 0, n_chans - 1 do begin
      P[*,i,*] = transpose(interpolate(temp[i,*,*], i_r_e, i_theta, /grid))
   endfor

   ncdf_close, nc_id
end


pro baum_read_field_channel, nc_id, name, chans, r_e_baum, r_e, i_r_e, x

   n_chans = n_elements(chans)

   var_id = ncdf_varid(nc_id, name)
   ncdf_varget, nc_id, var_id, temp

   for i = 0, n_chans - 1 do begin
      x[i,*] = interpolate(temp[chans[i] - 1, *], i_r_e, /grid)
   endfor
end
