; This library contains procedures for reading Anthony Baran's ice cloud optical
; properties. These properties are described by:
;
; A. J. Baran and S. Havemann. The dependence of retrieved cirrus ice-crystal
; effective dimension on assumed ice-crystal geometry and size-distribution
; function at solar wavelengths. Quarterly Journal of the Royal Meteorological
; Society, 130(601):21532167, 2004. doi:10.1256/qj.03.154.
;
; History:
; ??/??/15 G. McGarragh: Initial implementation.


; procedure read_baran
;
; Read and interpolate Anthony Baran's ice cloud optical properties for a set of
; wavelengths and effective radii.
;
; INPUT ARGUMENTS:
; path  (string) Path to the base directory containing the optical properties.
; wl    (array)  Wavelengths for which to obtain optical properties.
; r_e   (array)  Effective radii for which to obtain optical properties.
; theta (array)  Angles to which to interpolate the phase function output.
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; Bext (fltarr) Extinction coefficient   (n_elements(wl), n_elements(r_e))
; w    (fltarr) Single scattering albedo (n_elements(wl), n_elements(r_e))
; g    (fltarr) Asymmetry parameter      (n_elements(wl), n_elements(r_e))
; P    (fltarr) Phase function           (n_elements(theta),
;                                         n_elements(wl), n_elements(r_e))
;
; HISTORY:
; ??/??/15 G. McGarragh: Initial implementation.

pro read_baran, path, wl, r_e, Bext, w, g, theta, P

   n_wl    = n_elements(wl)
   n_r_e   = n_elements(r_e)
   n_theta = n_elements(theta)

   Bext    = fltarr(n_wl, n_r_e)
   w       = fltarr(n_wl, n_r_e)
   g       = fltarr(n_wl, n_r_e)
   P       = fltarr(n_theta, n_wl, n_r_e)

   theta1  = fltarr(1801)
   P1      = fltarr(1801)

   filename = path + '/' + 'solar' + '/' + $
              'solar_optical_phase_aggregate.dat'

   openr, lun, filename, /get_lun

   for i=0,n_r_e-1 do begin
      for j=0,n_wl-1 do begin
         if wl[j] lt 4. then begin
            point_lun, lun, 0
            baran_read_sw_interp, lun,  wl[j], 2.*r_e[i], $
                                  Bext1, w1, g1, theta1, P1
         endif else begin
            baran_read_lw_interp, path, wl[j],    r_e[i], $
                                  Bext1, w1,     theta1, P1
         endelse
         Bext[j,i] = Bext1
         w[j,i]    = w1
         g[j,i]    = g1

         P[*,j,i] = interpol(P1, theta1, theta)
      endfor
   endfor

   free_lun, lun
end


pro baran_combine_props, a, Bext1, w1, g1, P1, Bext2, w2, g2, P2, $
                         Bext, w, g, P

   Bext = (1. - a) * Bext1 + a * Bext2
   w    = (1. - a) * w1    + a * w2
   g    = (1. - a) * g1    + a * g2

   n = size(P, /dimensions)

   for i=0,n[0]-1 do begin
      P[i] = (1. - a) * P1[i] + a * P2[i]
   end
end


pro baran_read_sw_interp, lun, wl, De, Bext, w, g, theta, P

   wls = [0.2250000, 0.2750000, 0.3500000, 0.4250000, 0.4750000, $
          0.5500000, 0.6500000, 0.7500000, 0.8500000, 0.9500000, $
          1.050000,  1.150000,  1.250000,  1.350000,  1.450000,  $
          1.575000,  1.725000,  1.850000,  1.925000,  1.975000,  $
          2.050000,  2.150000,  2.250000,  2.350000,  2.500000,  $
          2.650000,  2.750000,  2.812000,  2.838000,  2.862000,  $
          2.888000,  2.912000,  2.938000,  2.963000,  2.987000,  $
          3.013000,  3.037000,  3.062000,  3.088000,  3.112000,  $
          3.138000,  3.162000,  3.188000,  3.225000,  3.275000,  $
          3.350000,  3.450000,  3.550000,  3.650000,  3.775000,  $
          3.925000,  4.100000,  4.300000,  4.500000,  4.700000,  $
          4.900000]
   n_wls = size(wls, /dimensions)

   Dms = [3.000000,  7.500000,  15.00000,  25.00000,  35.00000,  $
          45.00000,  60.00000,  80.00000,  100.0000,  130.0000,  $
          175.0000,  225.0000,  275.0000,  350.0000,  450.0000,  $
          550.0000,  650.0000,  750.0000,  900.0000,  1150.000,  $
          1400.000,  1750.000,  2500.000,  3500.000]
   n_Dms = size(Dms, /dimensions)

   Des = [0.6207641, 1.551911,  3.103821,  5.173034,  7.242249,  $
          9.311462,  12.41528,  16.55371,  20.69214,  26.89978,  $
          36.21124,  46.55731,  56.90338,  72.42247,  93.11462,  $
          113.8068,  134.4989,  155.1910,  186.2292,  237.9595,  $
          289.6899,  362.1125,  517.3034,  724.2249]
   n_Des = size(Des, /dimensions)

   i_wl1 = value_locate(wls, wl)
   i_wl1 = (n_wls - 2) < (0 > i_wl1)
   i_wl2 = i_wl1 + 1

   a_wl = (wl - wls[i_wl1]) / (wls[i_wl2] - wls[i_wl1])

   i_De1 = value_locate(Des, De)
   i_De1 = (n_Des - 2) < (0 > i_De1)
   i_De2 = i_De1 + 1

   a_De = (De - Des[i_De1]) / (Des[i_De2] - Des[i_De1])

   P_De1 = fltarr(1801)
   P_De2 = fltarr(1801)
   P_wl1 = fltarr(1801)
   P_wl2 = fltarr(1801)

   baran_read_sw_index, lun, i_wl1, i_De1, n_Dms, $
                        Bext_De1, w_De1, g_De1, theta, P_De1
   baran_read_sw_index, lun, i_wl1, i_De2, n_Dms, $
                        Bext_De2, w_De2, g_De2, theta, P_De2
   baran_combine_props, a_De, Bext_De1, w_De1, g_De1, P_De1, $
                              Bext_De2, w_De2, g_De2, P_De2, $
                              Bext_wl1, w_wl1, g_wl1, P_wl1

   baran_read_sw_index, lun, i_wl2, i_De1, n_Dms, $
                        Bext_De1, w_De1, g_De1, theta, P_De1
   baran_read_sw_index, lun, i_wl2, i_De2, n_Dms, $
                        Bext_De2, w_De2, g_De2, theta, P_De2
   baran_combine_props, a_De, Bext_De1, w_De1, g_De1, P_De1, $
                              Bext_De2, w_De2, g_De2, P_De2, $
                              Bext_wl2, w_wl2, g_wl2, P_wl2

   baran_combine_props, a_wl, Bext_wl1, w_wl1, g_wl1, P_wl1, $
                              Bext_wl2, w_wl2, g_wl2, P_wl2, $
                              Bext, w, g, P
end


pro baran_read_sw_index, lun, i_wl, i_Dm, n_Dms, Bext, w, g, theta, P

   point_lun, lun, 0

   n_skip = 1ul + (i_wl * n_Dms + i_Dm) * 1808ul

   skip_lun, lun, n_skip, /lines

   readf, lun, wl
   readf, lun, Dm
   readf, lun, Bext
   readf, lun, Bscat
   readf, lun, w
   readf, lun, g

   skip_lun, lun, 1, /lines

   row = fltarr(2)
   for i=0,1801-1 do begin
      readf, lun, row
      theta[i] = row[0]
      P    [i] = row[1]
   endfor
end


pro baran_read_lw_interp, path, wl, Re, Bext, w, theta, P

   wls =  [4.00000, 4.25000, 4.50000, 4.75000, 5.00000, $
           5.25000, 5.50000, 5.75000, 6.00000, 6.25000, $
           6.50000, 6.75000, 7.00000, 7.25000, 7.50000, $
           7.75000, 8.00000, 8.25000, 8.50000, 8.75000, $
           9.00000, 9.25000, 9.50000, 9.75000, 10.0000, $
           10.2500, 10.5000, 10.7500, 11.0000, 11.2500, $
           11.5000, 11.7500, 12.0000, 12.2500, 12.5000, $
           12.7500, 13.0000, 13.2500, 13.5000, 13.7500, $
           14.0000, 14.2500, 14.5000, 14.7500]
   n_wls = size(wls, /dimensions)

   tags = ['lv06',  'lv07',  'lv08',  'lv09',  'll01',  $
           'll02',  'll03',  'll04',  'll05',  'll06',  $
           'll07',  'll08',  'll09',  'll10',  'll11',  $
           'll12',  'l001',  'l002',  'l003',  'l004',  $
           'l005',  'l006',  'l007',  'l008',  'l009',  $
           'l010',  'l011',  'l012',  'l013',  'l014',  $
           'l015',  'l016',  'l017',  'l018',  'l019',  $
           'l020',  'l021',  'l022',  'l023',  'l024',  $
           'l025',  'l026',  'l027',  'l028']

   Res  = [3.6, 5.5, 13.0, 22.1, 33.6, 35.8, 46.1, 77.4, 92.2]
   n_Res = size(Res, /dimensions)

   i_wl1 = value_locate(wls, wl)
   i_wl1 = (n_wls - 2) < (0 > i_wl1)
   i_wl2 = i_wl1 + 1

   a_wl = (wl - wls[i_wl1]) / (wls[i_wl2] - wls[i_wl1])

   i_Re1 = value_locate(Res, Re)
   i_Re1 = (n_Res - 2) < (0 > i_Re1)
   i_Re2 = i_Re1 + 1

   a_Re = (Re - Res[i_Re1]) / (Res[i_Re2] - Res[i_Re1])

   P_Re1 = fltarr(1801)
   P_Re2 = fltarr(1801)
   P_wl1 = fltarr(1801)
   P_wl2 = fltarr(1801)

   baran_read_lw_file, path, tags[i_wl1], i_Re1, $
                       Bext_Re1, w_Re1, theta, P_Re1
   baran_read_lw_file, path, tags[i_wl1], i_Re2, $
                       Bext_Re2, w_Re2, theta, P_Re2
   g_Re1 = !Values.F_NAN
   g_Re2 = !Values.F_NAN
   baran_combine_props, a_Re, Bext_Re1, w_Re1, g_Re1, P_Re1, $
                              Bext_Re2, w_Re2, g_Re2, P_Re2, $
                              Bext_wl1, w_wl1, g_wl1, P_wl1

   baran_read_lw_file, path, tags[i_wl2], i_Re1, $
                       Bext_Re1, w_Re1, theta, P_Re1
   baran_read_lw_file, path, tags[i_wl2], i_Re2, $
                       Bext_Re2, w_Re2, theta, P_Re2
   g_Re1 = !Values.F_NAN
   g_Re2 = !Values.F_NAN
   baran_combine_props, a_Re, Bext_Re1, w_Re1, g_Re1, P_Re1, $
                              Bext_Re2, w_Re2, g_Re2, P_Re2, $
                              Bext_wl2, w_wl2, g_wl2, P_wl2

   baran_combine_props, a_wl, Bext_wl1, w_wl1, g_wl1, P_wl1, $
                              Bext_wl2, w_wl2, g_wl2, P_wl2, $
                              Bext, w, g, P
end


pro baran_read_lw_file, path, tag, i_sd, Bext, w, theta, P

   filename = path + '/' + 'nir' + '/' + 'sbxt_' + tag + '_sd' + $
              string(i_sd+1, format='(i0)') + '.dat'
   openr, lun, filename, /get_lun
   readf, lun, wl, Bext, w, fdelta

   row = fltarr(2)
   for i=0,1801-1 do begin
      readf, lun, row
      theta[i] = row[0]
      P    [i] = row[1]
   endfor

   free_lun, lun
end
