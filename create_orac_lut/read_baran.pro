; This library contains procedures for reading Anthony Baran's ice cloud optical
; properties. These properties are described by:
;
; A. J. Baran and S. Havemann. The dependence of retrieved cirrus ice-crystal
; effective dimension on assumed ice-crystal geometry and size-distribution
; function at solar wavelengths. Quarterly Journal of the Royal Meteorological
; Society, 130(601):21532167, 2004. doi:10.1256/qj.03.154.
;
; History:
; ??/??/15, G. McGarragh: Initial implementation.
; 12/10/16, G. McGarragh: Properly integrate the SW scattering properties over
;    size distribution and significant refactoring.
; 13/10/16, G. McGarragh: Add init_baran() so that the large solar files need
;    only be read once per session.


function baran_n_Dm
   return, 24
end


function baran_n_Re
   return, 9
end


function baran_n_theta
   return, 1801
end


function baran_Re
   return, [3.6, 5.5, 13.0, 22.1, 33.6, 35.8, 46.1, 77.4, 92.2]
end


function baran_map_i_9_to_i_30_subset, i

   index = [17, 7, 6, 8, 4, 9, 18, 12, 15]

   return, index[i]
end


; procedure init_baran
;
; Read and interpolate Anthony Baran's ice cloud optical properties for a set of
; wavelengths and effective radii.
;
; INPUT ARGUMENTS:
; path (string) Path to the base directory containing the optical properties.
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; baran (structure) Structure to be passed to each call of read_baran().

pro init_baran, path, baran

   n_wl_baran = 54
   n_Dm_baran = baran_n_Dm()
   n_Re_baran = baran_n_Re()

   n_theta_baran = baran_n_theta()

   part_V = [1.695980, 26.49968, 211.9975, 981.4697, 2693.153, 5723.932, $
             13567.84, 32160.80, 62814.06, 138002.5, 336644.1, 715491.4, $
             1306336., 2693153., 5723931., 1.0450690E+07, 1.7250314E+07, $
             2.6499678E+07, 4.5791448E+07, 9.5532328E+07, 1.7236178E+08, $
             3.3664413E+08, 9.8146970E+08, 2.6931530E+09]

   part_A = [2.732084, 17.07552, 68.30209, 189.7280, 371.8670, 614.7189, $
             1092.833, 1942.815, 3035.649, 5130.246, 9296.675, 15367.97, $
             22957.09, 37186.70, 61471.88, 91828.38, 128256.1, 170755.3, $
             245887.5, 401464.6, 594987.2, 929667.4, 1897281., 3718670.]

   theta1  = fltarr(n_theta_baran)
   theta2  = fltarr(n_theta_baran)
   P1      = fltarr(n_theta_baran)

   part_Cext = fltarr(n_Dm_baran, n_wl_baran)
   part_Csca = fltarr(n_Dm_baran, n_wl_baran)
   part_w    = fltarr(n_Dm_baran, n_wl_baran)
   part_g    = fltarr(n_Dm_baran, n_wl_baran)
   part_P    = fltarr(n_theta_baran, n_Dm_baran, n_wl_baran)

   filename = path + '/solar/solar_optical_phase_aggregate.dat'

   openr, lun, filename, /get_lun
   skip_lun, lun, 1ul, /lines

   flag = 1
   for i = 0, n_wl_baran - 1 do begin
      for j = 0, n_Dm_baran - 1 do begin
         baran_read_sw_opt_phase_agg, lun, i, j, Cext1, Csca1, w1, g1, theta2, P1

         part_Cext[j,i] = Cext1
         part_Csca[j,i] = Csca1
         part_w   [j,i] = w1
         part_g   [j,i] = g1
         part_P [*,j,i] = P1

         if flag then begin
            flag = 0
            theta1 = theta2
         endif
      endfor
   endfor

   free_lun, lun

   dist_Re   = fltarr(n_Re_baran)
   dist_Bext = fltarr(n_Re_baran, n_wl_baran)
   dist_Bsca = fltarr(n_Re_baran, n_wl_baran)
   dist_w    = fltarr(n_Re_baran, n_wl_baran)
   dist_g    = fltarr(n_Re_baran, n_wl_baran)
   dist_P    = fltarr(n_theta_baran, n_Re_baran, n_wl_baran)

   for i = 0, n_Re_baran - 1 do begin
      baran_read_sw_nsized, path, i, Dm, nx, Dl

      n_x_Dl = nx * Dl

      dist_Re[i] = total(part_V * n_x_Dl) / total(part_A * n_x_Dl) * 3. / 2. / 2.

      for j = 0, n_wl_baran - 1 do begin
         for k = 0, n_Dm_baran - 1 do begin
            dist_Bext[  i,j] += n_x_Dl[k] * part_Cext[  k,j]
            dist_Bsca[  i,j] += n_x_Dl[k] * part_Csca[  k,j]
            dist_g   [  i,j] += n_x_Dl[k] * part_g   [  k,j] * part_Csca[k,j]
            dist_P   [*,i,j] += n_x_Dl[k] * part_P   [*,k,j] * part_Csca[k,j]
         endfor
         dist_g [  i,j] /= dist_Bsca[i,j]
         dist_P [*,i,j] /= dist_Bsca[i,j]
      endfor
   endfor

   dist_w = dist_Bsca / dist_Bext

   baran = {baran, path  : path, $
                   Re    : dist_Re, $
                   Bext  : dist_Bext, $
                   Bsca  : dist_Bsca, $
                   w     : dist_w, $
                   g     : dist_g, $
                   theta : theta1, $
                   P     : dist_P}
end


; procedure read_baran
;
; Read and interpolate Anthony Baran's ice cloud optical properties for a set of
; wavelengths and effective radii.
;
; INPUT ARGUMENTS:
; baran (structure) Output from init_baran().
; wl    (array)     Wavelengths for which to obtain optical properties.
; Re    (array)     Effective radii for which to obtain optical properties.
; theta (array)     Angles to which to interpolate the phase function output.
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

pro read_baran, baran, wl, Re, Bext, w, g, theta, P

   n_theta_baran = baran_n_theta()

   n_wl    = n_elements(wl)
   n_Re    = n_elements(Re)
   n_theta = n_elements(theta)

   Bext    = fltarr(n_wl, n_Re)
   w       = fltarr(n_wl, n_Re)
   g       = fltarr(n_wl, n_Re)
   P       = fltarr(n_theta, n_wl, n_Re)

   theta1  = fltarr(n_theta_baran)
   P1      = fltarr(n_theta_baran)

   for i = 0, n_Re - 1 do begin
      for j = 0, n_wl - 1 do begin
         if wl[j] lt 4. then begin
            baran_interp_sw, wl[j], Re[i], baran.Bext, baran.w, baran.g, $
                             baran.P, Bext1, w1, g1, P1
         endif else begin
            baran_read_and_interp_lw, baran.path, wl[j], Re[i], Bext1, w1, $
                                      theta1, P1
            g1 = 0.
         endelse

         Bext[j,i] = Bext1
         w   [j,i] = w1
         g   [j,i] = g1
         P [*,j,i] = interpol(P1, baran.theta, theta)
      endfor
   endfor
end


pro baran_combine_props, a, Bext1, w1, g1, P1, Bext2, w2, g2, P2, $
                         Bext, w, g, P

   Bext = (1. - a) * Bext1 + a * Bext2
   w    = (1. - a) * w1    + a * w2
   g    = (1. - a) * g1    + a * g2

   n = n_elements(P)

   for i = 0, n[0] - 1 do begin
      P[i] = (1. - a) * P1[i] + a * P2[i]
   end
end


pro baran_interp_sw, wl, Re, Bext0, w0, g0, P0, Bext, w, g, P


   wl_baran = [0.2250000, 0.2750000, 0.3500000, 0.4250000, 0.4750000, 0.5500000, $
               0.6500000, 0.7500000, 0.8500000, 0.9500000, 1.050000,  1.150000,  $
               1.250000,  1.350000,  1.450000,  1.575000,  1.725000,  1.850000,  $
               1.925000,  1.975000,  2.050000,  2.150000,  2.250000,  2.350000,  $
               2.500000,  2.650000,  2.750000,  2.812000,  2.838000,  2.862000,  $
               2.888000,  2.912000,  2.938000,  2.963000,  2.987000,  3.013000,  $
               3.037000,  3.062000,  3.088000,  3.112000,  3.138000,  3.162000,  $
               3.188000,  3.225000,  3.275000,  3.350000,  3.450000,  3.550000,  $
               3.650000,  3.775000,  3.925000,  4.100000,  4.300000,  4.500000,  $
               4.700000,  4.900000]
   n_wl_baran = n_elements(wl_baran)

   Re_baran = baran_Re()
   n_Re_baran = n_elements(Re_baran)

   n_theta_baran = baran_n_theta()

   i_wl1 = value_locate(wl_baran, wl)
   i_wl1 = (n_wl_baran - 2) < (0 > i_wl1)
   i_wl2 = i_wl1 + 1

   a_wl = (wl - wl_baran[i_wl1]) / (wl_baran[i_wl2] - wl_baran[i_wl1])

   i_Re1 = value_locate(Re_baran, Re)
   i_Re1 = (n_Re_baran - 2) < (0 > i_Re1)
   i_Re2 = i_Re1 + 1

   a_Re = (Re - Re_baran[i_Re1]) / (Re_baran[i_Re2] - Re_baran[i_Re1])

   P_wl1 = fltarr(n_theta_baran)
   P_wl2 = fltarr(n_theta_baran)

   baran_combine_props, a_Re, Bext0[i_Re1, i_wl1], w0[i_Re1, i_wl1], $
                              g0[i_Re1, i_wl1], P0[*,i_Re1, i_wl1], $
                              Bext0[i_Re2, i_wl1], w0[i_Re2, i_wl1], $
                              g0[i_Re2, i_wl1], P0[*,i_Re2, i_wl1], $
                              Bext_wl1, w_wl1, g_wl1, P_wl1

   baran_combine_props, a_Re, Bext0[i_Re1, i_wl2], w0[i_Re1, i_wl2], $
                              g0[i_Re1, i_wl2], P0[*,i_Re1, i_wl2], $
                              Bext0[i_Re2, i_wl2], w0[i_Re2, i_wl2], $
                              g0[i_Re2, i_wl2], P0[*,i_Re2, i_wl2], $
                              Bext_wl2, w_wl2, g_wl2, P_wl2

   baran_combine_props, a_wl, Bext_wl1, w_wl1, g_wl1, P_wl1, $
                              Bext_wl2, w_wl2, g_wl2, P_wl2, $
                              Bext, w, g, P
end


pro baran_read_sw_nsized, path, i_dist, Dm, nx, Dl

   filename = path + '/solar/nsized' + $
              string(baran_map_i_9_to_i_30_subset(i_dist), format='(i0)') + '.dat'

   openr, lun, filename, /get_lun

   n = baran_n_Dm()

   Dm = fltarr(n)
   nx = fltarr(n)
   Dl = fltarr(n)

   row = fltarr(3)
   for i = 0, n - 1 do begin
      readf, lun, row
      Dm[i] = row[0]
      nx[i] = row[1]
      Dl[i] = row[2]
   endfor

   free_lun, lun
end


pro baran_skip_sw_opt_phase_agg, lun, i_wl, i_Dm

   point_lun, lun, 0

   skip_lun, lun, 1ul + (i_wl * baran_n_Dm() + i_Dm) * 1808ul, /lines
end


pro baran_read_sw_opt_phase_agg, lun, i_wl, i_Dm, Cext, Csca, w, g, theta, P

   readf, lun, wl
   readf, lun, Dm
   readf, lun, Cext
   readf, lun, Csca
   readf, lun, w
   readf, lun, g

   skip_lun, lun, 1, /lines

   n_theta = baran_n_theta()

   row = fltarr(2)
   for i = 0, n_theta - 1 do begin
      readf, lun, row
      theta[i] = row[0]
      P    [i] = row[1]
   endfor
end


pro baran_read_and_interp_lw, path, wl, Re, Bext, w, theta, P

   wl_baran =  [4.00000, 4.25000, 4.50000, 4.75000, 5.00000, 5.25000, $
                5.50000, 5.75000, 6.00000, 6.25000, 6.50000, 6.75000, $
                7.00000, 7.25000, 7.50000, 7.75000, 8.00000, 8.25000, $
                8.50000, 8.75000, 9.00000, 9.25000, 9.50000, 9.75000, $
                10.0000, 10.2500, 10.5000, 10.7500, 11.0000, 11.2500, $
                11.5000, 11.7500, 12.0000, 12.2500, 12.5000, 12.7500, $
                13.0000, 13.2500, 13.5000, 13.7500, 14.0000, 14.2500, $
                14.5000, 14.7500]
   n_wl_baran = n_elements(wl_baran)

   file_tags = ['lv06',  'lv07',  'lv08',  'lv09',  'll01',  'll02',  $
                'll03',  'll04',  'll05',  'll06',  'll07',  'll08',  $
                'll09',  'll10',  'll11',  'll12',  'l001',  'l002',  $
                'l003',  'l004',  'l005',  'l006',  'l007',  'l008',  $
                'l009',  'l010',  'l011',  'l012',  'l013',  'l014',  $
                'l015',  'l016',  'l017',  'l018',  'l019',  'l020',  $
                'l021',  'l022',  'l023',  'l024',  'l025',  'l026',  $
                'l027',  'l028']

   Re_baran = baran_Re()
   n_Re_baran = n_elements(Re_baran)

   n_theta_baran = baran_n_theta()

   i_wl1 = value_locate(wl_baran, wl)
   i_wl1 = (n_wl_baran - 2) < (0 > i_wl1)
   i_wl2 = i_wl1 + 1

   a_wl = (wl - wl_baran[i_wl1]) / (wl_baran[i_wl2] - wl_baran[i_wl1])

   i_Re1 = value_locate(Re_baran, Re)
   i_Re1 = (n_Re_baran - 2) < (0 > i_Re1)
   i_Re2 = i_Re1 + 1

   a_Re = (Re - Re_baran[i_Re1]) / (Re_baran[i_Re2] - Re_baran[i_Re1])

   P_Re1 = fltarr(n_theta_baran)
   P_Re2 = fltarr(n_theta_baran)
   P_wl1 = fltarr(n_theta_baran)
   P_wl2 = fltarr(n_theta_baran)

   baran_read_lw_file, path, file_tags[i_wl1], i_Re1, $
                       Bext_Re1, w_Re1, theta, P_Re1
   baran_read_lw_file, path, file_tags[i_wl1], i_Re2, $
                       Bext_Re2, w_Re2, theta, P_Re2
   g_Re1 = 0.
   g_Re2 = 0.
   baran_combine_props, a_Re, Bext_Re1, w_Re1, g_Re1, P_Re1, $
                              Bext_Re2, w_Re2, g_Re2, P_Re2, $
                              Bext_wl1, w_wl1, g_wl1, P_wl1

   baran_read_lw_file, path, file_tags[i_wl2], i_Re1, $
                       Bext_Re1, w_Re1, theta, P_Re1
   baran_read_lw_file, path, file_tags[i_wl2], i_Re2, $
                       Bext_Re2, w_Re2, theta, P_Re2
   g_Re1 = 0.
   g_Re2 = 0.
   baran_combine_props, a_Re, Bext_Re1, w_Re1, g_Re1, P_Re1, $
                              Bext_Re2, w_Re2, g_Re2, P_Re2, $
                              Bext_wl2, w_wl2, g_wl2, P_wl2

   baran_combine_props, a_wl, Bext_wl1, w_wl1, g_wl1, P_wl1, $
                              Bext_wl2, w_wl2, g_wl2, P_wl2, $
                              Bext, w, g, P

   Bext /= 1000.
end


pro baran_read_lw_file, path, tag, i_sd, Bext, w, theta, P

   filename = path + '/nir/sbxt_' + tag + '_sd' + $
              string(i_sd+1, format='(i0)') + '.dat'
   openr, lun, filename, /get_lun
   readf, lun, wl, Bext, w, fdelta

   n_theta = baran_n_theta()

   row = fltarr(2)
   for i = 0, n_theta - 1 do begin
      readf, lun, row
      theta[i] = row[0]
      P    [i] = row[1]
   endfor

   free_lun, lun
end
