; Capture the compact intermediate arrays already saved by the current cloud
; generator.  This is a validation-only post-processing helper: it does not
; participate in production LUT generation and does not alter any calculations.

function capture_shape_text, value
  dimensions = size(value, /dimensions)
  if n_elements(dimensions) eq 0 then return, 'scalar'
  return, strjoin(strtrim(dimensions, 2), ',')
end

pro capture_array_text, output_dir, name, value, units
  filename = output_dir + '/' + name + '.txt'
  openw, lun, filename, /get_lun
  printf, lun, 'name=' + name
  printf, lun, 'shape=' + capture_shape_text(value)
  printf, lun, 'count=' + strtrim(n_elements(value), 2)
  printf, lun, 'units=' + units
  flattened = reform(double(value), n_elements(value))
  for i = 0, n_elements(flattened) - 1 do $
    printf, lun, flattened[i], format='(ES24.16)'
  free_lun, lun
end

pro capture_legacy_test_reference, output_dir
  save_file = output_dir + '/scatfile.sav'
  if ~file_test(save_file, /read) then $
    message, 'Expected current-generator scattering save file was not found: ' + save_file

  restore, save_file
  capture_array_text, output_dir, 'bext550', bext550, 'legacy-defined'
  capture_array_text, output_dir, 'w550', w550, 'dimensionless'
  capture_array_text, output_dir, 'g550', g550, 'dimensionless'
  capture_array_text, output_dir, 'phs550', phs550, 'legacy phase-function units'
  capture_array_text, output_dir, 'amom550', amom550, 'Legendre-moment convention from legacy legpexp'
  capture_array_text, output_dir, 'bextrat', bextrat, 'relative to 0.55 micron'
  capture_array_text, output_dir, 'bext', bext, 'legacy-defined'
  capture_array_text, output_dir, 'w', w, 'dimensionless'
  capture_array_text, output_dir, 'g', g, 'dimensionless'
  capture_array_text, output_dir, 'vavg', vavg, 'legacy-defined'
  capture_array_text, output_dir, 'phs', phs, 'legacy phase-function units'
  capture_array_text, output_dir, 'amom', amom, 'Legendre-moment convention from legacy legpexp'
  openw, lun, output_dir + '/intermediate_manifest.txt', /get_lun
  printf, lun, 'source=scatfile.sav written by current create_orac_cloud_lut'
  printf, lun, 'reference_wavelength_microns=0.55'
  printf, lun, 'phase_function_and_moments=preserved with legacy array ordering'
  printf, lun, 'disort_layer_inputs=not exposed by the unmodified production routine'
  printf, lun, 'disort_outputs=represented by the generated V2 NetCDF operators'
  free_lun, lun
end
