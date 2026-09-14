; Record which source file IDL actually resolves each legacy cloud/aerosol routine
; from, before any calculation runs.
;
; This matters because the aerosol entry points exist only in the preserved
; working tree, while the repository holds the V2.1 versions of the routines
; they call.  Several of those routines differ between the two generations
; (notably load_srfstrarr, where the meaning of srf_quad 0 and 1 is swapped),
; so the provenance of every resolved routine is part of the reference
; capture, not an implementation detail.
;
; Each resolved routine is printed as
;     PROVENANCE <name> <file>
; and each unresolved one as
;     PROVENANCE <name> UNRESOLVED
; so the capture runner can record them verbatim.  Nothing is modified.

pro aerosol_legacy_provenance

  ; Whether each name is a procedure or a function is a property of the source
  ; generation in use, so both forms are tried rather than assumed.
  names = ['create_orac_aerosol_lut_wrapper', 'create_orac_aerosol_lut', $
           'create_orac_cloud_lut_wrapper', 'create_orac_cloud_lut', $
           'generate_scattering_properties', 'load_inststr', 'load_lutstr', $
           'load_mmdat', 'load_atmstr', 'load_gasstr', 'load_srfstrarr', $
           'read_srfstr', 'setup_disort', 'call_disort', 'write_v2_lut', $
           'create_bwgp', 'read_ri', 'segment', 'legpexp', 'lut_quadrature', $
           'quadrature', 'shift_quadrature', 'mie_size_dist_new', $
           'load_solar_spectrum', 'integrate_trapeziodal', 'bbconstants']

  for i = 0, n_elements(names) - 1 do begin
    name = names[i]
    path = 'UNRESOLVED'

    catch, error_status
    if (error_status eq 0) then begin
      resolve_routine, name, /no_recompile
      info = routine_info(name, /source)
      if (n_elements(info) gt 0) then path = info[0].path
    endif
    catch, /cancel

    if (path eq 'UNRESOLVED' or strlen(strtrim(path, 2)) eq 0) then begin
      catch, error_status
      if (error_status eq 0) then begin
        resolve_routine, name, /is_function, /no_recompile
        info = routine_info(name, /functions, /source)
        if (n_elements(info) gt 0) then path = info[0].path
      endif
      catch, /cancel
    endif

    if (strlen(strtrim(path, 2)) eq 0) then path = 'UNRESOLVED'
    print, 'PROVENANCE ' + name + ' ' + strtrim(path, 2)
  endfor

  help, /dlm, output = dlm_report
  for i = 0, n_elements(dlm_report) - 1 do begin
    line = dlm_report[i]
    if (strpos(strupcase(line), 'MIE') ge 0) or (strpos(strupcase(line), 'DISORT') ge 0) then $
      print, 'DLM ' + strtrim(line, 2)
  endfor

  print, 'PROVENANCE_COMPLETE'
end
