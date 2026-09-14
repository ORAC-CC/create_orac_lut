; Run one captured DISORT state through the current production IDL wrappers.
; This is a small diagnostic only; it does not generate a LUT.

pro read_float_file, filename, values
  openr, lun, filename, /get_lun
  readu, lun, values
  free_lun, lun
end

pro run_worst_r0v_idl
  root = '/home/g/grainger/project-oraclut/validation/diagnostics'
  nlyr = 45L
  nstr = 60L
  nmoments = 1000L
  nsatzen = 10L
  nrelazi = 11L

  dtauc = fltarr(nlyr)
  ssalb = fltarr(nlyr)
  pmom = fltarr(nmoments, nlyr)
  utau = fltarr(2)
  umu = fltarr(20)
  phi = fltarr(nrelazi)
  read_float_file, root + '/worst_r0v_disort_state_ifort_dtau.bin', dtauc
  read_float_file, root + '/worst_r0v_disort_state_ifort_ssalb.bin', ssalb
  read_float_file, root + '/worst_r0v_disort_state_ifort_pmom.bin', pmom
  read_float_file, root + '/worst_r0v_disort_state_ifort_utau.bin', utau
  read_float_file, root + '/worst_r0v_disort_state_ifort_umu.bin', umu
  read_float_file, root + '/worst_r0v_disort_state_ifort_phi.bin', phi

  setup_disort, nstr, nlyr, nsatzen, nrelazi, nmoments
  fbeam = 100.0
  umu0 = cos(50.0 * !dtor)
  fisot = 0.0
  call_disort, dtauc, ssalb, pmom, utau, umu, phi, fbeam, umu0, fisot, $
               rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed

  print, 'IDL targeted DISORT completed'
  print, 'IDL R_0v target = ', uu[11, 0, 10] * !pi * 0.01
  print, 'IDL RFldir = ', rfldir
  print, 'IDL RFldn = ', rfldn
  print, 'IDL Flup = ', flup
end

