; Capture the per-layer DISORT inputs the legacy aerosol generator assembles
; for one state, without running DISORT and without modifying legacy source.
;
; This repeats, line for line, the layer arithmetic of create_orac_aerosol_lut
; (profile interpolation at line 340-343, Rayleigh/gas/aerosol combination at
; lines 440-500 of the preserved working tree) using the legacy loaders and the
; legacy optical properties saved by the real run in scatfile.sav, so the
; captured arrays are exactly those the legacy code would pass to call_disort.
;
; Run from a clean directory that contains an input_files symlink and no .pro
; files, with the preserved tree on IDL_PATH:
;   idl -e "aerosol_layer_capture, 'luts/.../scatfile.sav', 'legacy_state.sav', channel=1, efr=0.01, opd=1.0, prs=950.0"

pro aerosol_layer_capture, scatfile, outfile, channel=channel, efr=efr, opd=opd, prs=prs

  in_path   = 'input_files'
  instfile  = 'meteosat-10_seviri_v1.inst'
  mmfile    = 'aerosol_a79.mm'
  lutfile   = 'aerosol_test.lut'
  atmospheres = '2'
  atmfile   = 'mls.atm'

  ; --- identical loader sequence to create_orac_aerosol_lut -----------------
  load_inststr, in_path + '/inst/' + instfile, inststr, RequestedChannelID = [channel]
  load_lutstr, in_path + '/lut/' + lutfile, inststr.max_sat_zenith, lutstr, /include_pressure
  load_srfstrarr, inststr, in_path + '/sun/Gueymard2018.sssi', srfstrarr, 1, nwvl_max
  load_mmdat, in_path + '/microphysics/' + mmfile, mmstr
  load_atmstr, in_path + '/atm/' + atmfile, atmospheres, atmstr
  load_gasstr, atmospheres, inststr.platform, inststr.instrument, inststr.channelid, in_path + '/gas/', gasstr

  ; --- legacy optics for this run (bextrat[m,l,r], w[m,l,r], amom[*,m,l,r]) --
  restore, scatfile

  ; --- create_orac_aerosol_lut lines 340-343 ---------------------------------
  nlayers = atmstr.nlevels - 1
  hlayers = (atmstr.height[0:nlayers-1] + atmstr.height[1:nlayers]) / 2
  scatreltau_raw = INTERPOL(mmstr.rext, mmstr.height, hlayers)
  scatreltau = scatreltau_raw / total(scatreltau_raw)

  ; --- indices of the requested state ----------------------------------------
  l = 0                                       ; single requested channel
  m = 0                                       ; srf_quad = 1: one spectral point
  r = (where(abs(lutstr.efr - efr) eq min(abs(lutstr.efr - efr))))[0]
  a = (where(abs(lutstr.opd - opd) eq min(abs(lutstr.opd - opd))))[0]
  k = (where(abs(lutstr.prs - prs) eq min(abs(lutstr.prs - prs))))[0]

  ; --- create_orac_aerosol_lut lines 440-500 ---------------------------------
  columntauray = (atmstr.pressure[atmstr.nlevels-1] / lutstr.prs(k)) / $
                 (117.03*srfstrarr[*].wvl_centre ^4 - 1.316*srfstrarr[*].wvl_centre ^2)
  GasIndx = (WHERE(Gasstr[*].ChannelID eq inststr.ChannelID[l])) [0]
  GasLvl = gasstr[GasIndx].Tau_Gas
  RayLvl = ColumnTauRay[l]*exp(-0.1188*atmstr.Height - 0.00116*atmstr.Height^2)
  TauGas = GasLvl[lindgen(NLayers)+1] - GasLvl[lindgen(NLayers)]
  TauRay = RayLvl[lindgen(NLayers)+1] - RayLvl[lindgen(NLayers)]

  tauscat = lutstr.opd[a] * scatreltau * bextrat[m,l,r]
  dtau = taugas + tauray + tauscat
  totaltau = total(dtau)
  SSALB = (TauRay + w[m,l,r]*tauscat) / DTau
  bd = WHERE(SSALB gt 1.0)
  IF bd[0] ge 0 then SSALB[bd] = 1.0
  ASYM = FLTARR(NLayers)
  nonzero = WHERE(tauscat gt 0.0)
  IF nonzero[0] ge 0 then ASYM[nonzero] = g[m,l,r]
  PMom = FLTARR(NMom, NLayers)
  FOR h=0,NLayers-1 do begin
     IF ASYM[h] eq 0.0 then GETMOM, 2, 0.0, NMom-1, PM $
     ELSE begin
        GETMOM, 2, 0.0, NMom-1, mPM
        PM = (mPM*TauRay[h] + AMom[*,m,l,r]*w[m,l,r]*tauscat[h]) / $
             (TauRay[h] + w[m,l,r]*tauscat[h])
     endelse
     bd = WHERE(PM gt 1.0)
     IF bd[0] ge 0 then PM[bd] = 1.0
     PMom[*,h] = PM
  ENDFOR
  GETMOM, 2, 0.0, NMom-1, molecular_moments
  particle_moments = reform(AMom[*,m,l,r])
  ; ----------------------------------------------------------------------------

  state_channel = inststr.channelid[l]
  state_efr = lutstr.efr[r]
  state_opd = lutstr.opd[a]
  state_prs = lutstr.prs[k]
  state_wavelength = srfstrarr[l].wvl_centre
  particle_ssa = w[m,l,r]
  particle_bextrat = bextrat[m,l,r]
  particle_g = g[m,l,r]
  level_height = atmstr.height
  level_pressure = atmstr.pressure
  layer_height = hlayers
  profile_height = mmstr.height
  profile_rext = mmstr.rext
  gas_level = GasLvl
  tau_gas = TauGas
  column_tau_rayleigh = columntauray[l]
  rayleigh_level = RayLvl
  tau_rayleigh = TauRay
  tau_aerosol_scattering = w[m,l,r]*tauscat
  tau_aerosol_absorption = (1.0 - w[m,l,r])*tauscat
  ssa = SSALB
  pmom = PMom

  save, filename = outfile, state_channel, state_efr, state_opd, state_prs, state_wavelength, $
        particle_ssa, particle_bextrat, particle_g, level_height, level_pressure, layer_height, $
        profile_height, profile_rext, scatreltau_raw, scatreltau, gas_level, tau_gas, $
        column_tau_rayleigh, rayleigh_level, tau_rayleigh, tauscat, tau_aerosol_scattering, $
        tau_aerosol_absorption, dtau, totaltau, ssa, pmom, molecular_moments, particle_moments, nlayers, nmom
  print, 'CAPTURED ', outfile, '  layers=', strtrim(nlayers,2), ' moments=', strtrim(nmom,2), $
         '  state: ch', strtrim(state_channel,2), ' efr=', strtrim(state_efr,2), ' opd=', strtrim(state_opd,2), ' prs=', strtrim(state_prs,2)
end
