;+
; pro write_orac_optprop_luts
; This procedure writes a NetCDF file containing LUTs for various
; microphysical and optical properties from an ORAC aerosol class, as a
; function of aerosol effective radius. These LUTs can then be used to 
; back-out the assumed aerosol properties corresponding to a given ORAC
; retrieval with knowledge of the aerosol type and retrieved effective
; radius alone.
; INPUT PARAMETERS
; outfile  The full file name of the NetCDF file to be generated
; Efr      An array specifying the effective radius LUT points
; Comp_name  The name of each component in the aerosol type
; Comp_fname  The filename for this component in the GADS dataset
; Comp_MR  The number mixing ratio of each component, as a function of
;          effective radius.
; Comp_RM  The mode radius of each component, as a function of effective
;          radius, in units of microns.
; Comp_S   The log-normal spread of each component, in units of microns.
;          This quantity is defined by (where N(r) is the number of
;          particles at radius r):
;              stdev(ln(N(r))) = ln(Comp_S)
; Wv       The wavelengths at which the refractive indexes have been
;          tabulated
; Comp_RI  The complex refractive indices of each component
; SSA      The single scatter albedo of the type as a whole at 0.55
;          microns, as a function of effective radius
; Bext     The extinction coefficent of the type as a whole at 0.55
;          microns, in units of km^-1, as a function of effective 
;          radius
; g        The asymmetry parameter of the type as a whole at 0.55
;          microns, as a function of effective radius
; SA       The scattering angles at which the phase function has been
;          tabulated, in units of radians.
; Phs      The phase function of the type as a whole at 0.55 microns,
;          as a function of effective radius
; Bext_c   The extinction coefficient for each component as a function
;          of effective radius and wavelength (optional)
; SSA_c    The single scatter albedo for each component as a function
;          of effective radius and wavelength (optional)
; HISTORY
; 04/06/2010 G Thomas, Original.
; 21/07/2014 G Thomas, Added Bext_c and SSA_c
; 20/03/2020 G Thomas, Added to ORAC-CC/create_orac_lut Git repository
;-
pro write_orac_optprop_luts, outfile, type_name, Efr, Comp_type,       $
                             Comp_name, Comp_name2, Comp_dist,         $
                             Comp_fname, Comp_MR, Comp_RM, Comp_S, Wv, $
                             Comp_RI, SSA, Bext, g, SA, Phs, Bext_c,   $
                             SSA_c

  nCmp = n_elements(Comp_name)
  nEfr = n_elements(Efr)
  nWv  = n_elements(Wv)
  nSA  = n_elements(SA)

  references = ['https://github.com/ORAC-CC/create_orac_lut']

  fid = ncdf_create(outfile,/clobber)
  ncdf_control,0,/verbose

; Define global attributes
  ncdf_attput, fid, /global, 'title', $
               'ORAC aerosol properties LUT'
  ncdf_attput, fid, /global, 'aerosol_type', $
               strtrim(type_name,2)
  ncdf_attput, fid, /global, 'history', $
               'LUT calculated on '+systime(/utc)
  ncdf_attput, fid, /global, 'references', $
               strjoin(references,', ')

; Define the dimensions to be used
  d1id = ncdf_dimdef(fid, 'component', nCmp)
  d2id = ncdf_dimdef(fid, 'effective_radius', nEfr)
  d3id = ncdf_dimdef(fid, 'wavelength', nWv)
  d4id = ncdf_dimdef(fid, 'scattering_angle', nSA)

; Define the "axis" variables corresponding to each of the above
; dimensions, and their attributes
  Cmpid = ncdf_vardef(fid, 'component', [d1id], /char)
     ncdf_attput, fid, Cmpid, 'long_name', /char, $
                  'Aerosol components which make up this type'
     ncdf_attput, fid, Cmpid, 'types', /char, $
                  strjoin(comp_type,', ')
     ncdf_attput, fid, Cmpid, 'names', /char, $
                  strjoin(Comp_name,', ')
     if max(comp_name2) ne '' then $
        ncdf_attput, fid, Cmpid, 'names', /char, $
                     strjoin(Comp_name2,', ')
     ncdf_attput, fid, Cmpid, 'size_distribution', /char, $
                  strjoin(Comp_dist,', ')
     ncdf_attput, fid, Cmpid, 'file_names', /char, $
                  strjoin(Comp_fname,', ')
  Efrid = ncdf_vardef(fid, 'effective_radius', [d2id], /float)
     ncdf_attput, fid, Efrid, 'long_name', /char, $
                  'aerosol effective radius LUT points'
     ncdf_attput, fid, Efrid, 'units', /char, 'microns'
  Wvid  = ncdf_vardef(fid, 'wavelength', [d3id], /float)
     ncdf_attput, fid, Wvid, 'long_name', /char, $
                  'Refractive index wavelength grid'
     ncdf_attput, fid, Wvid, 'units', /char, 'microns'
  SAid  = ncdf_vardef(fid, 'scattering_angle', [d4id], /float)
     ncdf_attput, fid, SAid, 'long_name', /char, $
                  'Phase function scattering angle grid'
     ncdf_attput, fid, SAid, 'units', /char, 'radians'

; Define the LUT variables
  MRid  = ncdf_vardef(fid, 'mixing_ratio', [d2id,d1id], /float)
     ncdf_attput, fid, MRid, 'long_name', /char, $
                  'Number mixing ratio of each component'
  RMid  = ncdf_vardef(fid, 'mode_radius', [d2id,d1id], /float)
     ncdf_attput, fid, RMid, 'long_name', /char, $
                  'Mode (median) radius of each component'
     ncdf_attput, fid, RMid, 'units', /char, 'microns'
  Sid   = ncdf_vardef(fid, 'sigma', [d1id], /float)
     ncdf_attput, fid, Sid, 'long_name', /char, $
                  'Log-normal spread of the size distribution of each component'
     ncdf_attput, fid, Sid, 'units', /char, 'microns'
  RInid = ncdf_vardef(fid, 'refractive_index', [d3id,d1id], /float)
     ncdf_attput, fid, RInid, 'long_name', /char, $
                  'Real part of the complex refractive index of each component'
  RIiid = ncdf_vardef(fid, 'absorption_coef', [d3id,d1id], /float)
     ncdf_attput, fid, RIiid, 'long_name', /char, $
                  'Imaginary part of the complex refractive index of each component'
  SSAid = ncdf_vardef(fid, 'single_scattering_albedo', [d3id,d2id], /float)
     ncdf_attput, fid, SSAid, 'long_name', /char, $
                  'Single scattering albedo of this type at 0.55 microns'
  Extid = ncdf_vardef(fid, 'extinction_coef', [d3id,d2id], /float)
     ncdf_attput, fid, Extid, 'long_name', /char, $
                  'Extinction coefficient of this type at 0.55 microns'
     ncdf_attput, fid, Extid, 'units', /char, 'km^-1'
  gid   = ncdf_vardef(fid, 'asymmetry_param', [d3id,d2id], /float)
     ncdf_attput, fid, gid, 'long_name', /char, $
                  'Asymmetry parameter of this type at 0.55 microns'
  Phsid = ncdf_vardef(fid, 'phase_func', [d4id,d3id,d2id], /float)
  ncdf_attput, fid, Phsid, 'long_name', /char, $
               'Scattering phase function of this type at 0.55 microns'
  if n_elements(Bext_c) gt 0 then begin
     Ext_cid = ncdf_vardef(fid, 'comp_extinction_coef', [d3id,d1id,d2id], /float)
     ncdf_attput, fid, Ext_cid, 'long_name', /char, $
                  'Extinction coefficient of each component at each wavelength'
     ncdf_attput, fid, Ext_cid, 'units', /char, 'km^-1'
  endif
  if n_elements(SSA_c) gt 0 then begin
     SSA_cid = ncdf_vardef(fid, 'comp_single_scattering_albedo', [d3id,d1id,d2id], /float)
     ncdf_attput, fid, SSA_cid, 'long_name', /char, $
                  'Single scatter albedo of each component at each wavelength'
  endif
  
; Switch the file to data-mode and write data into all variables
  ncdf_control, fid, /endef
  ncdf_varput, fid, Cmpid, strjoin(Comp_name,', ')
  ncdf_varput, fid, Efrid, Efr
  ncdf_varput, fid, Wvid, Wv
  ncdf_varput, fid, SAid, SA

  ncdf_varput, fid, MRid, transpose(Comp_MR)
  ncdf_varput, fid, RMid, transpose(Comp_RM)
  ncdf_varput, fid, Sid, Comp_S
  ncdf_varput, fid, RInid, float(Comp_RI)
  ncdf_varput, fid, RIiid, imaginary(Comp_RI)
  ncdf_varput, fid, SSAid, reform(SSA)
  ncdf_varput, fid, Extid, reform(Bext)
  ncdf_varput, fid, gid, reform(g)
  ncdf_varput, fid, Phsid, transpose(reform(Phs),[2,0,1])
  if n_elements(Bext_c) gt 0 then ncdf_varput, fid, Ext_cid, reform(Bext_c)
  if n_elements(SSA_c) gt 0 then ncdf_varput, fid, SSA_cid, reform(SSA_c)

; All done, close the file.
  ncdf_close, fid

end
