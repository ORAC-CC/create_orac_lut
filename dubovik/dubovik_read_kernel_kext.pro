pro dubovik_read_kernel_kext, file, dd, grid1=grid1

   if strmid(file,3,4,/reverse) eq '.fix' then fix = 1 else fix = 0


   if ~ file_test(file) then begin
       message,'File '+file+'not found',/continue
       dd = !values.f_nan
       return
   endif

   ;----------------------------------------------------
   ; Get values of radius
   if ~ keyword_set(grid1) then begin
;       gridpath = strsplit(file,'/',/extract) + '/'
;       gridpath = '/'+strjoin(gridpath[0:n_elements(gridpath)-2])
       gridpath = file_dirname(file,/mark_dir)
       grid1 = gridpath+'grid1.dat'
       if fix eq 1 then grid1 = grid1+'.fix'
   endif

   if ~ file_test(grid1) then begin
       message, 'File '+grid1+' not found',/continue
       dd = !values.f_nan
       return
   endif


   openr, lun, grid1, /get_lun
   readf, lun, NR_, WVL_
   NR_ = long(NR_)
   radius = dblarr(NR_)
   readf, lun, radius
   close, lun
   free_lun, lun


   ;----------------------------------------------------
   ; Read main file

   openr, lun, file, /get_lun

   if fix then begin
       ; key_RD value
       readf, lun, key_RD
       key_RD = long(key_RD)

       ; Some values I'm not sure about
       readf, lun, what_1, what_2

       ; Epsilon number
       readf, lun, neps
       neps = long(neps)

       ; Description line
       tmp = ''
       readf, lun, tmp
       eps = fltarr(neps)
       eps_distrb = fltarr(neps)
       ; EPS values
       for i_eps = 0, neps-1 do begin
           readf, lun,tmpeps, tmpepsd
           eps[i_eps] = tmpeps
           eps_distrb[i_eps] = tmpepsd
       endfor


   endif else begin
       ; Rmin, rmax, epsilon
       readf, lun, rmin, rmax, eps
   endelse


   ; Number of elements in parameters
   readf, lun, NR
   NR = -long(NR)

   ; Range of refractive indices
   tmp = ''
   readf, lun, tmp
   readf, lun, tmp
   ; Number of refractive index values
   readf, lun,  Nn, Nk
   Nn =  long(Nn)
   Nk = -long(Nk)

   ; Now start on the values of Kext and Kabs
   Kext = dblarr(NN, NK, NR)
   Kabs = dblarr(NN, NK, NR)
   RI_n = fltarr(NN)
   RI_k = fltarr(NK)

   tmp_NR  = dblarr(NR)
   for in = 0, NN-1 do begin
       for ik = 0, NK-1 do begin

           if fix then readf,lun,element else readf, lun, element, epstmp
           readf, lun, wvl, tmp1, tmp2

           RI_n[in] = tmp1
           RI_k[ik] = tmp2

           readf, lun, tmp
           readf, lun, tmp_NR
           Kext[in, ik, *] = tmp_NR

           readf, lun, tmp
           readf, lun, tmp_NR
           Kabs[in, ik, *] = tmp_NR

       endfor
   endfor

   close, lun
   free_lun, lun

   if fix then begin

       dd = create_struct('RI_real',   RI_n,      $
                          'RI_imag',   RI_k,      $
                          'radius',    radius,    $
                          'wvl',       wvl,       $
                          'epsilon',   eps,       $
                          'eps_distrb',eps_distrb,$
                          'Kext',      Kext,      $
                          'Kabs',      Kabs,      $
                          'key_RD',    key_RD      )

   endif else begin

       dd = create_struct('RI_real', RI_n,  $
                          'RI_imag', RI_k,  $
                          'radius',  radius,$
                          'wvl',     wvl,   $
                          'epsilon', eps,   $
                          'Kext',    Kext,  $
                          'Kabs',    Kabs    )

   endelse



   return


end
