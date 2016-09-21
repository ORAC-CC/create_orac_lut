;+
 pro dubovik_lognormal_multiple_eps, base_path, number_conc, rm, s, wn, ri, eps, neps, dqv=dqv, bext, bsca, w, g, ph, silent = silent, no_mie=no_mie, renorm_ph=renorm_ph
;
; As dubovik_lognormal, but eps and neps must be vectors with equal numbers
; of elements. If total(neps) != 1, then it will be renormalised before
; continuing.
;
; Andy Smith (smith@atm.ox.ac.uk)
;-

    n = n_elements(eps)
    nt= n_elements(dqv)

    ; Initial checks
    if n_elements(neps) ne n then message, $
      'EPS and NEPS must have the same number of elements'

    q = where(neps lt 0.0d0)
    if q[0] ne -1 then message,'NEPS values cannot be negative!'

    if total(neps) ne 1.0 and ~keyword_set(silent) then message,/continue,$
      'total(NEPS) ='+string(total(neps))+' has been renormalised'
    neps_norm = neps / total(neps)

    if n_elements(rm) ne 1 and n_elements(rm) ne n then message,$
      'Size of RM and EPS are not consistent'
    if n_elements(S)  ne 1 and n_elements(s)  ne n then message,$
      'Size of S and EPS are not consistent'
    if n_elements(RI) ne 1 and n_elements(RI) ne n then message,$
      'Size of RI and EPS are not consistent'

    if n_elements(rm) eq 1 then rm_ = replicate(rm,n) else rm_= rm
    if n_elements(S)  eq 1 then S_  = replicate(S, n) else S_ = S
    if n_elements(RI) eq 1 then RI_ = replicate(RI,n) else RI_= RI

    bext_arr = dblarr(n)
    bsca_arr = dblarr(n)
    w_arr    = dblarr(n)
    g_arr    = dblarr(n)
    if nt gt 0 then ph_arr   = dblarr(nt, n)

    for i = 0, n-1 do begin

        if neps_norm[i] gt 0.0d0 then begin

            dubovik_lognormal, base_path, number_conc, rm_[i], s_[i], wn, ri_[i], $
                               eps[i], dqv=dqv, silent=silent, no_mie=no_mie,$
                               bext_, bsca_, w_, g_, ph_, renorm_ph=renorm_ph

            bext_arr[i] = bext_
            bsca_arr[i] = bsca_
            w_arr[   i] = w_
            g_arr[   i] = g_
            if nt gt 0 then ph_arr[*,i] = ph_
        endif ; else begin ; Doesn't matter what this is since neps[i]=0

    endfor
;     fff='(a,1x,f7.3,23f6.3)'
;     print, 'n   ',neps_norm,format=fff
;     print, 'bsc ',bsca_arr,format=fff
;     print, 'g   ',g_arr,format=fff
;     print, 'w   ',w_arr,format=fff


    bext = total(bext_arr * neps_norm)
    bsca = total(bsca_arr * neps_norm)
;    w    = total(bsca_arr * neps_norm) / total(bext_arr * neps_norm)
    w    = bsca / bext
    g    = total(bsca_arr * neps_norm * g_arr) / bsca

    if nt gt 0 then begin
        ph   = dblarr(nt)
        for j = 0, n_elements(dqv) - 1 do $
          ph[j] = total(bsca_arr * neps_norm * ph_arr[j,*]) / bsca
    endif
    return
end
