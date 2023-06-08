; This library contains two routines which set-up and then call the DISORT DLM
; to do the actual radiative transfer calculations.
;
; HISTORY:
; 21/06/13, G Thomas: Original version.


; procedure create_orac_lut
;
; Define a common block containing the DISORT arguments which are common to all
; calls of DISORT.
;
; HISTORY:
; 21/06/13, G Thomas: Original version.
; 16/09/13, G Thomas: Added a keyword for setting the surface albedo to
;                     setup_disort.

pro setup_disort, NStreams, NLayers, NSatzen, NRelAzi, NMoments, alb=alb

;  Define a common block containing the fixed variables needed by DISORT
   common DISORT_VARS, NLYR, NMOM, NTAU, NUMU, NPHI, NSTR, MAXCLY, MAXULV, $
                       MAXUMU, MAXPHI, MAXMOM, USRTAU, USRANG, PHI0, IBCND, $
                       LAMBER, ALBEDO, ONLYFL, ACCUR, PRNT, HEADER

;  DISORT requires most array sizes to be passed to it twice. The "MAX***"
;  values define the size of arrays allocated by DISORT.
;  The "N***" define the number of values being passed to DISORT or requested
;  from it.
   NLYR   = NLayers            ; Actual number of layers
   NMOM   = NMoments-1         ; Actual number of phase moments
   NTAU   = 2                  ; Actual number of user output layers
   NUMU   = NSatZen*2          ; Actual number of output view zeniths
   NPHI   = NRelAzi            ; Actual number of output azimuths
   NSTR   = NStreams           ; Number of streams to calculate
   MAXCLY = NLayers            ; Maximum number of layers
   MAXULV = 2                  ; Maximum number of output layers
   MAXUMU = NSatZen*2          ; Maximum number of output view zeniths
   MAXPHI = NRelAzi            ; Maximum number of output azimuths
   MAXMOM = NMoments-1         ; Maximum number of phase moments

;  WVNMLO = 0.0                ; Dummy wavelength interval
;  WVNMHI = 0.0                ;   "        "        "
   USRTAU = 1                  ; Flag for user defined output levels
   USRANG = 1                  ; Flag for user defined output angles
   PHI0   = 0.0                ; Azimuth angle of incident beam
   IBCND  = 0                  ; Boundary condition flag
   LAMBER = 1                  ; Flag for Lambertian bottom boundary
   ALBEDO = 0.0                ; Albedo of bottom boundary
   if keyword_set(alb) then ALBEDO = alb
;  BTEMP  = 0.0                ; Bottom boundary temperature
;  TTEMP  = 0.0                ; Top boundary temperature
;  TEMIS  = 0.0                ; Emissivity of top boundary
;  TEMPER = fltarr(NLYR)       ; Temperature of each layer
;  PLANK  = 0                  ; Flag for enabling emission calculations
   ONLYFL = 0                  ; Flag for disabling intensity output
   ACCUR  = 1e-8               ; Convergence criteria
   PRNT   = [0,0,0,0,0]        ; Printing flags
   HEADER = string(replicate(' ',127),format='(127a)')
                               ; Header string placeholder
end