; This is an example call to DISORT from IDL.

; Radiances are returned to matrix UU where
; UU((0:NUMU-1),(0:NTAU-1),(0:NPHI-1))
;
; To run from this directory set the IDL_DLM_PATH environmental variable, for
; example for bash, with
;
; > set IDL_DLM_PATH=src:$IDL_DLM_PATH
;
; or within IDL with
;
; IDL> PREF_SET, 'IDL_DLM_PATH', 'src:<IDL_DEFAULT>', /COMMIT
;
; and run with
; IDL> .r ExampleDISORT.pro
;
; PREF_SET, 'IDL_DLM_PATH', 'src:<IDL_DEFAULT>', /COMMIT
;
; HISTORY:
; 05/02/13, S Marsh: Original implementation.
;

;********** DISORT VARIABLES *******************************************

NLYR = 1                ; NO OF LAYERS
DTAUC = FLTARR(NLYR)    ; NLYR LAYERS OPTICAL DEPTH
SSALB = FLTARR(NLYR)    ; NLYR LAYERS SSA
NSTR = 16               ; NUMBER OF STREAMS
NMOM = NSTR             ; NSTR PHASE MOMENTS
PMOM = FLTARR(NMOM+1,NLYR); PHASE MOMENTS
TEMPER = FLTARR(NLYR)   ; TEMPERATURE FIELD (LEFT ZERO)
WVNMLO = 0.0            ; ONLY USED IF PLANK TRUE
WVNMHI = 0.0            ; ONLY USED IF PLANK TRUE
USRTAU = 1              ; RADIANT QUANTITY OPTICAL DEPTHS SPECIFIED
NTAU = 2                ; NUMBER OPTICAL DEPTHS SPECIFIED (TOP AND BOTTOM)
UTAU = FLTARR(NTAU)     ; USER DEFINED OPTICAL DEPTHS
USRANG = 1              ; RETURN POLAR ANGLES TO BE SPECIFIED
NUMU = 2                ; NUMBER OF POLAR ANGLES
UMU = FLTARR(NUMU)      ; COMPUTATIONAL POLAR ANGLES
NPHI = 1                ; NUMBER OF RETURN AZIMUTH ANGLES
PHI = fltarr(NPHI)      ; AZIMUTHAL OUTPUT ANGLES
IBCND = 0               ; GENERAL CASE BOUNDARY CONDITIONS
FBEAM = 1.0             ; INTENSITY OF INCIDENT BEAM
UMU0 = 1.0              ; POLAR ANGLE INCIDENT BEAM
PHI0 = 0.0              ; AZIMUTH OF INCIDENT BEAM
FISOT = 0               ; TOP BOUNDARY INTENSITY
LAMBER = 1              ; LAMBERTIAN BOTTOM SURFACE
ALBEDO = 0.0            ; BOTTOM SURFACE ALBEDO
BTEMP = 0.0             ; ONLY USED IF PLANK TRUE (1)
TTEMP = 0.0             ; ONLY USED IF PLANK TRUE (1)
TEMIS = 0.0             ; ONLY USED IF PLANK TRUE (1)
PLANK = 0               ; IGNORE THERMAL SOURCES
ONLYFL = 0              ; RETURN FLUXES AND INTENSITIES
ACCUR = 0.000001        ; CONVERGENCE ACCURACY
PRNT = [0,0,0,0,0]      ; LOGICAL PRINT FLAGS
MAXCLY = NLYR           ; MAX DIM OF DTAUC SSALB TEMPER
MAXULV = NTAU           ; DIM OF UTAU RFLDIR RFLDN FLUP DFDT
MAXUMU = NUMU           ; DIM OF UMU
MAXPHI = NPHI           ; DIM OF PHI
MAXMOM = NMOM           ; FIRST DIM OF PMOM

;********** END OF DISORT VARIABLES ************************************

g = fltarr(NLYR)
dtauc(0) = 1
ssalb(0) = 0.9
g(0) = 0.61
UMU(0) = -1.0
UMU(1) =  1.0
UTAU(0) = 0.0           ; TOP OF ATMOSPHERE
UTAU(1) = TOTAL(DTAUC)  ; BOTTOM OF ATMOSPHERE

; Compute the expansion of the phase functions in terms of Legendre
; polynomials for each layer.  If G(I) equal to zero then a Rayleigh
; phase function is used other a Henyey-Greenstein phase function is
; used with an asymmetry factor equal to G(I).
FOR i=0,NLYR-1 DO BEGIN
    IF (G(I) EQ 0) THEN BEGIN
        GETMOM,2,G(I),NMOM,PM
    ENDIF
    IF (G(I) NE 0) THEN BEGIN
        GETMOM,3,G(I),NMOM,PM
    ENDIF
    PMOM(*,I) = PM
ENDFOR

; Call DISORT
DISORT, NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO, WVNMHI, $
        USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, UMU, NPHI, PHI, IBCND, $
        FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, $
        PLANK, ONLYFL, ACCUR, PRNT, MAXCLY, MAXULV, MAXUMU, MAXPHI, $
        MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED

print,uu

END
