; Define a common block containing the DISORT arguments which are common to all
; calls of DISORT.


pro setup_disort4, NStreams, NLayers, NSatzen, NRelAzi, NMoments
;  Define a common block containing the fixed variables needed by DISORT
   common DISORT_VARS, NLYR, NMOM, NTAU, NUMU, NPHI, NSTR, $
     USRANG, USRTAU, IBCND, ONLYFL, PRNT, PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE, WVNMLO, WVNMHI, PHI0,  ALBEDO, BTEMP, TTEMP, TEMIS, EARTH_RADIUS, ACCUR, HEADER 

   print,'Setting-up DISORT 4'
   NLYR   = NLayers            ; Actual number of layers
   NMOM   = NMoments-1         ; Actual number of phase moments
   NTAU   = 2                  ; Actual number of user output layers
   NUMU   = NSatZen*2          ; Actual number of output view zeniths
   NPHI   = NRelAzi            ; Actual number of output azimuths
   NSTR   = NStreams           ; Number of streams to calculate
 
   USRTAU    = 1L 
   USRANG    = 1L ; ensure 4 bytes
   IBCND     = 0L ; ensure 4 bytes
   ONLYFL    = 0L ; ensure 4 bytes
   PRNT =  [0l,0,0,0,0] ; ensure 5 element array of 4 bytes
   PLANK     = 0L ; ensure 4 bytes    
   LAMBER    = 0L ; ensure 4 bytes
   DELTAMPLUS = 0L ; ensure 4 bytes  
   DO_PSEUDO_SPHERE = 0L ; ensure 4 bytes    
   WVNMLO = 0.0 
   WVNMHI = 0.0
   UMU0 = 0.0 
   PHI0 = 0.0 
   FBEAM = 0.0
   FISOT = 0.0 
   ALBEDO = 0.0 
   BTEMP = 0.0 
   TTEMP = 0.0 
   TEMIS = 0.0
   EARTH_RADIUS=6371.0 
   ACCUR = 0.0
   HEADER = string(replicate(' ',127),format='(127a)') ; Header string placeholder
end