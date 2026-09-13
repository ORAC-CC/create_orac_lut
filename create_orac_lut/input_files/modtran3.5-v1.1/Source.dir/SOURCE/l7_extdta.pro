function l7_extdta
;+
;      BLOCK DATA EXTDTA                            
;      INCLUDE 'PARAM.LST'                          
;C>    BLOCK DATA                                   
;CCC                                                
;CCC   ALTITUDE REGIONS FOR AEROSOL EXTINCTION COEFFICIENTS
;CCC                                                
;CCC                                                
;CCC         0-2KM                                  
;CCC           RUREXT=RURAL EXTINCTION   RURABS=RURAL ABSORPTION
;CCC           RURSYM=RURAL ASYMMETRY FACTORS       
;CCC           URBEXT=URBAN EXTINCTION   URBABS=URBAN ABSORPTION
;CCC           URBSYM=URBAN ASYMMETRY FACTORS       
;CCC           OCNEXT=MARITIME EXTINCTION  OCNABS=MARITIME ABSORPTION
;CCC           OCNSYM=MARITIME ASYMMETRY FACTORS    
;CCC           TROEXT=TROPSPHER EXTINCTION  TROABS=TROPOSPHER ABSORPTION
;CCC           TROSYM=TROPSPHERIC ASYMMETRY FACTORS 
;CCC           FG1EXT=FOG1 .2KM VIS EXTINCTION  FG1ABS=FOG1 ABSORPTION
;CCC           FG1SYM=FOG1 ASYMMETRY FACTORS        
;CCC           FG2EXT=FOG2 .5KM VIS EXTINCTION  FG2ABS=FOG2 ABSORPTION
;CCC           FG2SYM=FOG2 ASYMMETRY FACTORS        
;CCC         >2-10KM                                
;CCC           TROEXT=TROPOSPHER EXTINCTION  TROABS=TROPOSPHER ABSORPTION
;CCC           TROSYM=TROPOSPHERIC ASYMMETRY FACTORS
;CCC         >10-30KM                               
;CCC           BSTEXT=BACKGROUND STRATOSPHERIC EXTINCTION
;CCC           BSTABS=BACKGROUND STRATOSPHERIC ABSORPTION
;CCC           BSTSYM=BACKGROUND STRATOSPHERIC ASYMMETRY FACTORS
;CCC           AVOEXT=AGED VOLCANIC EXTINCTION      
;CCC           AVOABS=AGED VOLCANIC ABSORPTION      
;CCC           AVOSYM=AGED VOLCANIC ASYMMETRY FACTORS
;CCC           FVOEXT=FRESH VOLCANIC EXTINCTION     
;CCC           FVOABS=FRESH VOLCANIC ABSORPTION     
;CCC           FVOSYM=FRESH VOLCANIC ASYMMETRY FACTORS
;CCC         >30-100KM                              
;CCC           DMEEXT=METEORIC DUST EXTINCTION      
;CCC           DMEABS=METEORIC DUST ABSORPTION      
;CCC           DMESYM=METEORIC DUST ASYMMETRY FACTORS
;C                                                  
;C     AEROSOL EXTINCTION AND ABSORPTION DATA       
;C                                                  
;C     MODIFIED TO INCLUDE ASYMMETRY DATA - JAN 1986 (A.E.R. INC.)
;-                                                  
;      COMMON/EXTD/VX0(NWAVLN),                     
;     1  RUREXT(NWAVLN,4),RURABS(NWAVLN,4),RURSYM(NWAVLN,4),
;     2  URBEXT(NWAVLN,4),URBABS(NWAVLN,4),URBSYM(NWAVLN,4),
;     3  OCNEXT(NWAVLN,4),OCNABS(NWAVLN,4),OCNSYM(NWAVLN,4),
;     4  TROEXT(NWAVLN,4),TROABS(NWAVLN,4),TROSYM(NWAVLN,4),
;     5  FG1EXT(NWAVLN),FG1ABS(NWAVLN),FG1SYM(NWAVLN),
;     6  FG2EXT(NWAVLN),FG2ABS(NWAVLN),FG2SYM(NWAVLN),
;     7  BSTEXT(NWAVLN),BSTABS(NWAVLN),BSTSYM(NWAVLN),
;     8  AVOEXT(NWAVLN),AVOABS(NWAVLN),AVOSYM(NWAVLN),
;     9  FVOEXT(NWAVLN),FVOABS(NWAVLN),FVOSYM(NWAVLN),
;     &  DMEEXT(NWAVLN),DMEABS(NWAVLN),DMESYM(NWAVLN),
;     1  CCUEXT(NWAVLN),CCUABS(NWAVLN),CCUSYM(NWAVLN),
;     2  CALEXT(NWAVLN),CALABS(NWAVLN),CALSYM(NWAVLN),
;     3  CSTEXT(NWAVLN),CSTABS(NWAVLN),CSTSYM(NWAVLN),
;     4  CSCEXT(NWAVLN),CSCABS(NWAVLN),CSCSYM(NWAVLN),
;     5  CNIEXT(NWAVLN),CNIABS(NWAVLN),CNISYM(NWAVLN)
;C                                                  
;C         CI64--    STANDARD  CIRRUS  CLOUD  MODEL 
;C              ICE 64 MICRON MODE RADIUS CIRRUS CLOUD MODEL
;C                                                  
;C         CIR4--    OPTICALLY  THIN  CIRRUS  MODEL 
;C              ICE  4 MICRON MODE RADIUS CIRRUS CLOUD MODEL
;C                                                  
;      COMMON/CIRR/CI64XT(NWAVLN),CI64AB(NWAVLN),CI64G(NWAVLN),
;     1            CIR4XT(NWAVLN),CIR4AB(NWAVLN),CIR4G(NWAVLN)

      NWAVLN=47         
      VX0=fltarr(NWAVLN)
      RUREXT=fltarr(NWAVLN,4) & RURABS=fltarr(NWAVLN,4) & RURSYM=fltarr(NWAVLN,4)
      URBEXT=fltarr(NWAVLN,4) & URBABS=fltarr(NWAVLN,4) & URBSYM=fltarr(NWAVLN,4)
      OCNEXT=fltarr(NWAVLN,4) & OCNABS=fltarr(NWAVLN,4) & OCNSYM=fltarr(NWAVLN,4)
      TROEXT=fltarr(NWAVLN,4) & TROABS=fltarr(NWAVLN,4) & TROSYM=fltarr(NWAVLN,4)
      FG1EXT=fltarr(NWAVLN) &  FG1ABS=fltarr(NWAVLN) & FG1SYM=fltarr(NWAVLN)
      FG2EXT=fltarr(NWAVLN) & FG2ABS=fltarr(NWAVLN) & FG2SYM=fltarr(NWAVLN)
      BSTEXT=fltarr(NWAVLN) & BSTABS=fltarr(NWAVLN) & BSTSYM=fltarr(NWAVLN)
      AVOEXT=fltarr(NWAVLN) & AVOABS=fltarr(NWAVLN) & AVOSYM=fltarr(NWAVLN)
      FVOEXT=fltarr(NWAVLN) & FVOABS=fltarr(NWAVLN) & FVOSYM=fltarr(NWAVLN)
      DMEEXT=fltarr(NWAVLN) & DMEABS=fltarr(NWAVLN) & DMESYM=fltarr(NWAVLN)
      CCUEXT=fltarr(NWAVLN) & CCUABS=fltarr(NWAVLN) & CCUSYM=fltarr(NWAVLN)
      CALEXT=fltarr(NWAVLN) & CALABS=fltarr(NWAVLN) & CALSYM=fltarr(NWAVLN)
      CSTEXT=fltarr(NWAVLN) & CSTABS=fltarr(NWAVLN) & CSTSYM=fltarr(NWAVLN)
      CSCEXT=fltarr(NWAVLN) & CSCABS=fltarr(NWAVLN) & CSCSYM=fltarr(NWAVLN)
      CNIEXT=fltarr(NWAVLN) & CNIABS=fltarr(NWAVLN) & CNISYM=fltarr(NWAVLN)

	get_lun,lun
	openr,lun,'l7_extdta.dat'
	on_ioerror,close_file
	line=';'
	while strmid(line,0,1) eq ';' do begin
		point_lun,-lun,pos
		readf,lun,line
	endwhile
	point_lun,lun,pos
	readf,lun,vx0
	readf,lun,line
	readf,lun,RUREXT
	readf,lun,line
	readf,lun,RURABS
	readf,lun,line
	readf,lun,RURSYM
	readf,lun,line
	readf,lun,URBEXT
	readf,lun,line
	readf,lun,URBABS
	readf,lun,line
	readf,lun,URBSYM
	readf,lun,line
	readf,lun,OCNEXT
	readf,lun,line
	readf,lun,OCNABS
	readf,lun,line
	readf,lun,OCNSYM
	readf,lun,line
	readf,lun,TROEXT
	readf,lun,line
	readf,lun,TROABS
	readf,lun,line
	readf,lun,TROSYM
	readf,lun,line
	readf,lun,FG1EXT
	readf,lun,line
	readf,lun,FG1ABS
	readf,lun,line
	readf,lun,FG1SYM
	readf,lun,line
	readf,lun,FG2EXT
	readf,lun,line
	readf,lun,FG2ABS
	readf,lun,line
	readf,lun,FG2SYM
	readf,lun,line
	readf,lun,BSTEXT
	readf,lun,line
	readf,lun,BSTABS
	readf,lun,line
	readf,lun,BSTSYM
	readf,lun,line
	readf,lun,AVOEXT
	readf,lun,line
	readf,lun,AVOABS
	readf,lun,line
	readf,lun,AVOSYM
	readf,lun,line
	readf,lun,FVOEXT
	readf,lun,line
	readf,lun,FVOABS
	readf,lun,line
	readf,lun,FVOSYM
	readf,lun,line
	readf,lun,DMEEXT
	readf,lun,line
	readf,lun,DMEABS
	readf,lun,line
	readf,lun,DMESYM
	readf,lun,line
	readf,lun,CCUEXT
	readf,lun,line
	readf,lun,CCUABS
	readf,lun,line
	readf,lun,CCUSYM
	readf,lun,line
	readf,lun,CALEXT
	readf,lun,line
	readf,lun,CALABS
	readf,lun,line
	readf,lun,CALSYM
	readf,lun,line
	readf,lun,CSTEXT
	readf,lun,line
	readf,lun,CSTABS
	readf,lun,line
	readf,lun,CSTSYM
	readf,lun,line
	readf,lun,CSCEXT
	readf,lun,line
	readf,lun,CSCABS
	readf,lun,line
	readf,lun,CSCSYM
	readf,lun,line
	readf,lun,CNIEXT
	readf,lun,line
	readf,lun,CNIABS
	readf,lun,line
	readf,lun,CNISYM
close_file:
	close,lun
	free_lun,lun
;
; form and return structure
;
	stru={vx0:vx0,$
		rurext:RUREXT,$
		rurabs:RURABS,$
		rursym:RURSYM,$
		urbext:URBEXT,$
		urbabs:URBABS,$
		urbsym:URBSYM,$
		ocnext:OCNEXT,$
		ocnabs:OCNABS,$
		ocnsym:OCNSYM,$
		troext:TROEXT,$
		troabs:TROABS,$
		trosym:TROSYM,$
		fg1ext:FG1EXT,$
		fg1abs:FG1ABS,$
		fg1sym:FG1SYM,$
		fg2ext:FG2EXT,$
		fg2abs:FG2ABS,$
		fg2sym:FG2SYM,$
		bstext:BSTEXT,$
		bstabs:BSTABS,$
		bstsym:BSTSYM,$
		avoext:AVOEXT,$
		avoabs:AVOABS,$
		avosym:AVOSYM,$
		fvoext:FVOEXT,$
		fvoabs:FVOABS,$
		fvosym:FVOSYM,$
		dmeext:DMEEXT,$
		dmeabs:DMEABS,$
		dmesym:DMESYM,$
		ccuext:CCUEXT,$
		ccuabs:CCUABS,$
		ccusym:CCUSYM,$
		calext:CALEXT,$
		calabs:CALABS,$
		calsym:CALSYM,$
		cstext:CSTEXT,$
		cstabs:CSTABS,$
		cstsym:CSTSYM,$
		cscext:CSCEXT,$
		cscabs:CSCABS,$
		cscsym:CSCSYM,$
		cniext:CNIEXT,$
		cniabs:CNIABS,$
		cnisym:CNISYM}
	return,stru
end                      
