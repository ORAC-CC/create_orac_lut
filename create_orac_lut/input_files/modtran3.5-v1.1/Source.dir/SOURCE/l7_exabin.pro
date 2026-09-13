;==============================================================================
;+
; L7_EXABIN
; 
; IDL version of modtran routine to load aerosol parameters.
;C                                                                       EXA 0004
;C      LOADS EXTINCTION, ABSORPTION AND ASYMMETRY COEFFICIENTS          EXA 0005
;C      FOR THE FOUR AEROSOL ALTITUDE REGIONS                            EXA 0006
;C                                                                       EXA 0007
;C      MODIFIED FOR ASYMMETRY - JAN 1986 (A.E.R. INC.)                  EXA 0008
;C                                                                       EXA 0009
;
; PARAMETERS
;	EXDTA	Extinction data as read by l7_extdta
;    On return...
;	EXTC	Aerosol extinction coefficient
;	ABSC	absorption coefficient
;	ASYM	Asymmetry parameters
;
; KEYWORDS
;	ISEASN
;	IVULCN
;	ICSTL
;	ICLD	All as modtran card 2 parameters
;	IVSA
;	VIS
;	WSS
;	WHH
;	RAINRT
;
; R.S. 13/11/97
;-
;==============================================================================
pro l7_exabin,EXDTA,CH,EXTC,ABSC,ASYM,$
	ihaze=IHAZE,iseasn=ISEASN,ivulvn=IVULCN,$
	icstl=ICSTL,icld=ICLD,ivsa=IVSA,vis=VIS,wss=WSS,whh=WHH,rainrt=RAINRT
	if not keyword_set(ihaze) then ihaze=0
	if not keyword_set(iseasn) then iseasn=1
        if not keyword_set(ivulcn) then ivulcn=1
	if not keyword_set(icstl) then icstl=3
	if not keyword_set(icld) then icld=0
	if not keyword_set(ivsa) then ivsa=0
	if not keyword_set(vis) then vis=0
	if not keyword_set(wss) then wss=0
	if not keyword_set(whh) then whh=0
	if not keyword_set(rainrt) then rainrt=0.
	nwavln=n_elements(exdta.vx0)
        IF IVSA EQ 1 AND IHAZE EQ 3   then message,'Cannot cope with given conditions !'
        ICH=intarr(4)
        ICH(0)=IHAZE
        ICH(1)=6 
        ICH(2)=9+IVULCN
        ICH(3)=18
        IF ICH(0) LE 0 then ICH(0)=1
        IF ICH(2) LE 9 then ICH(2)=10
        IF ICLD EQ 11 THEN begin
          ICH(3)=ICH(2)
          ICH(2)=ICH(1)
          ICH(1)=ICLD
        ENDIF
;      COMMON /CARD2/ IHAZE,ISEASN,IVULCN,ICSTL,ICLD,IVSA,VIS,WSS,WHH,   
;     1    RAINRT                                                        
;      COMMON /CARD2D/ IREG(4),ALTB(4),IREGC(4)                          
;C                                                                       
;      INTEGER KPOINT                                                    
;      REAL WAVLEN,EXTC,ABSC,ASYM,AWCCON,TX,W,RELHUM,PATM,TBBY,WPATH     
;      COMMON/BASE/WAVLEN(MXWVLN),EXTC(NAER,MXWVLN),ABSC(NAER,MXWVLN),   
;     1  ASYM(NAER,MXWVLN),AWCCON(NAER),KPOINT(NSPC),TX(KMAX),W(KMAX),   
;     2  RELHUM(LAYDIM),PATM(LAYTHR),TBBY(LAYTHR),WPATH(LAYTHR,KMAX)     
      RHZONE=[0.,70.,80.,99.]                                       
      ELWCR=[3.517E-04,3.740E-04,4.439E-04,9.529E-04]               
      ELWCM=[4.675E-04,6.543E-04,1.166E-03,3.154E-03]               
      ELWCU=[3.102E-04,3.802E-04,4.463E-04,9.745E-04]               
      ELWCT=[1.735E-04,1.820E-04,2.020E-04,2.408E-04]               
      AFLWC=1.295E-02 & RFLWC=1.804E-03 & CULWC=7.683E-03           
      ASLWC=4.509E-03 & STLWC=5.272E-03 & SCLWC=4.177E-03           
      SNLWC=7.518E-03 & BSLWC=1.567E-04 & FVLWC=5.922E-04           
      AVLWC=1.675E-04 & MDLWC=4.775E-04                            
      awccon=fltarr(4)
;C                                                                       
;C     "NWAVLN" VALUES CALCULATED IN AEREXT USING ROUTINE GAMFOG         
      NEW=NWAVLN-1                                                       
;C                                                                       
      for M=0 ,3 do begin
      	ITA=ICH(M)                                                        
      	ITC=ICH(M)-7                                                      
      	ITAS = ITA                                                        
        WRH=W(15)                                                         
      	if ICH(M) eq 6 AND M NE 1 then WRH=70.                               
;C     THIS CODING  DOES NOT ALLOW TROP RH DEPENDENT  ABOVE EH(7,I)      
;C     DEFAULTS TO TROPOSPHERIC AT 70. PERCENT                           
      	for I=1,3 do if WRH lt RHZONE(I) then goto,exabin15                                    
      	I=4                                                               
exabin15:
        II=I-1                                                            
      	if WRH GT 0.0 AND WRH lt 99. then X=ALOG(100.0-WRH)                    
      	X1=ALOG(100.0-RHZONE(II))                                         
      	X2=ALOG(100.0-RHZONE(I))                                          
      	if WRH GE 99.0 then X=X2                                             
      	if WRH LE 0.0 then X=X1                                              
        for N=0,NEW-1 do begin
      		ITA = ITAS                                                        
      		if ITA eq 3 AND M EQ 0 then goto,exabin18                                
      		ABSC(M,N)=0.                                                      
      		EXTC(M,N)=0.                                                      
      		ASYM(M,N)=0.                                                      
      		if ITA GT 6 then goto,exabin45                                             
      		if ITA LE 0 then goto,exabin80                                             
exabin18:
      		if N GE 40  AND  ITA eq 3 then ITA = 4                                
;C     RH DEPENDENT AEROSOLS                                             
      		if ITA eq 1 then goto,exabin20
      		if ITA eq 2 then goto,exabin20
      		if ITA eq 3 then goto,exabin22
      		if ITA eq 4 then goto,exabin25
      		if ITA eq 5 then goto,exabin30
      		if ITA eq 6 then goto,exabin35
exabin20:
	 	Y2=ALOG(RUREXT(N,I))                                              
      		Y1=ALOG(RUREXT(N,II))                                             
      		Z2=ALOG(RURABS(N,I))                                              
      		Z1=ALOG(RURABS(N,II))                                             
      		A2=ALOG(RURSYM(N,I))                                              
      		A1=ALOG(RURSYM(N,II))                                             
      		E2=ALOG(ELWCR(I))                                                 
      		E1=ALOG(ELWCR(II))                                                
      		goto,exabin40                                                          
exabin22:
	  	if M GT 1 then goto,exabin25                                               
      		A2=ALOG(OCNSYM(N,I))                                              
      		A1=ALOG(OCNSYM(N,II))                                             
      		A=A1+(A2-A1)*(X-X1)/(X2-X1)                                       
      		ASYM(M,N)=EXP(A)                                                  
      		E2=ALOG(ELWCM(I))                                                 
      		E1=ALOG(ELWCM(II))                                                
;C                                                                       
;C     NAVY MARITIME AEROSOL CHANGES TO MARINE IN MICROWAVE              
;C     NO NEED TO DEFINE EQUIVALENT WATER                                
;C                                                                       
      		goto,exabin80                                                          
exabin25:
	 	Y2=ALOG(OCNEXT(N,I))                                              
      		Y1=ALOG(OCNEXT(N,II))                                             
      		Z2=ALOG(OCNABS(N,I))                                              
      		Z1=ALOG(OCNABS(N,II))                                             
      		A2=ALOG(OCNSYM(N,I))                                              
      		A1=ALOG(OCNSYM(N,II))                                             
      		E2=ALOG(ELWCM(I))                                                 
      		E1=ALOG(ELWCM(II))                                                
      		goto,exabin40                                                          
exabin30:
		Y2=ALOG(URBEXT(N,I))                                              
      		Y1=ALOG(URBEXT(N,II))                                             
      		Z2=ALOG(URBABS(N,I))                                              
      		Z1=ALOG(URBABS(N,II))                                             
      		A2=ALOG(URBSYM(N,I))                                              
      		A1=ALOG(URBSYM(N,II))                                             
      		E2=ALOG(ELWCU(I))                                                 
      		E1=ALOG(ELWCU(II))                                                
      		goto,exabin40                                                          
exabin35:
		Y2=ALOG(TROEXT(N,I))                                              
      		Y1=ALOG(TROEXT(N,II))                                             
      		Z2=ALOG(TROABS(N,I))                                              
      		Z1=ALOG(TROABS(N,II))                                             
      		A2=ALOG(TROSYM(N,I))                                              
      		A1=ALOG(TROSYM(N,II))                                             
      		E2=ALOG(ELWCT(I))                                                 
      		E1=ALOG(ELWCT(II))                                                
exabin40:
		Y=Y1+(Y2-Y1)*(X-X1)/(X2-X1)                                       
      		ZK=Z1+(Z2-Z1)*(X-X1)/(X2-X1)                                      
      		A=A1+(A2-A1)*(X-X1)/(X2-X1)                                       
      		ABSC(M,N)=EXP(ZK)                                                 
      		EXTC(M,N)=EXP(Y)                                                  
      		ASYM(M,N)=EXP(A)                                                  
      		if N eq 0 then EC=E1+(E2-E1)*(X-X1)/(X2-X1)                           
      		if N eq 0 then AWCCON(M)=EXP(EC)                                      
      		goto,exabin80                                                          
exabin45:
		if ITA GT 19 then goto,exabin75                                           
      		if ITC lt 1 then goto,exabin80                                            
      		if ITC eq 1 then goto,exabin50
      		if ITC eq 2 then goto,exabin55
      		if ITC eq 3 then goto,exabin80
      		if ITC eq 4 then goto,exabin60
      		if ITC eq 5 then goto,exabin65
      		if ITC eq 6 then goto,exabin70
      		if ITC eq 7 then goto,exabin65
      		if ITC eq 8 then goto,exabin70
      		if ITC eq 9 then goto,exabin60
      		if ITC eq 10 then goto,exabin60
      		if ITC eq 11 then goto,exabin70
      		if ITC eq 12 then goto,exabin75
exabin50:
 		ABSC(M,N)=FG1ABS(N)                                               
      		EXTC(M,N)=FG1EXT(N)                                               
      		ASYM(M,N)=FG1SYM(N)                                               
      		if N eq 0 then AWCCON(M)=AFLWC                                        
      		goto,exabin80                                                          
exabin55:
 		ABSC(M,N)=FG2ABS(N)                                               
      		EXTC(M,N)=FG2EXT(N)                                               
      		ASYM(M,N)=FG2SYM(N)                                               
      		if N eq 0 then AWCCON(M)=RFLWC                                        
      		goto,exabin80                                                          
exabin60:
 		ABSC(M,N)=BSTABS(N)                                               
      		EXTC(M,N)=BSTEXT(N)                                               
      		ASYM(M,N)=BSTSYM(N)                                               
      		if N eq 0 then AWCCON(M)=BSLWC                                        
      		goto,exabin80                                                          
exabin65:
 		ABSC(M,N)=AVOABS(N)                                               
      		EXTC(M,N)=AVOEXT(N)                                               
      		ASYM(M,N)=AVOSYM(N)                                               
      		if N eq 0 then AWCCON(M)=AVLWC                                        
      		goto,exabin80                                                          
exabin70:
 		ABSC(M,N)=FVOABS(N)                                               
      		EXTC(M,N)=FVOEXT(N)                                               
      		ASYM(M,N)=FVOSYM(N)                                               
      		if N eq 0 then AWCCON(M)=FVLWC                                        
      		goto,exabin80                                                          
exabin75:
 		ABSC(M,N)=DMEABS(N)                                               
      		EXTC(M,N)=DMEEXT(N)                                               
      		ASYM(M,N)=DMESYM(N)                                               
      		if N eq 0 then AWCCON(M)=MDLWC                                        
exabin80:
      	endfor
      	goto,exabin185                                                         
exabin100:
;CCC                                                                     
;CCC       SECTION TO LOAD EXTINCTION AND ABSORPTION COEFFICIENTS        
;CCC       FOR CLOUD AND OR RAIN MODELS                                  
;CCC                                                                     
      	for N=0,NEW-1 do begin
      		ABSC(M,N)=0.0                                                     
      		EXTC(M,N)=0.0                                                     
      		ASYM(M,N)=0.0                                                     
      		IC=ICLD                                                           
      		if ic eq 1 then goto,exabin125
      		if ic eq 2 then goto,exabin130
      		if ic eq 3 then goto,exabin135
      		if ic eq 4 then goto,exabin140
      		if ic eq 5 then goto,exabin145
      		if ic eq 6 then goto,exabin135
      		if ic eq 7 then goto,exabin145
      		if ic eq 8 then goto,exabin145
      		if ic eq 9 then goto,exabin125
      		if ic eq 10 then goto,exabin125
      		if ic eq 11 then goto,exabin125
exabin125:
      		ABSC(M,N)=CCUABS(N)                                               
      		EXTC(M,N)=CCUEXT(N)                                               
      		ASYM(M,N)=CCUSYM(N)                                               
      		if N eq 0 then AWCCON(M)=CULWC                                        
      		goto,exabin150                                                         
exabin130:
      		ABSC(M,N)=CALABS(N)                                               
      		EXTC(M,N)=CALEXT(N)                                               
      		ASYM(M,N)=CALSYM(N)                                               
      		if N eq 0 then AWCCON(M)=ASLWC                                        
      		goto,exabin150                                                         
exabin135:
      		ABSC(M,N)=CSTABS(N)                                               
      		EXTC(M,N)=CSTEXT(N)                                               
      		ASYM(M,N)=CSTSYM(N)                                               
      		if N eq 0 then AWCCON(M)=STLWC                                        
      		goto,exabin150                                                         
exabin140:
      		ABSC(M,N)=CSCABS(N)                                               
      		EXTC(M,N)=CSCEXT(N)                                               
      		ASYM(M,N)=CSCSYM(N)                                               
      		if N eq 0 then AWCCON(M)=SCLWC                                        
      		goto,exabin150                                                         
exabin145:
      		ABSC(M,N)=CNIABS(N)                                               
      		EXTC(M,N)=CNIEXT(N)                                               
      		ASYM(M,N)=CNISYM(N)                                               
      		if N eq 0 then  AWCCON(M)=SNLWC                                        
exabin150:
	endfor
exabin185:
      endfor
      for N=0,NWAVLN-1 do begin
      	ABSC(5,N)=0.                                                      
      	EXTC(5,N)=0.                                                      
      	ASYM(5,N)=0.                                                      
      	AWCCON(5)=0.                                                      
      	if ICLD  eq  18 then begin                                             
           	EXTC(5,N)= CI64XT(N)                                         
           	ABSC(5,N)= CI64AB(N)                                         
           	ASYM(5,N)= CI64G(N)                                          
           	AWCCON(5)=3.446E-3                                           
      	ENDIF                                                             
      	if ICLD  eq  19 then begin                                             
           	EXTC(5,N)= CIR4XT(N)                                         
           	ABSC(5,N)= CIR4AB(N)                                         
           	ASYM(5,N)= CIR4G(N)                                          
           	AWCCON(5)=5.811E-2                                           
      	ENDIF                                                             
      endfor
END                                                               
