; procedure load_inststr
;
; Read an intrument definition file.
;
; INPUT ARGUMENTS:
; file (string) input file.
;
; INPUT KEYWORDS:
; Channelid defines which channels to load otherwise load all.
;
; OUTPUT ARGUMENTS:
; inststr (structure) Structure with instrument information.
;
; HISTORY:
; 20/07/20, RGG: Created

pro load_inststr, file, inststr, RequestedChannelID = RequestedChannelID

  instrument_filename = FILE_BASENAME(file)
; Read in the instrument description
  Lines = strarr(File_lines(file))
  openr, lun, file, /get_lun
    readf, lun, lines
  close, lun
  free_lun, lun  
; Remone Leading and trailing spaces
  Lines = strtrim(Lines,2)  
; Remove commented lines (those starting with "*" or "#".
  Q = WHERE(Lines NE '' AND STRMID(Lines,0,1) NE '*' AND STRMID(lines,0,1) NE '#', Count )
  Lines  = Lines[Q]      
; first find number of channels
; Now deconstruct line and set output values
  LHS = strarr(Count)
  RHS = strarr(Count)
  Variable = strarr(Count)
  FOR I = 0, Count - 1 do begin  
    Result = strsplit(Lines[I],'=',/extract) 
    LHS[I] = Result[0]
    RHS[I] = Result[1]
  ENDFOR
  LSB = strpos(LHS,'[')
  RSB = strpos(LHS,']')
  Q = where(LSB  eq -1, Singles, Complement = NQ, NComplement = Multiples)
 
  FOR I = 0, Singles - 1 do begin    
    Case strlowcase(strcompress(LHS[Q[I]],/remove_all)) of
     'platform'                 : platform            = strlowcase(strtrim(RHS[Q[I]],2))
     'instrument'               : instrument          = strlowcase(strtrim(RHS[Q[I]],2))
     'instrumentversion'        : instrument_version  = strlowcase(strtrim(RHS[Q[I]],2))
     'view'                     : view                = fix(RHS[Q[I]])
     'availablechannels'        : available_channelid = fix(strsplit(RHS[Q[I]],' ',/extract))
     'solarchannels'            : solar_channelid     = fix(strsplit(RHS[Q[I]],' ',/extract))
     'thermalchannels'          : thermal_ChannelID   = fix(strsplit(RHS[Q[I]],' ',/extract))
     'maximumsatellitezenith'   : max_Sat_Zenith      = float(RHS[Q[I]])
    ELSE: message,'Warning: Unknown variable in instrument definition file: '+file     
    ENDCASE 
  ENDFOR
   
  number_of_nadir_channels = N_elements(available_channelid)
  SRF_File = strarr(number_of_nadir_channels)  
;
; TEMPORARY FOR BACK COMPATIBILITY
   oldF0   = fltarr(number_of_nadir_channels)
   oldF1   = fltarr(number_of_nadir_channels)
   oldNEFR = fltarr(number_of_nadir_channels)
   oldWVN  = fltarr(number_of_nadir_channels)
   oldb1   = fltarr(number_of_nadir_channels)
   oldb2   = fltarr(number_of_nadir_channels)
   oldt1   = fltarr(number_of_nadir_channels)
   oldt2   = fltarr(number_of_nadir_channels)
   oldnebt = fltarr(number_of_nadir_channels)
; *******************************
  snr   = fltarr(number_of_nadir_channels)
  rgu   = fltarr(number_of_nadir_channels)
  rou   = fltarr(number_of_nadir_channels)
  rua   = fltarr(number_of_nadir_channels)
  rub   = fltarr(number_of_nadir_channels)
  ruc   = fltarr(number_of_nadir_channels)
  refbt = fltarr(number_of_nadir_channels)
  nedt  = fltarr(number_of_nadir_channels)
  Solar_Channel_Flag   = replicate(0,number_of_nadir_channels)
  Thermal_Channel_Flag = replicate(0,number_of_nadir_channels)
  FOR I = 0,N_Elements(solar_channelid) - 1 do begin
    J = Where(solar_channelid(I) Eq available_channelid, Count)
    IF (Count Eq -1) then $        
      Stop, 'Solar channel not included in available channels'        $
    Else $     
      Solar_Channel_Flag[J] = 1 
  ENDFOR
  FOR I = 0,N_Elements(Thermal_ChannelID) - 1 do begin
    J = where (Thermal_ChannelID(I) eq available_channelid, Count)
    IF (Count Eq -1) then $
      Stop, 'Thermal channel not included in available channels' $
    Else $
      Thermal_Channel_Flag[J] = 1
  ENDFOR

  FOR i = 0, multiples - 1 do begin
    channelid  = string(fix(strmid(LHS[NQ[I]],LSB[NQ[I]]+1,RSB[NQ[I]]-LSB[NQ[I]]-1)),FORMAT='(I2.2)')
    J = where (channelid eq available_channelid, Count)
    IF (Count Eq -1) then $
      Stop, 'channel not included in available channels' $
    Else begin
      Case strlowcase(strcompress(strmid(LHS[NQ[I]],0,LSB[NQ[I]]),/remove_all)) of
        'srf'   : SRF_File[J] = strtrim(RHS[NQ[I]],2)
;       TEMPORARY FOR BACK COMPATIBILITY
        'oldf0'   : oldf0[J]    = strtrim(RHS[NQ[I]],2)
        'oldf1'   : oldf1[J]    = strtrim(RHS[NQ[I]],2)
        'oldnefr' : oldnefr[J] = strtrim(RHS[NQ[I]],2) 
        'oldwvn'  : oldwvn[J]   = strtrim(RHS[NQ[I]],2) 
        'oldb1'   : oldb1[J]   = strtrim(RHS[NQ[I]],2) 
        'oldb2'   : oldb2[J]   = strtrim(RHS[NQ[I]],2) 
        'oldt1'   : oldt1[J]   = strtrim(RHS[NQ[I]],2) 
        'oldt2'   : oldt2[J]   = strtrim(RHS[NQ[I]],2) 
        'oldnebt' : oldnebt[J] = strtrim(RHS[NQ[I]],2) 
;       *******************************

        'snr'  : snr[J]   = strtrim(RHS[NQ[I]],2)
        'rgu'  : rgu[J]   = strtrim(RHS[NQ[I]],2)
        'rou'  : rou[J]   = strtrim(RHS[NQ[I]],2)
        'rua'  : rua[J]   = strtrim(RHS[NQ[I]],2)
        'rub'  : rub[J]   = strtrim(RHS[NQ[I]],2)
        'ruc'  : ruc[J]   = strtrim(RHS[NQ[I]],2)
        'refbt': refbt[J] = strtrim(RHS[NQ[I]],2)
        'nedt' : nedt[J]  = strtrim(RHS[NQ[I]],2)
      Else: print,'Warning: Unknown variable in instrument definition file: '+   strlowcase(strcompress(LHS[NQ[I]],/remove_all)) 
      EndCase 
    Endelse
  ENDFOR
    
; If the RequestedChannelID keyword has been set, then build the instrument structure for  the selected channels.
; Keeping only those channels contained in the keyword.
  IF KEYWORD_SET(RequestedChannelID) gt 0 then begin
      match = REPLICATE(0, number_of_nadir_channels)
      FOR i=0,N_ELEMENTS(Requestedchannelid)-1 do begin
        matchch = WHERE(Requestedchannelid[i] eq available_channelid,count)
        Case Count of
          0: message,/info, 'Warning: channel '+strtrim(Requestedchannelid[i])+ ' not found in '+file 
          1: match[matchch[0]] = 1
          ELSE: message,/info, 'Warning: channel '+strtrim(Requestedchannelid[i])+ ' found multiple times in '+file
        EndCase
      ENDFOR
      matchi = WHERE(match,number_of_nadir_channels)
      available_channelid  = available_channelid[matchi] 
      solar_channel_flag   = solar_channel_flag[matchi]       
      thermal_channel_flag = thermal_channel_flag[matchi]  
      srf_file             = srf_file[matchi]
;     TEMPORARY FOR BACK COMPATIBILITY
      oldf0    =  oldf0[matchi]      
      oldf1    =  oldf1[matchi]      
      oldnefr = oldnefr[matchi]
      oldwvn  =  oldwvn[matchi] 
      oldb1   =   oldb1[matchi] 
      oldb2   =   oldb2[matchi] 
      oldt1   =   oldt1[matchi]
      oldt2   =   oldt2[matchi]
      oldnebt = oldnebt[matchi] 
;     *******************************
      snr   = snr[matchi]
      rgu   = rgu[matchi]
      rou   = rou[matchi]
      rua   = rua[matchi]
      rub   = rub[matchi]
      ruc   = ruc[matchi]
      refbt = refbt[matchi] 
      nedt  = nedt[matchi] 
  Endif     
  
  mixed_channel_flag = solar_channel_flag and thermal_channel_flag

;  If dual view than nadir results are replicated for forward view
   If (View Gt 0) Then $
     Number_of_channels = 2 * number_of_nadir_channels $
   Else $
     Number_of_channels = number_of_nadir_channels   
        
  inststr = {instrument_filename: instrument_filename      , $
                        platform: platform                 , $
                      instrument: instrument               , $
              instrument_version: instrument_version       , $
                  max_sat_zenith: max_sat_zenith           , $
                            view: view                     , $
              number_of_channels: number_of_channels       , $        
        number_of_nadir_channels: number_of_nadir_channels , $            
                       channelid: available_channelid      , $
              solar_channel_flag: solar_channel_flag       , $
              mixed_channel_flag: mixed_channel_flag       , $       
            thermal_channel_flag: thermal_channel_flag     , $       
                        srf_file: srf_file                 , $                 
;     TEMPORARY FOR BACK COMPATIBILITY
                           oldf0: oldf0                    , $
                           oldf1: oldf1                    , $
                         oldnefr: oldnefr                  , $
                          oldwvn: oldwvn                   , $
                           oldb1: oldb1                    , $
                           oldb2: oldb2                    , $
                           oldt1: oldt1                    , $ 
                           oldt2: oldt2                    , $
                         oldnebt: oldnebt                  , $
;     *******************************                  
                             snr: snr                      , $
                             rgu: rgu                      , $               
                             rou: rou                      , $               
                             rua: rua                      , $               
                             rub: rub                      , $               
                             ruc: ruc                      , $               
                           refbt: refbt                    , $               
                            nedt: nedt}              
end


