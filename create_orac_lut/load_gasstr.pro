
; procedure load
;
; Read a gas optical depth profile file.
;

; INPUT ARGUMENTS:
; file (string) Path to file.
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; gasstr (structure) Structure with gas optical depth profile information.
;
; HISTORY:
; 21/06/13, G Thomas: Original version.

pro read_gasstr, Filename, Onegasstr
  line = ''
  openr,lun,filename,/get_lun
    for i = 0,1 do readf,lun, line
    for i = 0,3 do begin
      readf,lun, line
      case strlowcase(strtrim(strmid(line,1),2)) of
        'atmosphere': begin
                        readf,lun, line
                        atmosphere = strtrim(line,2)
                      end
        'instrument': begin
               readf,lun, line
               instrument = strtrim(line,2)
            end
        'channelid': begin
               readf,lun, line
               channelid = strtrim(line,2)
            end
        'nlevels': begin
               readf,lun, line
               nlevels = long(line)
            end
      endcase    
    endfor
    Values = Fltarr(2,Nlevels)
    readf,lun, Values
  free_lun, lun
  Onegasstr = {atmosphere: atmosphere,  $
               instrument: instrument,  $
                channelid: channelid ,  $
                  nlevels: nlevels ,    $
                   height: reform(values[0,*]), $
                  tau_gas: reform(values[1,*])}
end

pro load_gasstr, atmospheres, platform, instrument, channelid,  gasdir, gasstr
  NChannels = N_Elements(channelid)
; readin first gas profile
  gasfile = 'ModtranGasOpd_A'+atmospheres+'_'+platform+'_'+instrument+'_'+'ch'+string(channelid[0],Format='(I2.2)')+'.gas'   
  read_gasstr,gasdir+gasfile,onegasstr
  gasstr = replicate(onegasstr,NChannels)
  
  For I = 1, NChannels -1 do begin
    channels = string(channelid[i],Format='(I2.2)')   
    gasfile = 'ModtranGasOpd_A'+atmospheres+'_'+platform+'_'+instrument+'_'+'ch'+channels+'.gas' 
    read_gasstr,gasdir+gasfile,onegasstr
    gasstr(i) = onegasstr
  EndFor
end