;+
; This code runs MODTRAN v3.5 to produce gas absorption profiles for a
; climatological atmosphere, across the UV - SW-IR spectral range,
; with a 1 cm-1 spectral resolution. The code writes the "tape" config
; files needed for MODTRAN, before calling MODTRAN itself and storing
; the "tape7" output data.
;
; The code can be run for specific instruments, but it is probably
; more efficient to just run for the entire spectral range and use the
; same tape7 files for all instruments (by using "all" as the
; instrument name).
;
; Once this code has been run, the gas optical depth files can be
; generated using the calc_gas_opd.pro procedure.
; 
; The code assumes it has been called from within the modtran directory.
; modtran_exe specifies the path to the Modtran executable. This is assumed to be modtran_dir/Source.dir/mod3p5.exe
; HISTORY
; 2021/02/15 Don Grainger: Added new atmospheres and kept old (which just has fewer heights) for backwards compatibility.
; 2018/02/09 Gareth Thomas: Modified from the write_tape5 routines (of various versions which had been kicking around the AOPP servers for 15 years or so).
;-

pro create_orac_lut_run_modtran        

;  Below array is old height array           
  oldheights=[ '0.000',  '1.000',  '2.000',  '3.000',  '4.000',  '5.000',  '6.000',  '7.000',  '8.000',  '9.000', $
              '10.000', '11.000', '12.000', '13.000', '14.000', '15.000', '16.000', '17.000', '18.000', '19.000', $
              '20.000', '21.000', '22.000', '23.000', '24.000', '25.000', '30.000', '35.000', '40.000', '45.000', $
              '50.000', '70.000','100.000']
; NewHeights 
  newheights=[ '0.000',  '1.000',  '2.000',  '3.000',  '4.000',   '5.000',   '6.000',   '7.000',   '8.000',   '9.000', $
           '10.000', '11.000', '12.000', '13.000', '14.000',  '15.000',  '16.000',  '17.000',  '18.000',  '19.000', $
           '20.000', '21.000', '22.000', '23.000', '24.000',  '25.000',  '27.500',  '30.000',  '32.500',  '35.000', $
           '37.500', '40.000', '42.500', '45.000', '47.500',  '50.000',  '55.000',  '60.000',  '65.000',  '70.000', $
           '75.000', '80.000', '85.000', '90.000', '95.000', '100.000'];, '105.000', '110.000', '115.000', '120.000']   

  cd, current=modtran_dir
  modtran_exe = 'Source.dir/mod3p5.exe'

 For A = 0, 6 do begin
 
    If (A Eq 0) Then begin
      AID = 2 
      Heights = OldHeights
    EndIf Else Begin
      AID = A
      Heights = NewHeights
    EndElse
   AS = string(A,Format='(I1)')
   AIDS = string(AID,Format='(I1)')

 
    Case AID of
      1: Atmospheres = 'Tropical Atmosphere'
      2: Atmospheres = 'Midlatitude Summer'
      3: Atmospheres = 'Midlatitude Winter'
      4: Atmospheres = 'Subarctic Summer'
      5: Atmospheres = 'Subarctic Winter'
      6: Atmospheres = '1976 US Standard'
    EndCase

    NumberOfHeights = n_elements(heights)

;   Output channel data which can then be read by calcgasopd
    save, filename=modtran_dir+'/output/modtran_height_A'+As+'.sav', atmospheres,heights

    w1 =    '625'
    w2 =  '27000'
    h2 = '500.000'
    
    for j=0,NumberOfHeights -1 do begin        ; height loop stop at 100 km as that is as high as this version of MODTRAN can do
      h1=heights(j)
;     Write the input tape files for Modtran
      openw,fid,'tape5',/get_lun    
        printf,fid,'t   '+AIDS+'    2    0    0    0    0    0    0    0    0    0    0    0    290.000    .00'
        printf,fid,'f   0f   0'
        printf,fid,'    0    0    0    0    0    0      .000      .000      .000      .000      .000'
        printf,fid,'   ',h1,'   ',h2,'    0.000      .000      .000      .000    0'
        printf,fid,'    ',w1,'     ',w2,'         1         0'
        printf,fid,'    0'
        close,fid
      free_lun,fid
      print,'run MODTRAN '+Atmospheres+ ' h1 = '+h1
;     Run Modtran and copy output file
     spawn, modtran_exe
      print,'Copying tape7 to output/tape7_A'+As+'_'+strtrim(fix(h1),2)
      file_copy, 'tape7', modtran_dir+'/output/tape7_'+As+'_'+strtrim(fix(h1),2), /overwrite
;     Tidy up
      file_delete, 'tape5'
      file_delete, 'tape6'
      file_delete, 'tape7'
      file_delete, 'tape8'
    endfor
  endfor
END
