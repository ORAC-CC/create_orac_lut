pro read_modtran,filename, wvn1, nrows, wvn, trans
; wvn wavenumber
; trans transmission
; (12 july 2005) Elisa: gas trans = total trans/molecular scattering

  header=strarr(12)             ; 12 lines of Header in file.
  dataline=' '                  ; last line is -9999
  data = fltarr(22,nrows)
  openr,fID,filename,/GET_LUN
    readf,fID,header
    wvnval=-1
    while wvnval lt floor(wvn1)-1 do begin
       readf,fID,dataline
       vals = strsplit(dataline,' ')
       wvnval = long(strmid(dataline,vals[0],vals[1]-vals[0]))
    endwhile   
    readf,fID,data
  close,fid
  free_lun,fID
  wvn      = data(0,*)
  tottrans = data(1,*)
  molsca   = data(8,*)
  trans=  ((tottrans/molsca <1.0) >0.0)
END
