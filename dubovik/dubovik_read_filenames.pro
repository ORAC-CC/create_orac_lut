pro dubovik_read_filenames, file, kext = file_ke, $
                            k11 = file_11, k12  = file_12, $
                            k22 = file_22, k33  = file_33, $
                            k34 = file_34, k44  = file_44


if ~ file_test(file) then message, file + ' not found'


openr, lun, file,/get_lun

readf, lun, n_files

file_11 = strarr(n_files)
file_12 = strarr(n_files)
file_22 = strarr(n_files)
file_33 = strarr(n_files)
file_34 = strarr(n_files)
file_44 = strarr(n_files)
file_ke = strarr(n_files)

tmp = ''

for i = 0, n_files-1 do begin

    readf, lun, tmp, form='(a)'
    file_11[i] = tmp
    readf, lun, tmp, form='(a)'
    file_12[i] = tmp
    readf, lun, tmp, form='(a)'
    file_22[i] = tmp
    readf, lun, tmp, form='(a)'
    file_33[i] = tmp
    readf, lun, tmp, form='(a)'
    file_34[i] = tmp
    readf, lun, tmp, form='(a)'
    file_44[i] = tmp
    readf, lun, tmp, form='(a)'
    file_ke[i] = tmp

end

close, lun
free_lun, lun

return

end
