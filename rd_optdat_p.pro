; Reads data in the "optdat" aerosol optical properties file format
; used by OPAC (Optical Properties of Aerosols and Clouds) programme
; and GADS (Global Aerosol Data Set) from Michael Hess et. al.,
; Meteorologisches Institut der Universitaet Muenchen
; History:
; 13/11/XX G Thomas: Original version.

pro rd_optdat_p, fname, distname, rmodd, rmodw, rmin, rmax, sigma, $
                 lambda, Bext, Bsca, Babs, g, w, n, k, theta, phasef

; rmodd  : Mode radius (dry)
; rmodw  : Mode radius (wet - at specified humidity)
; rmin   : Minimum radius of distribution
; rmax   : Maximum radius of distribution
; sigma  : S = ln(std dev of alog(r))
; lambda : Wavelengths
; Bext   : Extinction coefficient
; Bsca   : Scattering coefficient
; Babs   : Absorption coefficient
; g      : Asymmetry parameter
; w      : Single scattering albedo
; n      : Refractive index (real)
; k      : Refractive index (imaginary)
; theta  : Scattering angles for phase function
; phasef : Phase function

openr,lun,fname,/get_lun

line = ''

for i=1,3 do readf,lun,line
words = strsplit(line,/extract)
distname = words(n_elements(words)-1)
for i=1,3 do readf,lun,line
words = strsplit(line,/extract)
Rmin = float(words(n_elements(words)-1))
readf,lun,line
words = strsplit(line,/extract)
Rmax = float(words(n_elements(words)-1))
readf,lun,line
words = strsplit(line,/extract)
sigma = float(words(n_elements(words)-1))
readf,lun,line
words = strsplit(line,/extract)
Rmodw = float(words(n_elements(words)-1))
readf,lun,line
words = strsplit(line,/extract)
Rmodd = float(words(n_elements(words)-1))

for i=1,7 do readf,lun,line

readf,lun,line
words = strsplit(line,/extract)
lambda = float(words(1))
Bext = float(words(2))
Bsca = float(words(3))
Babs = float(words(4))
w = float(words(5))
g = float(words(6))
ext_nor = float(words(7))
m = complex(float(words(8)),float(words(9)))
readf,lun,line
words = strsplit(line,/extract)
while words(0) eq '#' and n_elements(words) eq 10 do begin
    lambda = [lambda, float(words(1))]
    Bext = [Bext, float(words(2))]
    Bsca = [Bsca, float(words(3))]
    Babs = [Babs, float(words(4))]
    w = [w, float(words(5))]
    g = [g, float(words(6))]
    ext_nor = [ext_nor, float(words(7))]
    m = [m, complex(float(words(8)),float(words(9)))]
    readf,lun,line
    words = strsplit(line,/extract)
endwhile

for i=1,5 do readf,lun,line

readf,lun,line
words = strsplit(line,/extract)
lambda1 = float(words(2:*))
readf,lun,line
readf,lun,line
words = strsplit(line,/extract)
theta = words(1)
PhaseF = float(words(2:*))
PhaseF = transpose(PhaseF)
while not eof(lun) do begin
    readf,lun,line
    words = strsplit(line,/extract)
    theta = [theta, words(1)]
    PhaseF = [PhaseF, transpose(float(words(2:*)))]
endwhile

if max(abs(lambda1 - lambda)) gt 0.0 then $
  print,"Warning: Wavelength values don't seem to match!"

n = float(m) & k = imaginary(m)

close, lun
free_lun, lun

end
