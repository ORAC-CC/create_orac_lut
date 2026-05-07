Quadtype=['T','G','R','L']

npts = 100
for i = 0,n_elements(Quadtype)-1 do begin
  Quadrature, Quadtype[i], NPts, Abscissa, Weight
  print, Quadtype[i],total(Weight),total(Weight*abscissa^2)
endfor
end
