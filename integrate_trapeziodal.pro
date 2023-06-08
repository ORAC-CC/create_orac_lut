FUNCTION integrate_trapeziodal,x,y
; Trapezoidal rule integration
  dx = x(1:*)-x
  ym = 0.5*(y+y(1:*))
  RETURN, TOTAL(dx*ym)
END