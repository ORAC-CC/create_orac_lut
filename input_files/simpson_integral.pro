function simpson_integral,x,y
  deltax = x(1:*)-x
  Yave = (y+y(1:*))/2
  return,total(x*y)
end