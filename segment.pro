pro segment, oldx, oldy, x, y, fe, minn =minn
; This routine takes a function descibed by oldx,oldy and replaces it with a similar shaped function x,y whose integral is within
; a fraction fe of the original function.  It does this by fitting line segments to the original function.  Each line segment must
; not intrduce an absolute unceeratiny bigger then 'tol'. As 'tol' is not known to achieve fe it is slowly reduced. 
; using the keyword minn you can force the fit to have a minn number of new points.
; The integrate_trapeziodal test catches  when the old function is made up of two lines 

  If (oldx[0] gt oldx[-1]) then swap = 1 else swap = 0
  
  If (swap) then begin
    oldx=reverse(oldx)
    oldy=reverse(oldy)
  Endif
  
  TrueTotalArea = int_tabulated (oldx,oldy)
  TrueTrapeziodalArea = integrate_trapeziodal(oldx,oldy)
  
  n = n_elements(oldx)-1 ; last indice in input array
   
; initial guess at line segment tolerance
  tol = 0.1

  MaxSegmentError = Tol * TrueTotalArea  

  repeat begin

;   Perform segmentation
;   set first point in new function
    x = [oldx[0]] 
    y = [oldy[0]]

;   old function indices  
    i = 0
    j = 0 ; note is incremented to 1 at start of repeat loop

;   new function indices
    k = 0
    l = 1

    repeat begin
      j = j+1
      areaold = int_tabulated (oldx(i:j),oldy(i:j)) ; area of current function over this segment
;     build trial segment by adding  old value at j onto new function       
      xt = [x,oldx[j]] 
      yt = [y,oldy[j]]
      areanew = int_tabulated (xt(k:l),yt(k:l))     ; area of new function over this segment

 ;    if the new area is too different from the old area  then have added on one point too many
      if abs(areanew - areaold) gt MaxSegmentError then begin
 ;      accept j-1 th point as the new function value 
        x = [x,oldx[j-1]]
        y = [y,oldy[j-1]]
 ;      increment new function indices       
        k = l
        l = l+1
 ;      reset old function indices so that i is the last accepted point       
        i = j - 1
        j = j - 1  ; will be incremented by 1 at the start of loop
      endif
    endrep until j eq n ; stop when j is last point
;   add last point to the new function  
    x = [x,oldx[n]]
    y = [y,oldy[n]]

    TrapeziodalArea = integrate_trapeziodal(x,y)
    If (TrapeziodalArea Eq TrueTrapeziodalArea) then LS = 1 else LS = 0
       
;   evaluate fractional area in new function      
    currentfe = (abs(TrueTotalArea-int_tabulated(x,y))/TrueTotalArea) < (abs(TrapeziodalArea-TrueTrapeziodalArea)/TrueTrapeziodalArea)

;   reduce MaxSegmentError by 10% in case the current fit is not good enough    
    MaxSegmentError = MaxSegmentError*.9   
 
    If keyword_set(minn) then number_OK = n_elements(X) ge minn else number_OK = 1
 ;   print,n_elements(X),(abs(TrueTotalArea-int_tabulated(x,y))/TrueTotalArea), (abs(TrapeziodalArea-TrueTrapeziodalArea)/TrueTrapeziodalArea),currentfe,fe
   
  endrep until (((currentfe  le fe) and number_OK) OR LS); stop when the desired fractional error in area is satisfied.   
  If (swap) then begin
    x=reverse(x)
    y=reverse(y)
  Endif
end