pro idl_total_reference
  ; Small, single-threaded IDL 8.9 reference cases for the V22 TOTAL mirror.
  ; TPOOL_NOTHREAD makes the required accumulation order explicit.
  f = [1.0e8, 1.0, 1.0, -1.0e8, 3.0]
  d = [1.0d16, 1.0d0, 1.0d0, -1.0d16, 3.0d0]
  a = reform([1.0e8, 1.0, 1.0, -1.0e8, 3.0, 4.0], 2, 3)
  i = [1s, 2s, 3s, 4s]

  tf = total(f, /tpool_nothread)
  tfd = total(f, /double, /tpool_nothread)
  td = total(d, /tpool_nothread)
  ti = total(i, /tpool_nothread)
  ta1 = total(a, 1, /tpool_nothread)
  ta2 = total(a, 2, /tpool_nothread)
  cf = total(f, /cumulative, /tpool_nothread)
  ca2 = total(a, 2, /cumulative, /tpool_nothread)

  print, 'float_total_type=', size(tf, /type), ' value=', tf, format='(A,I0,A,G0.9)'
  print, 'float_double_type=', size(tfd, /type), ' value=', tfd, format='(A,I0,A,G0.17)'
  print, 'double_total_type=', size(td, /type), ' value=', td, format='(A,I0,A,G0.17)'
  print, 'int16_total_type=', size(ti, /type), ' value=', ti, format='(A,I0,A,G0.9)'
  print, 'dimension_1_type=', size(ta1, /type), ' values=', ta1, format='(A,I0,A,3(G0.9,1X))'
  print, 'dimension_2_type=', size(ta2, /type), ' values=', ta2, format='(A,I0,A,2(G0.9,1X))'
  print, 'cumulative_type=', size(cf, /type), ' values=', cf, format='(A,I0,A,5(G0.9,1X))'
  print, 'cumulative_dimension_2_type=', size(ca2, /type), ' values=', ca2, format='(A,I0,A,6(G0.9,1X))'
end
