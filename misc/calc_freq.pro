;Calculate frequency array
;Taken from https://www.l3harrisgeospatial.com/docs/fft.html
;James Ruppert
;7/1/21
function calc_freq,N,T

  X = FINDGEN((N - 1)/2) + 1

  is_N_even = (N MOD 2) EQ 0

  if (is_N_even) then $
    freq = [0.0, X, N/2, -N/2 + X]/(N*T) $
  else $
    freq = [0.0, X, -(N/2 + 1) + X]/(N*T)

  return,freq

end
